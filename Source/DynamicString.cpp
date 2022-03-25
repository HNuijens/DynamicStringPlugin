/*
  ==============================================================================

    DynamicString.cpp
    Created: 25 Mar 2022 11:33:06am
    Author:  Helmer Nuijens

  ==============================================================================
*/
#include "DynamicString.h"

DynamicString::DynamicString()
{

}

DynamicString::~DynamicString()
{

}

void DynamicString::setFs(double Fs)
{
    this->Fs = Fs;
    k = 1. / Fs;
}

void DynamicString::setGrid(NamedValueSet& parameters)
{
    // Fetch parameters
    f0 = *parameters.getVarPointer("f0");
    L = *parameters.getVarPointer("L");
    sig0 = *parameters.getVarPointer("sig0");
    sig1 = *parameters.getVarPointer("sig1");

    c = 2.0f * L * f0;
    h = sqrt(c * c * k * k + 4.0f * sig1 * k);
    N = L / h;

    if (floor(N) - floor(N1) > 1 || floor(N) - floor(N1) < -1)
    {
        resetGrid();
    }

    getSchemeWeights();
}
void DynamicString::resetGrid()
{
    M = ceil(0.5 * L / h);
    Mw = floor(0.5 * L / h);

    uStates = vector<vector<double>>(3, vector<double>(M, 0));
    wStates = vector<vector<double>>(3, vector<double>(Mw, 0));

    u.resize(3, nullptr);
    w.resize(3, nullptr);

    for (int n = 0; n < 3; n++)
    {
        u[n] = &uStates[n][0];
        w[n] = &wStates[n][0];
    }

    N1 = N;

    wPrevMin1 = wMin1 = 0;
    uPrevM1 = uM1 = 0;
}

void DynamicString::exciteSystem(double width, double excitationLoc)
{
    //// Raised cosine excitation ////
    // width (in grid points) of the excitation
    // make sure we’re not going out of bounds at the left boundary
    int switchVal = 0;
    int start = max(floor((N + 1) * excitationLoc) - floor(width * 0.5), 1.0);
    int startw = 0;
    if (start > M)
    {
        startw = start - M;
    }

    // DBG(start);
    for (int l = 0; l < width; ++l)
    {
        // make sure we’re not going out of bounds
        // at the right boundary (this does ’cut off’
        // the raised cosine)
        if (l + startw - switchVal > Mw - 1)
            break;
        else if (l + start < M)
        {
            u[1][l + start] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width - 1.0)));
            u[2][l + start] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width - 1.0)));
        }
        else if (l + start == M)
        {
            switchVal = l + 1;
        }
        else if (l + start >= M)
        {
            w[1][l + startw - switchVal] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width - 1.0)));
            w[2][l + startw - switchVal] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width - 1.0)));
        }
    }
}

double DynamicString::getNextSample(float outputPos)
{
    alpha = N - floor(N); // calc distance between two systems

    // Check if there is a need to add/delete grid points
    if (floor(N) > floor(N1)) addPoint();
    else if (floor(N) < floor(N1)) removePoint();

    getVirtualGridPoints();
    //getSchemeWeights();
    getConnectionForce();
    calculateScheme();
    
    double out;

    // Fetch output from right system
    int pos = static_cast<int>(round(outputPos * N));
    if (pos >= M)
    {
        pos = pos - M;
        out = w[1][pos];
    }

    else out = u[1][pos];

    updateStates();

    return out; 
}

void DynamicString::addPoint()
{
    // Calculate cubic interpolation
    I3[0] = -(alpha * (alpha + 1)) / ((alpha + 2) * (alpha + 3));
    I3[1] = 2 * alpha / (alpha + 2);
    I3[2] = 2 / (alpha + 2);
    I3[3] = -(2 * alpha) / ((alpha + 3) * (alpha + 2));
    I3Flipped = I3;
    reverse(I3Flipped.begin(), I3Flipped.end()); // Flip vector

    // Floored N is even, add to right system
    if (static_cast<int>(floor(N)) % 2 == 0)
    {
        // Add to the start of w
        wStates[0].insert(wStates[0].begin(), 0);
        wStates[1].insert(wStates[1].begin(), 0);
        wStates[2].insert(wStates[2].begin(), 0);

        for (int n = 0; n < 3; n++)
        {
            w[n] = &wStates[(((index + n) % 3) + 3) % 3][0];
        }

        w[1][0] = I3Flipped[0] * u[1][M - 2] + I3Flipped[1] * u[1][M - 1] + I3Flipped[2] * w[1][1] + I3Flipped[3] * w[1][2];
        w[2][0] = I3Flipped[0] * u[2][M - 2] + I3Flipped[1] * u[2][M - 1] + I3Flipped[2] * w[2][1] + I3Flipped[3] * w[2][2];

        Mw++;
    }
    else // Floored N is odd, add the to left system
    {
        // Add to the end of u
        uStates[0].push_back(0);
        uStates[1].push_back(0);
        uStates[2].push_back(0);

        for (int n = 0; n < 3; n++)
        {
            u[n] = &uStates[(((index + n) % 3) + 3) % 3][0];
        }

        u[1][M] = I3[0] * u[1][M - 2] + I3[1] * u[1][M - 1] + I3[2] * w[1][0] + I3[3] * w[1][1];
        u[2][M] = I3[0] * u[2][M - 2] + I3[1] * u[2][M - 1] + I3[2] * w[2][0] + I3[3] * w[2][1];

        M++;
    }
}

void DynamicString::removePoint()
{
    // Floored N is even, remove from left system
    if (static_cast<int>(floor(N)) % 2 == 0)
    {
        // Remove from right side of u
        uStates[0].pop_back();
        uStates[1].pop_back();
        uStates[2].pop_back();

        M--;
    }
    else // Floored N is odd, remove from left system
    {
        // Remove from left side of w
        wStates[0].erase(wStates[0].begin());
        wStates[1].erase(wStates[1].begin());
        wStates[2].erase(wStates[2].begin());

        Mw--;
    }
}

void DynamicString::getVirtualGridPoints()
{
    I2[0] = -(alpha - 1) / (alpha + 1);
    I2[1] = 1;
    I2[2] = (alpha - 1) / (alpha + 1);
    I2Flipped = I2;
    reverse(I2Flipped.begin(), I2Flipped.end()); // Flip vector

    uM1 = u[1][M - 1] * I2Flipped[0] + w[1][0] * I2Flipped[1] + w[1][1] * I2Flipped[2];
    wMin1 = u[1][M - 2] * I2[0] + u[1][M - 1] * I2[1] + w[1][0] * I2[2];
}

void DynamicString::getSchemeWeights()
{
    D = 1. + sig0 * k;
    C1 = 2.;
    C2 = sig0 * k - 1.;
    C3 = (c * c * k * k) / (h * h);
    C4 = (2 * sig1 * k) / (h * h);
    C5 = (k * k / h);

    C1 = C1 / D;
    C2 = C2 / D;
    C3 = C3 / D;
    C4 = C4 / D;
    C5 = C5 / D;
}

void DynamicString::getConnectionForce()
{
    r1 = (omegaS * omegaS - (sigmaS / k)) / (omegaS * omegaS + (sigmaS / k));
    r2 = (h * (1 + sig0 * k) * (1 - alpha) * (omegaS * omegaS + (sigmaS / k))) / (2 * h * (1 + sig0 * k) * alpha + 2 * k * k * (1 - alpha) * (omegaS * omegaS + (sigmaS / k)));
    uMI = C1 * u[1][M - 1] + C2 * u[2][M - 1] + C3 * (uM1 - 2 * u[1][M - 1] + u[1][M - 2]) + C4 * (uM1 - 2 * u[1][M - 1] + u[1][M - 2] - uPrevM1 + 2 * u[2][M - 1] - u[2][M - 2]);
    w0I = C1 * w[1][0] + C2 * w[2][0] + C3 * (w[1][1] - 2 * w[1][0] + wMin1) + C4 * (w[1][1] - 2 * w[1][0] + wMin1 - w[2][1] + 2 * w[2][0] - wPrevMin1);
    Fc = r2 * (w0I - uMI + r1 * (w[2][0] - u[2][M - 1]));
}

void DynamicString::calculateScheme()
{
    for (int m = 1; m < M; m++)
    {
        if (m == M - 1)
        {
            u[0][m] = C1 * u[1][m] + C2 * u[2][m] + C3 * (uM1 - 2 * u[1][m] + u[1][m - 1]) + C4 * (uM1 - 2 * u[1][m] + u[1][m - 1] - uPrevM1 + 2 * u[2][m] - u[2][m - 1]) + C5 * Fc;
        }
        else
        {
            u[0][m] = C1 * u[1][m] + C2 * u[2][m] + C3 * (u[1][m + 1] - 2 * u[1][m] + u[1][m - 1]) + C4 * (u[1][m + 1] - 2 * u[1][m] + u[1][m - 1] - u[2][m + 1] + 2 * u[2][m] - u[2][m - 1]);
        }
    }

    for (int m = 0; m < Mw - 1; m++)
    {
        if (m == 0)
        {
            w[0][m] = C1 * w[1][m] + C2 * w[2][m] + C3 * (w[1][m + 1] - 2 * w[1][m] + wMin1) + C4 * (w[1][m + 1] - 2 * w[1][m] + wMin1 - w[2][m + 1] + 2 * w[2][m] - wPrevMin1) - C5 * Fc;
        }
        else
        {
            w[0][m] = C1 * w[1][m] + C2 * w[2][m] + C3 * (w[1][m + 1] - 2 * w[1][m] + w[1][m - 1]) + C4 * (w[1][m + 1] - 2 * w[1][m] + w[1][m - 1] - w[2][m + 1] + 2 * w[2][m] - w[2][m - 1]);
        }
    }
}

void DynamicString::updateStates()
{
    if (index < 3) index++;
    else index = 0;

    double* uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;

    double* wTmp = w[2];
    w[2] = w[1];
    w[1] = w[0];
    w[0] = wTmp;

    N1 = N;

    wPrevMin1 = wMin1;
    uPrevM1 = uM1;
}