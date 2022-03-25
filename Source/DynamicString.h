/*
  ==============================================================================

    DynamicString.h
    Created: 25 Mar 2022 11:33:06am
    Author:  Helmer Nuijens

  ==============================================================================
*/

#pragma once

#define _USE_MATH_DEFINES
#include <JuceHeader.h>
#include <cmath>
#include <vector>
using namespace std;

class DynamicString
{
public:

    DynamicString(); // Constructor

    ~DynamicString();                 // Destructor

    // Public Methods:
    void setFs(double Fs);
    void setGridSize(double frequency);  // set grid size N
    void exciteSystem(double width, double excitationLoc);
    void processSample();
    double getOutput(double ratio);              //
    void setDamping(double sigma0, double sigma1);
    vector<double> getStringState(int maxLength);

    // Damping
    double sigma0 = 0.1; // Frequency independent damping
    double sigma1 = 0.01; // Frequency dependent damping

private:

    // Private Methods:
    void addPoint();
    void removePoint();
    void getVirtualGridPoints();
    void getSchemeWeights();
    void getConnectionForce();
    void calculateScheme();
    void updateStates();
    void resetGrid();

    // System defined:
    double Fs = 48000.f;         // Sample rate 
    double k;                    // Descrete time grid spacing

    // User defined:
    double L = 1.;               // String length
    double c;                    // Dynamic wave speed
    double h;                    // Grid spacing
    double N = 0;                // Current grid size
    double N1 = 0;               // Previous grid size
    double alpha;                // Distance between the systems

    int M;                       // Length of left system
    int Mw;                      // lengtg of right system

    // State vectors:
    vector<vector<double>> uStates;      // State vectors of left system
    vector<double*> u;                   // Pointer to the left state vectors

    vector<vector<double>> wStates;      // State vectors of right system
    vector<double*> w;                   // Pointer to the right state vectors 

    //vector<double> U;                    // Full string
    //vector<double> newU;                 // Resized string
    // Interpolation vectors:
    vector<double> I2 = { 0,0,0 };       // Kwadratic interpolation vector 
    vector<double> I3 = { 0,0,0,0 };     // Cubic interpolation vector
    vector<double> I2Flipped = I2;       // Kwadratic interpolation vector flipped 
    vector<double> I3Flipped = I3;       // Cubic interpolation vector flipped

    // Virtual grid points:
    double uM1 = 0.f;                  // u(M+1)
    double wMin1 = 0.f;                // w(-1)
    double uPrevM1 = 0.f;              // u(M+1)^{n-1}
    double wPrevMin1 = 0.f;            // w(-1)^{n-1}

    // Connection force:
    double r1;                         // Coefficient 1
    double r2;                         // Coefficient 2
    double uMI;                        // Intermediate scheme u
    double w0I;                        // Intermediate scheme w
    double Fc;                         // Connection force 
    double omegaS = 1.;                // Frequency connection spring
    double sigmaS = 1.;                // Damping connection spring

    // Scheme weights:
    double D;
    double C1;
    double C2;
    double C3;
    double C4;
    double C5;

    // Output
    double out;


    int index = 0;  // Keeping track of pointer switching
    int steps;      // steps for resizing U
};

