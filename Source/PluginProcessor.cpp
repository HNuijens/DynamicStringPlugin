/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
DynamicStringPluginAudioProcessor::DynamicStringPluginAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
    addParameter(fundFreq = new AudioParameterFloat("fundamentalFreq", // parameter ID
        "Fundamental Frequency", // parameter name
        20.0f,          // minimum value
        10000.0f,       // maximum value
        220.0f));          // default value

    addParameter(modulation = new AudioParameterFloat("modulation", // parameter ID
        "modulation", // parameter name
        -12.0f,          // minimum value
        12.0f,       // maximum value
        0.0f));          // default value


    addParameter(excited = new AudioParameterBool("excited", // parameter ID
        "excited", // parameter name
        false   // default value
    )); // default value
}

DynamicStringPluginAudioProcessor::~DynamicStringPluginAudioProcessor()
{
}

//==============================================================================
const juce::String DynamicStringPluginAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool DynamicStringPluginAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool DynamicStringPluginAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool DynamicStringPluginAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double DynamicStringPluginAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int DynamicStringPluginAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int DynamicStringPluginAudioProcessor::getCurrentProgram()
{
    return 0;
}

void DynamicStringPluginAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String DynamicStringPluginAudioProcessor::getProgramName (int index)
{
    return {};
}

void DynamicStringPluginAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
}

//==============================================================================
void DynamicStringPluginAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    
    dynamicString.setFs(sampleRate);
    
    rootNote = *fundFreq;
    mod = *modulation; 
    f0 = rootNote * powf(2.f, mod / 12.0);

    parameters.set("L", 1.0);
    parameters.set("f0", f0);
    parameters.set("sig0", sig0);
    parameters.set("sig1", sig1);

    dynamicString.setGrid(parameters);

}

void DynamicStringPluginAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool DynamicStringPluginAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void DynamicStringPluginAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();


    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

 
    if (*excited)
    {
        dynamicString.exciteSystem(15.0f, 0.3f);
        *excited = false;
    }

    if (mod != *modulation || *fundFreq != rootNote)
    {
        rootNote = *fundFreq;
        mod = *modulation;
        f0 = rootNote * powf(2.f, mod / 12.0);
        parameters.set("f0", f0);
    }
   


    for (int n = 0; n < buffer.getNumSamples(); ++n)
    {
        auto outL = buffer.getWritePointer(0);
        auto outR = buffer.getWritePointer(1);

        float out = 0.f;

        dynamicString.setGrid(parameters);
        out = dynamicString.getNextSample(0.2f);

        out = limit(out, -1.0f, 1.0f);
        outL[n] = out;
        outR[n] = out;
    }
    
}

//==============================================================================
bool DynamicStringPluginAudioProcessor::hasEditor() const
{
    return false; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* DynamicStringPluginAudioProcessor::createEditor()
{
    return new DynamicStringPluginAudioProcessorEditor (*this);
}

//==============================================================================
void DynamicStringPluginAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void DynamicStringPluginAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new DynamicStringPluginAudioProcessor();
}

double DynamicStringPluginAudioProcessor::limit(double sample, double min, double max)
{
    if (sample < min) return min;
    else if (sample > max) return max;
    else return sample; 
}