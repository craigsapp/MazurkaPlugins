//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Thu Jan  4 09:38:06 PST 2007
// Last Modified: Sat Jan  6 03:56:55 PST 2007
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzPitchPower.cpp
// URL:           http://sv.mazurka.org.uk/src/MzPitchPower.cpp
// Documentation: http://sv.mazurka.org.uk/MzPitchPower
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Measure the power of a particular note by combining
//                the individual power in each harmonic.
//

// Defines used in getPluginVersion():
#define P_VER    "200701060"
#define P_NAME   "MzPitchPower"

#include "MzPitchPower.h" 

#include <stdio.h>
#include <math.h>

#include <string>

#define METHOD_ADDITION 0
#define METHOD_MULTIPLICATION 1

using namespace std;        // avoid stupid std:: prefixing


///////////////////////////////////////////////////////////////////////////
//
// Vamp Interface Functions
//

///////////////////////////////
//
// MzPitchPower::MzPitchPower -- class constructor.  The values
//   for the mz_* variables are just place holders demonstrating the
//   default value.  These variables will be set in the initialise()
//   function from the user interface.
//

MzPitchPower::MzPitchPower(float samplerate) : MazurkaPlugin(samplerate) {
   mz_harmonics     = 5;            // how many harmonics to consider
   mz_method        = METHOD_ADDITION;  // add harmonic powers together
   mz_transformsize = 16384;     // size of the FFT transformations
}



///////////////////////////////
//
// MzPitchPower::~MzPitchPower -- class destructor.
//

MzPitchPower::~MzPitchPower() {
   // do nothing
}


////////////////////////////////////////////////////////////
//
// parameter functions --
//

//////////////////////////////
//
// MzPitchPower::getParameterDescriptors -- return a list of
//      the parameters which can control the plugin.
//

MzPitchPower::ParameterList 
MzPitchPower::getParameterDescriptors(void) const {

   ParameterList       pdlist;
   ParameterDescriptor pd;

   // first parameter: Number of samples in the audio window
   pd.identifier   = "windowsamples";
   pd.name         = "Window Size";
   pd.unit         = "samples";
   pd.minValue     = 2.0;
   pd.maxValue     = 10000;
   pd.defaultValue = 2048.0;
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // second parameter: Step size between analysis windows
   pd.identifier   = "stepsamples";
   pd.name         = "Step Size";
   pd.unit         = "samples";
   pd.minValue     = 2.0;
   pd.maxValue     = 30000.0;
   pd.defaultValue = 441.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // third parameter: Number of harmonics to consider
   pd.identifier   = "harmonics";
   pd.name         = "Harmonics";
   pd.unit         = "";
   pd.minValue     = 0.0;
   pd.maxValue     = 100.0;
   pd.defaultValue = 5.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // fourth parameter: Search Pitch
   pd.identifier   = "pitch";
   pd.name         = "Pitch";
   pd.unit         = "";
   pd.minValue     = 0.0;
   pd.maxValue     = 127.0;
   generateMidiNoteList(pd.valueNames, 0, 127);
   pd.defaultValue = 67.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // fifth parameter: Cents deviation
   pd.identifier   = "cents";
   pd.name         = "Cents";
   pd.unit         = "cent";
   pd.minValue     = -100.0;
   pd.maxValue     =  100.0;
   pd.defaultValue = 0.0; 
   pd.isQuantized  = false;
   // pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // sixth parameter: Tuning
   pd.identifier   = "tuning";
   pd.name         = "A4 Tuning";
   pd.unit         = "Hz";
   pd.minValue     = 20.0;
   pd.maxValue     = 1000.0;
   pd.defaultValue = 440.0; 
   pd.isQuantized  = false;
   // pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // seventh parameter: Fundamental frequency
   pd.identifier   = "freq";
   pd.name         = "or\nFundamental\nFrequency";
   pd.unit         = "Hz";
   pd.minValue     = -1.0;
   pd.maxValue     = 10000.0;
   pd.defaultValue = -1.0; 
   pd.isQuantized  = false;
   // pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // eighth parameter: Mixing method
   pd.identifier   = "method";
   pd.name         = "Method";
   pd.unit         = "";
   pd.minValue     = 0.0;
   pd.maxValue     = 1.0;
   pd.defaultValue = 0.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pd.valueNames.push_back("Addition");
   pd.valueNames.push_back("Multiplication");
   pdlist.push_back(pd);
   pd.valueNames.clear();

   return pdlist;
}


////////////////////////////////////////////////////////////
//
// optional polymorphic functions inherited from PluginBase:
//

//////////////////////////////
//
// MzPitchPower::getPreferredStepSize -- overrides the 
//     default value of 0 (no preference) returned in the 
//     inherited plugin class.
//

size_t MzPitchPower::getPreferredStepSize(void) const {
   return getParameterInt("stepsamples");
}



//////////////////////////////
//
// MzPitchPower::getPreferredBlockSize -- overrides the 
//     default value of 0 (no preference) returned in the 
//     inherited plugin class.
//

size_t MzPitchPower::getPreferredBlockSize(void) const {
   return getParameterInt("windowsamples");
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzPitchPower::getIdentifier(void) const
{ return "mzpitchpower"; }

std::string MzPitchPower::getName(void) const
{ return "Pitch Power"; }

std::string MzPitchPower::getDescription(void) const
{ return "Pitch Power"; }

std::string MzPitchPower::getMaker(void) const
{ return "The Mazurka Project"; }

std::string MzPitchPower::getCopyright(void) const
{ return "2006 Craig Stuart Sapp"; }

int MzPitchPower::getPluginVersion(void) const {
   const char *v = "@@VampPluginID@" P_NAME "@" P_VER "@" __DATE__ "@@";
   if (v[0] != '@') { std::cerr << v << std::endl; return 0; }
   return atol(P_VER);
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from Plugin:
//

//////////////////////////////
//
// MzPitchPower::getInputDomain -- the host application needs
//    to know if it should send either:
//
// TimeDomain      == Time samples from the audio waveform.
// FrequencyDomain == Spectral frequency frames which will arrive
//                    in an array of interleaved real, imaginary
//                    values for the complex spectrum (both positive 
//                    and negative frequencies). Zero Hz being the
//                    first frequency sample and negative frequencies
//                    at the far end of the array as is usually done.
//                    Note that frequency data is transmitted from
//                    the host application as floats.  The data will
//                    be transmitted via the process() function which
//                    is defined further below.
//

MzPitchPower::InputDomain MzPitchPower::getInputDomain(void) const { 
   return TimeDomain; 
}



//////////////////////////////
//
// MzPitchPower::getOutputDescriptors -- return a list describing
//    each of the available outputs for the object.  OutputList
//    is defined in the file vamp-sdk/Plugin.h:
//
// .identifier       == short name of output for computer use.  Must not
//                      contain spaces or punctuation.
// .name             == long name of output for human use.
// .unit             == the units or basic meaning of the data in the 
//                      specified output.
// .hasFixedBinCount == true if each output feature (sample) has the 
//                      same dimension.
// .binCount         == when hasFixedBinCount is true, then this is the 
//                      number of values in each output feature.  
//                      binCount=0 if timestamps are the only features,
//                      and they have no labels.
// .binNames         == optional description of each bin in a feature.
// .hasKnownExtent   == true if there is a fixed minimum and maximum
//                      value for the range of the output.
// .minValue         == range minimum if hasKnownExtent is true.
// .maxValue         == range maximum if hasKnownExtent is true.
// .isQuantized      == true if the data values are quantized.  Ignored
//                      if binCount is set to zero.
// .quantizeStep     == if isQuantized, then the size of the quantization,
//                      such as 1.0 for integers.
// .sampleType       == Enumeration with three possibilities:
//   OD::OneSamplePerStep   -- output feature will be aligned with
//                             the beginning time of the input block data.
//   OD::FixedSampleRate    -- results are evenly spaced according to 
//                             .sampleRate (see below).
//   OD::VariableSampleRate -- output features have individual timestamps.
// .sampleRate       == samples per second spacing of output features when
//                      sampleType is set toFixedSampleRate.
//                      Ignored if sampleType is set to OneSamplePerStep
//                      since the start time of the input block will be used.
//                      Usually set the sampleRate to 0.0 if VariableSampleRate
//                      is used; otherwise, see vamp-sdk/Plugin.h for what
//                      positive sampleRates would mean.
//

MzPitchPower::OutputList 
MzPitchPower::getOutputDescriptors(void) const {

   OutputList       odlist;
   OutputDescriptor od;

   // First output channel: Pitch Power
   od.identifier       = "pitchpower";
   od.name             = "Pitch Power";
   od.unit             = "db";
   od.hasFixedBinCount = true;
   od.binCount         = 1;
   od.hasKnownExtents  = false;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::OneSamplePerStep;
   // od.sampleRate    = 0.0;
   odlist.push_back(od);
   #define OUTPUT_PITCH_POWER 0
   od.binNames.clear();

   // second output is an array of the loudness of each harmonic
   // which are summed together in the first output channel.

   // Second output channel: Power of individual harmonics
   od.identifier       = "harmonicpowers";
   od.name             = "Harmonic Powers";
   od.unit             = "db";
   od.hasFixedBinCount = true;
   od.binCount         = mz_harmonics;
   od.hasKnownExtents  = false;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::OneSamplePerStep;
   // od.sampleRate    = 0.0;
   odlist.push_back(od);
   #define OUTPUT_HARMONIC_POWER 1
   od.binNames.clear();

   return odlist; 
}



//////////////////////////////
//
// MzPitchPower::initialise -- this function is called once
//     before the first call to process().
//

bool MzPitchPower::initialise(size_t channels, size_t stepsize, 
      size_t blocksize) {

   if (channels < getMinChannelCount() || channels > getMaxChannelCount()) {
      return false;
   }

   // step size and block size should never be zero
   if (stepsize <= 0 || blocksize <= 0) {
      return false;
   }

   setStepSize(stepsize);
   setBlockSize(blocksize);
   setChannelCount(channels);

   mz_harmonics  = getParameterInt("harmonics");
   mz_method     = getParameterInt("method");

   double midi   = getParameterDouble("pitch");
   double cents  = getParameterDouble("cents");
   double a4tune = getParameterDouble("tune");
   double freq   = getParameterDouble("freq");
   double a4midi = 6.0;

   if (freq < 0.0) {
      freq = a4tune * pow(2.0, (a4midi - midi + cents/100.0) / 12.0);
      cerr << "Pitch Fundamental Frequency: " << freq << endl;
   }

   if (getBlockSize() > mz_transformsize) {
      mz_transformsize = getBlockSize();
   }

   double basebin = freq * mz_transformsize / getSrate();
   mz_bins.resize(mz_harmonics);
   int i;
   cerr << "Transform size: " << mz_transformsize << endl;
   cerr << "Bins:\t";
   for (i=0; i<mz_harmonics; i++) {
      mz_bins[i] = (int)(basebin * (i+1) + 0.5);
      cerr << mz_bins[i];
      if (i < mz_harmonics-1) {
         cerr << "\t";
      }
      cerr << endl;
   }
   
   mz_transformer.setSize(mz_transformsize);
   mz_transformer.zeroSignal();
   mz_windower.setSize(getBlockSize());
   // Note that Hann window or other high side-lobe windows
   // will not work so well because the side-lobes are noise
   // which will over-estimate the energy in harmonics.  However,
   // windows with wide main-lobes (Such as the 4-term Blackman-Harris window)
   // have a similar but less serious problem.
   //
   // One idea to add is to have the size of the window be equal to the
   // size of the fundamental (or a sub harmonic).
   mz_windower.makeWindow("BlackmanHarris4_92");

   return true;
}



//////////////////////////////
//
// MzPitchPower::process -- This function is called sequentially on the 
//    input data, block by block.  After the sequence of blocks has been
//    processed with process(), the function getRemainingFeatures() will 
//    be called.
//
// Here is a reference chart for the Feature struct:
//
// .hasTimestamp   == If the OutputDescriptor.sampleType is set to
//                    VariableSampleRate, then this should be "true".
// .timestamp      == The time at which the feature occurs in the time stream.
// .values         == The float values for the feature.  Should match
//                    OD::binCount.
// .label          == Text associated with the feature (for time instants).
//

MzPitchPower::FeatureSet MzPitchPower::process(AUDIODATA inputbufs, 
      Vamp::RealTime timestamp) {

   if (getStepSize() <= 0) {
      cerr << "ERROR: MzPitchPower::process: "
                << "MzPitchPower has not been initialized" << std::endl;
      return FeatureSet();
   }

   int i;
   Feature    feature;
   FeatureSet returnFeatures;

   mz_windower.windowNonCausal(mz_transformer, inputbufs[0], getBlockSize());
   mz_transformer.doTransform();

   vector<double> harmonic_power;
   extractHarmonicPowers(harmonic_power, mz_bins, mz_transformer);

   double ppower = 0.0;

   switch (mz_method) {


   case METHOD_MULTIPLICATION:
   {
      int counter = 0;
      ppower = 1.0;
      for (i=0; i<(int)harmonic_power.size(); i++) {
         if (harmonic_power[i] != 0.0) {
            ppower *= harmonic_power[i];
            counter++;
         }
      }
      if (counter >= 2) {
         ppower = pow(ppower, 1.0/counter);
      }
   }
   break;


   case METHOD_ADDITION:
   default:
      for (i=0; i<(int)harmonic_power.size(); i++) {
         ppower += harmonic_power[i];
      }
   }

   // convert value to decibels
   if (ppower > 0) {
      ppower = 20.0 * log10(ppower);
   } else {
      ppower = -120.0;
   }

   ////////////////////////////////////////////////////////////////////////
   ///// store the plugin's FIRST output: pitch power /////////////////////
   ////////////////////////////////////////////////////////////////////////

   feature.values.clear();
   feature.values.push_back(ppower);
   feature.hasTimestamp = false;
   returnFeatures[OUTPUT_PITCH_POWER].push_back(feature);


   // convert the harmonic_power values to decibels:

   for (i=0; i<(int)harmonic_power.size(); i++) {
      if (harmonic_power[i] > 0.0) {
         harmonic_power[i] = 20.0 * log10(harmonic_power[i]);
      } else {
         harmonic_power[i] = -120.0;
      }
   }

   ////////////////////////////////////////////////////////////////////////
   ///// store the plugin's SECOND output: individual harmonic powers /////
   ////////////////////////////////////////////////////////////////////////

   feature.values.resize(harmonic_power.size());
   for (i=0; i<(int)harmonic_power.size(); i++) {
      feature.values[i] = harmonic_power[i];
   }
   feature.hasTimestamp = false;
   returnFeatures[OUTPUT_HARMONIC_POWER].push_back(feature);

   return returnFeatures;
}



//////////////////////////////
//
// MzPitchPower::getRemainingFeatures -- This function is called
//    after the last call to process() on the input data stream has 
//    been completed.  Features which are non-causal can be calculated 
//    at this point.  See the comment above the process() function
//    for the format of output Features.
//

MzPitchPower::FeatureSet MzPitchPower::getRemainingFeatures(void) {
   return FeatureSet();
}



//////////////////////////////
//
// MzPitchPower::reset -- This function may be called after data 
//    processing has been started with the process() function.  It will 
//    be called when processing has been interrupted for some reason and 
//    the processing sequence needs to be restarted (and current analysis 
//    output thrown out).  After this function is called, process() will 
//    start at the beginning of the input selection as if initialise() 
//    had just been called.  Note, however, that initialise() will NOT 
//    be called before processing is restarted after a reset().
//

void MzPitchPower::reset(void) {
   // do nothing
}


///////////////////////////////////////////////////////////////////////////
//
// Non-Interface Functions 
//

//////////////////////////////
//
// MzPitchPower::generateMidiNoteList -- Create a list of pitch names 
//   for the specified MIDI key number range.
//

void MzPitchPower::generateMidiNoteList(vector<std::string>& alist,
	int minval, int maxval) {

   alist.clear();
   if (maxval < minval) {
      std::swap(maxval, minval);
   }

   int i;
   int octave;
   int pc;
   char buffer[32] = {0};
   for (i=minval; i<=maxval; i++) {
      octave = i / 12;
      pc = i - octave * 12;
      octave = octave - 1;  // Make middle C (60) = C4
      switch (pc) {
         case 0:   sprintf(buffer, "C%d",  octave); break;
         case 1:   sprintf(buffer, "C#%d", octave); break;
         case 2:   sprintf(buffer, "D%d",  octave); break;
         case 3:   sprintf(buffer, "D#%d", octave); break;
         case 4:   sprintf(buffer, "E%d",  octave); break;
         case 5:   sprintf(buffer, "F%d",  octave); break;
         case 6:   sprintf(buffer, "F#%d", octave); break;
         case 7:   sprintf(buffer, "G%d",  octave); break;
         case 8:   sprintf(buffer, "G#%d", octave); break;
         case 9:   sprintf(buffer, "A%d",  octave); break;
         case 10:  sprintf(buffer, "A#%d", octave); break;
         case 11:  sprintf(buffer, "B%d",  octave); break;
         default:  sprintf(buffer, "x%d", i);
      }
      alist.push_back(buffer);
   }
}



//////////////////////////////
//
// MzPitchPower::extractHarmonicPower -- 
//

void MzPitchPower::extractHarmonicPowers(vector<double>& harmonic_power, 
      vector<int>& bins, MazurkaTransformer& transformer) { 
   harmonic_power.resize(bins.size());

   for (int i=0; i<(int)bins.size(); i++) {
      if (bins[i] >= 0) {
         harmonic_power[i] = transformer.getSpectrumMagnitude(bins[i]);
      } else {
         harmonic_power[i] = 0.0;
      }
   }
}



