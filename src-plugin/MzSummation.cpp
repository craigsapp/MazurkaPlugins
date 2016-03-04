//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Wed Mar  7 23:58:43 PST 2007
// Last Modified: Wed Mar  7 23:58:47 PST 2007
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzSummation.cpp
// URL:           http://sv.mazurka.org.uk/src/MzSummation.cpp
// Documentation: http://sv.mazurka.org.uk/MzSummation
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Calculate the power of an audio signal as it changes 
//                over time.
//

// Defines used in getPluginVersion():
#define P_VER    "200703080"
#define P_NAME   "mzsummation"

#include "MzSummation.h" 
#include "MazurkaWindower.h"

#include <math.h>
#include <stdlib.h>


///////////////////////////////////////////////////////////////////////////
//
// Vamp Interface Functions
//

///////////////////////////////
//
// MzSummation::MzSummation -- class constructor.
//

MzSummation::MzSummation(float samplerate) : MazurkaPlugin(samplerate) {
   // do nothing
}



///////////////////////////////
//
// MzSummation::~MzSummation -- class destructor.
//

MzSummation::~MzSummation() {
   // do nothing
}



////////////////////////////////////////////////////////////
//
// optional polymorphic parameter functions inherited from PluginBase:
//
// Note that the getParameter() and setParameter() polymorphic functions
// are handled in the MazurkaPlugin class.
//

//////////////////////////////
//
// MzSummation::getParameterDescriptors -- return a list of
//      the parameters which can control the plugin.
//

MzSummation::ParameterList MzSummation::getParameterDescriptors(void) const {

   ParameterList       pdlist;
   ParameterDescriptor pd;

   // first parameter: The size of the analysis window in samples
   pd.identifier   = "windowsize";
   pd.name         = "Window size";
   pd.unit         = "samples";
   pd.minValue     = 1.0;
   pd.maxValue     = 100000.0;
   pd.defaultValue = 2048.0;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   // second parameter: The hop size between windows in samples
   pd.identifier   = "hopsize";
   pd.name         = "Window hop size";
   pd.unit         = "samples";
   pd.minValue     = 1.0;
   pd.maxValue     = 10000.0;
   pd.defaultValue = 441.0;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   return pdlist;
}


////////////////////////////////////////////////////////////
//
// optional polymorphic functions inherited from Plugin:
//

/////////////////////////////
//
// MzSummation::getPreferredStepSize --
//

size_t MzSummation::getPreferredStepSize(void) const { 
   return size_t(getParameter("hopsize"));
}



/////////////////////////////
//
// MzSummation::getPreferredBlockSize --
//

size_t MzSummation::getPreferredBlockSize(void) const { 
   return size_t(getParameter("windowsize"));
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzSummation::getIdentifier(void) const
   { return "mzsummation"; }

std::string MzSummation::getName(void) const
   { return "Summation"; }

std::string MzSummation::getDescription(void) const
   { return "Summation"; }

std::string MzSummation::getMaker(void) const
   { return "The Mazurka Project"; }

std::string MzSummation::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

int MzSummation::getPluginVersion(void) const {
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
// MzSummation::getInputDomain -- the host application needs
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

MzSummation::InputDomain MzSummation::getInputDomain(void) const { 
   return TimeDomain; 
}



//////////////////////////////
//
// MzSummation::getOutputDescriptors -- return a list describing
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

MzSummation::OutputList MzSummation::getOutputDescriptors(void) const {

   OutputList       list;
   OutputDescriptor od;

   // First output channel: waveform summation
   od.identifier       = "timedomain";
   od.name             = "Time domain summation";
   od.unit             = "sum";
   od.hasFixedBinCount = true;
   od.binCount         = 1;
   od.hasKnownExtents  = false;
   // od.minValue      = 0.0;
   // od.maxValue      = 0.0;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::VariableSampleRate;
   // od.sampleRate    = 0.0;
   list.push_back(od);

   // Second output channel: frequency summation
   od.identifier       = "freqdomain";
   od.name             = "Frequency domain summation";
   od.unit             = "sum";
   od.hasFixedBinCount = true;
   od.binCount         = 1;
   od.hasKnownExtents  = false;
   // od.minValue      = 0.0;
   // od.maxValue      = 0.0;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::VariableSampleRate;
   // od.sampleRate    = 0.0;
   list.push_back(od);

   return list; 
}



//////////////////////////////
//
// MzSummation::initialise -- this function is called once
//     before the first call to process().
//

bool MzSummation::initialise(size_t channels, size_t stepsize, 
      size_t blocksize) {

   if (channels < getMinChannelCount() || channels > getMaxChannelCount()) {
      return false;
   }

   // step size and block size should never be zero
   if (stepsize <= 0 || blocksize <= 0) {
      return false;
   }

   setChannelCount(channels);
   setStepSize(stepsize);
   setBlockSize(blocksize);

   mz_transformer.setSize(getBlockSize());

   return true;
}



//////////////////////////////
//
// MzSummation::process -- This function is called sequentially on the 
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

MzSummation::FeatureSet MzSummation::process(AUDIODATA inputbufs, 
      Vamp::RealTime timestamp) {

   if (getChannelCount() <= 0) {
      std::cerr << "ERROR: MzSummation::process: "
                << "MzSummation has not been initialized"
                << std::endl;
      return FeatureSet();
   }

   FeatureSet returnFeatures;
   Feature    feature;
   int        i;
   double     sum;

   sum = 0.0;
   for (i=0; i<getBlockSize(); i++) {
      sum += inputbufs[0][i];
   }

   feature.values.push_back(sum);
   feature.hasTimestamp = true;
   feature.timestamp    = timestamp;
   returnFeatures[0].push_back(feature);

   // calculate the spectrum and sum the magnitude
   for (i=0; i<getBlockSize(); i++) {
      mz_transformer.signalNonCausal(i) = double(inputbufs[0][i]);
   }
   mz_transformer.doTransform();
   sum = 0.0;
   for (i=0; i<getBlockSize(); i++) {
      sum += mz_transformer.getSpectrumMagnitude(i);
   }

   feature.values.clear();
   feature.values.push_back(sum);
   returnFeatures[1].push_back(feature);

   return returnFeatures;
}



//////////////////////////////
//
// MzSummation::getRemainingFeatures -- This function is called
//    after the last call to process() on the input data stream has 
//    been completed.  Features which are non-causal can be calculated 
//    at this point.  See the comment above the process() function
//    for the format of output Features.
//

MzSummation::FeatureSet MzSummation::getRemainingFeatures(void) {
   return FeatureSet();
}



//////////////////////////////
//
// MzSummation::reset -- This function may be called after data processing
//    has been started with the process() function.  It will be called when
//    processing has been interrupted for some reason and the processing
//    sequence needs to be restarted (and current analysis output thrown out).  
//    After this function is called, process() will start at the beginning
//    of the input selection as if initialise() had just been called.
//    Note, however, that initialise() will NOT be called before processing 
//    is restarted after a reset().
//

void MzSummation::reset(void) {
   // do nothing
}



///////////////////////////////////////////////////////////////////////////
//
// Non-Interface Functions 
//

