//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Tue May  9 05:25:27 PDT 2006
// Last Modified: Sat May 20 05:41:31 PDT 2006 (added parameters)
// Last Modified: Mon May 22 23:54:11 PDT 2006
// Filename:      MzChronogram.cpp
// URL:           http://sv.mazurka.org.uk/src/MzChronogram.cpp
// Documentation: http://sv.mazurka.org.uk/MzChronogram
// Syntax:        ANSI99 C++; vamp plugin
//
// Description:   Display audio signal in two dimensions.
//

#include "MzChronogram.h" 

#include <math.h>

#include <iostream>


///////////////////////////////////////////////////////////////////////////
//
// Vamp Interface Functions
//

///////////////////////////////
//
// MzChronogram::MzChronogram -- class constructor.
//

MzChronogram::MzChronogram(float samplerate) : MazurkaPlugin(samplerate) {
   mz_blocksize = 0;
   mz_stepsize  = 0;
   mz_channels  = 0;
}



///////////////////////////////
//
// MzChronogram::~MzChronogram -- class destructor.
//

MzChronogram::~MzChronogram() {
   // do nothing
}


////////////////////////////////////////////////////////////
//
// Parameter functions --
//

//////////////////////////////
//
// MzChronogram::getParameterDescriptors -- return a list of
//      the parameters which can control the plugin.
//

MzChronogram::ParameterList MzChronogram::getParameterDescriptors(void) const {

   ParameterList       pdlist;
   ParameterDescriptor pd;

   // first parameter: The number of samples on the vertical axis
   pd.name         = "verticalperiod";
   pd.description  = "Vertical period";
   pd.unit         = "samples";
   pd.minValue     = 1.0;
   pd.maxValue     = 10000;
   pd.defaultValue = 551.0;     // close to A440 / 5
   pd.isQuantized  = 1;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);

   // second parameter: The Frequency of the period on the vertical axis
   pd.name         = "verticalfrequency";
   pd.description  = "Vertical frequency";
   pd.unit         = "Hz";
   pd.minValue     = 1.0;
   pd.maxValue     = 10000.0;
   pd.defaultValue = 440.0;   
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   return pdlist;
}


////////////////////////////////////////////////////////////
//
// optional polymorphic functions inherited from PluginBase:
//

//////////////////////////////
//
// MzChronogram::getPreferredStepSize -- overrides the 
//     default value of 0 (no preference) returned in the 
//     inherited plugin class.
//

size_t MzChronogram::getPreferredStepSize(void) const {
   return getPreferredBlockSize();
}



//////////////////////////////
//
// MzChronogram::getPreferredBlockSize -- overrides the 
//     default value of 0 (no preference) returned in the 
//     inherited plugin class.
//

size_t MzChronogram::getPreferredBlockSize(void) const {
   double output = 0.0;
   if (isParameterAtDefault("verticalperiod")) {
      output = getSrate() / getParameter("verticalfrequency");
   } else {
      output = getParameter("verticalperiod");
   }
   return int(output + 0.5);
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzChronogram::getName(void) const
   { return "mzchronogram"; }

std::string MzChronogram::getMaker(void) const
   { return "Craig Stuart Sapp <craig@ccrma.stanford.edu>"; }

std::string MzChronogram::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

std::string MzChronogram::getDescription(void) const
   { return "Chronogram"; }

int MzChronogram::getPluginVersion(void) const
   { return 200605221; }


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from Plugin:
//

//////////////////////////////
//
// MzChronogram::getInputDomain -- the host application needs
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

MzChronogram::InputDomain MzChronogram::getInputDomain(void) const { 
   return TimeDomain; 
}



//////////////////////////////
//
// MzChronogram::getOutputDescriptors -- return a list describing
//    each of the available outputs for the object.  OutputList
//    is defined in the file vamp-sdk/Plugin.h:
//
// .name             == short name of output for computer use.  Must not
//                      contain spaces or punctuation.
// .description      == long name of output for human use.
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

MzChronogram::OutputList MzChronogram::getOutputDescriptors(void) const {

   OutputList       odlist;
   OutputDescriptor od;

   // First and only output channel:
   od.name             = "chronogram";
   od.description      = "Chronogram";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = getBlockSize() * getChannels();
   od.hasKnownExtents  = false;
   // od.minValue      = 0.0;
   // od.maxValue      = 0.0;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::OneSamplePerStep;
   // od.sampleRate    = 0.0;
   odlist.push_back(od);

   return odlist; 
}



//////////////////////////////
//
// MzChronogram::initialise -- this function is called once
//     before the first call to process().
//

bool MzChronogram::initialise(size_t channels, size_t stepsize, 
      size_t blocksize) {

   if (channels < getMinChannelCount() || channels > getMaxChannelCount()) {
      return false;
   }

   // step size and block size should never be zero
   if (stepsize <= 0 || blocksize <= 0) {
      return false;
   }

   mz_channels  = channels;
   mz_stepsize  = stepsize;
   mz_blocksize = blocksize;

   // Only one copy of a particular sample should be displayed.
   // If the step size is smaller than the block size, pretend
   // that the block size is the same as the step size.
   mz_blocksize = std::min(mz_stepsize, mz_blocksize);

   return true;
}



//////////////////////////////
//
// MzChronogram::process -- This function is called sequentially on the 
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

MzChronogram::FeatureSet MzChronogram::process(float **inputbufs, 
      Vamp::RealTime timestamp) {

   if (getStepSize() <= 0) {
      std::cerr << "ERROR: MzChronogram::process: "
                << "MzChronogram has not been initialized"
                << std::endl;
      return FeatureSet();
   }

   FeatureSet returnFeatures;
   Feature    feature;

   feature.hasTimestamp = false;

   // The Chronogram display has to be turned "upside-down" so that 
   // steeper downward slopes indicate flatter notes, and steeper 
   // higher slopes indicate sharper notes.
   int chan, samp;
   for (chan=getChannels()-1; chan>=0; chan--) {
      for (samp=getBlockSize()-1; samp>=0; samp--) {
         feature.values.push_back(inputbufs[chan][samp]);
      }
   }
   returnFeatures[0].push_back(feature);

   return returnFeatures;
}



//////////////////////////////
//
// MzChronogram::getRemainingFeatures -- This function is called
//    after the last call to process() on the input data stream has 
//    been completed.  Features which are non-causal can be calculated 
//    at this point.  See the comment above the process() function
//    for the format of output Features.
//

MzChronogram::FeatureSet MzChronogram::getRemainingFeatures(void) {
   // no remaining features, so return a dummy feature
   return FeatureSet();
}



//////////////////////////////
//
// MzChronogram::reset -- This function may be called after data processing
//    has been started with the process() function.  It will be called when
//    processing has been interrupted for some reason and the processing
//    sequence needs to be restarted (and current analysis output thrown out).  
//    After this function is called, process() will start at the beginning
//    of the input selection as if initialise() had just been called.
//    Note, however, that initialise() will NOT be called before processing 
//    is restarted after a reset().
//

void MzChronogram::reset(void) {
   // no actions necessary to reset this plugin
}


///////////////////////////////////////////////////////////////////////////
//
// Non-Interface Functions 
//

// no non-interface functions present here



