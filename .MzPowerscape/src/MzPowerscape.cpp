//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat Jul  8 00:22:00 PDT 2006 (derived from MzPowerCurve)
// Last Modified: Sat Jul  8 00:22:13 PDT 2006
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzPowerscape.cpp
// URL:           http://sv.mazurka.org.uk/src/MzPowerscape.cpp
// Documentation: http://sv.mazurka.org.uk/MzPowerscape
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Calculate the power of an audio signal as it changes 
//                over time.
//

// Defines used in getPluginVersion():
#define P_VER    "200607080"
#define P_NAME   "MzPowerscape"

#include "MzPowerscape.h" 

#include <math.h>
#include <stdlib.h>

#define ZEROLOG -120.0

#define FILTER_NONE		0
#define FILTER_SYMMETRIC	1
#define FILTER_FORWARD		2
#define FILTER_BACKWARD		3


///////////////////////////////////////////////////////////////////////////
//
// Vamp Interface Functions
//

///////////////////////////////
//
// MzPowerscape::MzPowerscape -- class constructor.
//

MzPowerscape::MzPowerscape(float samplerate) : MazurkaPlugin(samplerate) {
   mz_levels = 0;
}



///////////////////////////////
//
// MzPowerscape::~MzPowerscape -- class destructor.
//

MzPowerscape::~MzPowerscape() {
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
// MzPowerscape::getParameterDescriptors -- return a list of
//      the parameters which can control the plugin.
//

MzPowerscape::ParameterList MzPowerscape::getParameterDescriptors(void) const {

   ParameterList       pdlist;
   ParameterDescriptor pd;

   // first parameter: The size of the analysis window in milliseconds
   // the hop size for the window is equivalent to its size.
   pd.identifier   = "windowsize";
   pd.name         = "Window size";
   pd.unit         = "seconds";
   pd.minValue     = 0.01;
   pd.maxValue     = 10.0;
   pd.defaultValue = 1.0;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   // second parameter: The hierarchy levels in the display
   pd.identifier   = "levels";
   pd.name         = "Vertical Levels";
   pd.unit         = "";
   pd.minValue     = 1;
   pd.maxValue     = 10000.0;
   pd.defaultValue = 100.0;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   // third parameter: Windowing method
   // for the hierarchical part, not for the weighting of 
   // instantaneous power in individual cells as is done in MzPowerCurve.
   pd.identifier   = "window";
   pd.name         = "Weighting window";
   pd.unit         = "";
   pd.minValue     = 1.0;
   //MazurkaWindower::getWindowList(pd.valueNames);
   //pd.maxValue     = pd.valueNames.size();
   pd.maxValue     = 5.0;
   pd.defaultValue = 1.0;
   pd.isQuantized  = 1;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // fourth parameter: Factor for exponential smoothing filter
   pd.identifier   = "smoothingfactor";
   pd.name         = "Smoothing";
   pd.unit         = "";
   pd.minValue     = -1.0;
   pd.maxValue     = 1.0;
   pd.defaultValue = 0.2;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   // fifth parameter: Filtering method
   pd.identifier   = "filtermethod";
   pd.name         = "Filter method";
   pd.unit         = "";
   pd.minValue     = 0.0;
   pd.maxValue     = 3.0;
   pd.defaultValue = 0.0;
   pd.isQuantized  = 1;
   pd.quantizeStep = 1.0;
   pd.valueNames.push_back("None");
   pd.valueNames.push_back("Symmetric");
   pd.valueNames.push_back("Forward");
   pd.valueNames.push_back("Reverse");
   pdlist.push_back(pd);
   pd.valueNames.clear();

   return pdlist;
}


////////////////////////////////////////////////////////////
//
// optional polymorphic functions inherited from Plugin:
//

/////////////////////////////
//
// MzPowerscape::getPreferredStepSize --
//

size_t MzPowerscape::getPreferredStepSize(void) const { 
   return getPreferredBlockSize();
}



/////////////////////////////
//
// MzPowerscape::getPreferredBlockSize --
//

size_t MzPowerscape::getPreferredBlockSize(void) const { 
   return size_t(getParameter("windowsize")*getSrate() + 0.5);
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzPowerscape::getIdentifier(void) const
   { return "mzpowerscape"; }

std::string MzPowerscape::getName(void) const
   { return "Powerscape"; }

std::string MzPowerscape::getDescription(void) const
   { return "Powerscape"; }

std::string MzPowerscape::getMaker(void) const
   { return "The Mazurka Project"; }

std::string MzPowerscape::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

int MzPowerscape::getPluginVersion(void) const {
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
// MzPowerscape::getInputDomain -- the host application needs
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

MzPowerscape::InputDomain MzPowerscape::getInputDomain(void) const { 
   return TimeDomain; 
}



//////////////////////////////
//
// MzPowerscape::getOutputDescriptors -- return a list describing
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

MzPowerscape::OutputList MzPowerscape::getOutputDescriptors(void) const {

   OutputList       list;
   OutputDescriptor od;

   // First output channel:
   od.identifier       = "powerscape";
   od.name             = "Powerscape";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = mz_levels;
   od.hasKnownExtents  = false;
   // od.minValue      = 0.0;
   // od.maxValue      = 0.0;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::FixedSampleRate;
   od.sampleRate       = getStepSize() / getSrate();
   list.push_back(od);
   
   return list; 
}



//////////////////////////////
//
// MzPowerscape::initialise -- this function is called once
//     before the first call to process().
//

bool MzPowerscape::initialise(size_t channels, size_t stepsize, 
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

   mz_levels = getParameterInt("levels");

/*
   for (int i=0; i<mz_levels; i++) {
      mz_windows[i].makeWindow(getParameterString("window"),i*2+3);
   }

   mz_window.makeWindow(getParameterString("window"), getBlockSize());
*/

   switch (getParameterInt("filtermethod")) {
      case FILTER_FORWARD:
         mz_filterforward  = 1;
         mz_filterbackward = 0;
         break;
      case FILTER_BACKWARD:
         mz_filterforward  = 0;
         mz_filterbackward = 1;
         break;
      case FILTER_SYMMETRIC:
         mz_filterforward  = 1;
         mz_filterbackward = 1;
      case FILTER_NONE:
      default:
         mz_filterforward  = 0;
         mz_filterbackward = 0;
   }

   surfacepower.clear();

   return true;
}



//////////////////////////////
//
//
// MzPowerscape::process -- This function is called sequentially on the 
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

MzPowerscape::FeatureSet MzPowerscape::process(AUDIODATA inputbufs, 
      Vamp::RealTime timestamp) {

   if (getChannelCount() <= 0) {
      std::cerr << "ERROR: MzPowerscape::process: "
                << "MzPowerscape has not been initialized"
                << std::endl;
      return FeatureSet();
   }

   // calculate the raw power for the given input block of audio.
   // Further processing of the raw power data will be done in
   // the getRemainingFeatures() function.

   int    i;
   double value;
   double sum = 0.0;

   for (i=0; i<getBlockSize(); i++) {
      value = inputbufs[0][i];
      sum += value * value;
   }

   float power = sum / getBlockSize(); 
   // conversion to decibels will be done later

   // store the power measurement for later processing in
   // getRemainingFeatures():
   surfacepower.push_back(power);

   // return a dummy feature since features will be calculated
   // later.

   return FeatureSet();
}



//////////////////////////////
//
// MzPowerscape::getRemainingFeatures -- This function is called
//    after the last call to process() on the input data stream has 
//    been completed.  Features which are non-causal can be calculated 
//    at this point.  See the comment above the process() function
//    for the format of output Features.
//

MzPowerscape::FeatureSet MzPowerscape::getRemainingFeatures(void) {
//return FeatureSet();
   int i;
   double filk = getParameter("smoothingfactor");
   double oneminusk = 1.0 - filk;
   int size = surfacepower.size();
   std::vector<double> smoothpower(size, true);

   // Difference equation for smoothing: y[n] = k * x[n] + (1-k) * y[n-1]

   // do reverse smoothing: time symmetric filtering to remove filter delays
   if (mz_filterbackward && mz_filterforward) {
      // reverse filtering first 
      smoothpower[size-1] = surfacepower[size-1];
      for (i=size-2; i>=0; i--) {
         smoothpower[i] = filk*surfacepower[i] + oneminusk*smoothpower[i+1];
      }
      // then forward filtering
      for (i=1; i<size; i++) {
         smoothpower[i] = filk*smoothpower[i] + oneminusk*smoothpower[i-1];
      }
   } else if (mz_filterbackward) {
      smoothpower[size-1] = surfacepower[size-1];
      for (i=size-2; i>=0; i--) {
         smoothpower[i] = filk * surfacepower[i] + oneminusk * smoothpower[i+1];
      }
   } else if (mz_filterforward) {
      // do forward smoothing:
      smoothpower[0] = surfacepower[0];
      for (i=1; i<size; i++) {
         smoothpower[i] = filk * surfacepower[i] + oneminusk * smoothpower[i-1];
      }
   } else {
      smoothpower = surfacepower;
   }

   FeatureSet returnFeatures;
   Feature feature;
   feature.values.resize(mz_levels);
   feature.hasTimestamp = true;

   double value;
   int j;
   for (i=0; i<(int)smoothpower.size(); i++) {
      for (j=0; j<mz_levels; j++) {
         value = getPowerLevel(i, j, smoothpower);   // change later
         feature.values[j] = value;
	 feature.timestamp = Vamp::RealTime::fromSeconds(i * getStepSize() / 
			                                           getSrate());
      }
      returnFeatures[0].push_back(feature);
   }

   return returnFeatures;
}



//////////////////////////////
//
// MzPowerscape::reset -- This function may be called after data processing
//    has been started with the process() function.  It will be called when
//    processing has been interrupted for some reason and the processing
//    sequence needs to be restarted (and current analysis output thrown out).  
//    After this function is called, process() will start at the beginning
//    of the input selection as if initialise() had just been called.
//    Note, however, that initialise() will NOT be called before processing 
//    is restarted after a reset().
//

void MzPowerscape::reset(void) {
   surfacepower.clear();
}



///////////////////////////////////////////////////////////////////////////
//
// Non-Interface Functions 
//


//////////////////////////////
//
// getPowerLevel -- return the power for a given level.
//

double MzPowerscape::getPowerLevel(int i, int j, std::vector<double>& power) {
   int size = j * 2 + 1;
   int starti = i - size / 2;
   if (starti < 0) {
      //return -i*1000 - j;               // cell serial number
      return ZEROLOG;
   }
   if (starti + size > (int)power.size()) {
      //return -i*1000 - j;               // cell serial number
      return ZEROLOG;
   }

   double sum = 0.0;
   for (int ii=0; ii<size; ii++) {
      sum += power[ii + starti];

   }

   //return i*1000 + j;                  // cell serial number
   return 10.0 * log10(sum / size);
}



