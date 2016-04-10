//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat May 13 12:19:13 PDT 2006
// Last Modified: Sat May 20 15:24:51 PDT 2006 (added parameter control)
// Filename:      MzPowerCurve.cpp
// URL:           http://sv.mazurka.org.uk/MzPowerCurve/src/MzPowerCurve.cpp
// Documentation: http://sv.mazurka.org.uk/MzPowerCurve
// Syntax:        ANSI C++; vamp 0.9 plugin
//
// Description:   Calculate the power of an audio signal.
//

#include "MzPowerCurve.h" 
#include <math.h>

using std::vector;
using std::list;


///////////////////////////////////////////////////////////////////////////
//
// Vamp Interface Functions
//

///////////////////////////////
//
// MzPowerCurve::MzPowerCurve -- class constructor.
//

MzPowerCurve::MzPowerCurve(float samplerate) : MazurkaPlugin(samplerate) {
   initializeParameterDatabase(getParameterDescriptors());
   mz_blocksize = 0;
   mz_stepsize  = 0;
   mz_channels  = 0;
}



///////////////////////////////
//
// MzPowerCurve::~MzPowerCurve -- class destructor.
//

MzPowerCurve::~MzPowerCurve() {
   // do nothing
}


////////////////////////////////////////////////////////////
//
// Initialization parameter functions --
//

//////////////////////////////
//
// MzTimeogram::getParameterDescriptors -- return a list of
//      the parameters which can control the plugin.
//

MzPowerCurve::ParameterList MzPowerCurve::getParameterDescriptors(void) const {

   ParameterList       pdlist;
   ParameterDescriptor pd;

   // first parameter: The size of the analysis window in milliseconds
   pd.name         = "windowsize";
   pd.description  = "Window size";
   pd.unit         = "ms";
   pd.minValue     = 1.0;
   pd.maxValue     = 1000000.0;
   pd.defaultValue = 10.0;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   // second parameter: The hop size between windows in milliseconds
   pd.name         = "hopsize";
   pd.description  = "Window hop size";
   pd.unit         = "ms";
   pd.minValue     = 1.0;
   pd.maxValue     = 1000000.0;
   pd.defaultValue = 10.0;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   // third parameter: Cut-off threshold for scaled power slope
   pd.name         = "cutoffthreshold";
   pd.description  = "Cut-off threshold";
   pd.unit         = "dB";
   pd.minValue     = -100.0;
   pd.maxValue     = 10.0;
   pd.defaultValue = -50.0;
   pd.isQuantized  = 0;
   // pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

   // fourth parameter: The width of the cut-off transition region
   pd.name         = "cutoffwidth";
   pd.description  = "Cut-off transition width";
   pd.unit         = "dB";
   pd.minValue     = 1.0;
   pd.maxValue     = 100.0;
   pd.defaultValue = 10.0;
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
// MzPowerCurve::getPreferredStepSize --
//

size_t MzPowerCurve::getPreferredStepSize(void) const { 
   return size_t(getParameter("hopsize")*getSrate() + 0.5);
}



/////////////////////////////
//
// MzPowerCurve::getPreferredBlockSize --
//

size_t MzPowerCurve::getPreferredBlockSize(void) const { 
   return size_t(getParameter("windowsize")*getSrate() + 0.5);
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzPowerCurve::getName(void) const
   { return "mzpowercurve"; }

std::string MzPowerCurve::getMaker(void) const
   { return "Craig Stuart Sapp <craig@ccrma.stanford.edu>"; }

std::string MzPowerCurve::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

std::string MzPowerCurve::getDescription(void) const
   { return "Power Curve"; }

int MzPowerCurve::getPluginVersion(void) const
   { return 200605210; }


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from Plugin:
//

//////////////////////////////
//
// MzPowerCurve::getInputDomain -- the host application needs
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

MzPowerCurve::InputDomain MzPowerCurve::getInputDomain(void) const { 
   return TimeDomain; 
}



//////////////////////////////
//
// MzPowerCurve::getOutputDescriptors -- return a list describing
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

MzPowerCurve::OutputList MzPowerCurve::getOutputDescriptors(void) const {

   OutputList       list;
   OutputDescriptor od;

   // First output channel:
   od.name             = "rawpower";
   od.description      = "Raw Power";
   od.unit             = "dB";
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

   // Second output channel: (smoothed data)
   od.name             = "smoothpower";
   od.description      = "Smoothed Power";
   od.unit             = "dB";
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

   // Third output channel: (smoothed data slope)
   od.name             = "smoothpowerslope";
   od.description      = "Smoothed Power Slope";
   od.unit             = "dB";
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

   // Fourth output channel: (smoothed data slope product)
   od.name             = "powerslopeproduct";
   od.description      = "Power Slope Product";
   od.unit             = "decibelish";
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
// MzPowerCurve::initialise -- this function is called once
//     before the first call to process().
//

bool MzPowerCurve::initialise(size_t channels, size_t stepsize, 
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

   return true;
}



//////////////////////////////
//
//
// MzPowerCurve::process -- This function is called sequentially on the 
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

MzPowerCurve::FeatureSet MzPowerCurve::process(float **inputbufs, 
      Vamp::RealTime timestamp) {

   if (getStepSize() <= 0) {
      std::cerr << "ERROR: MzPowerCurve::process: "
                << "MzPowerCurve has not been initialized"
                << std::endl;
      return FeatureSet();
   }

   float power = float(getPower(inputbufs[0], getBlockSize()));

   FeatureSet returnFeatures;
   Feature    feature;

   // center the location of the power measurement at the
   // middle of the analysis region rather than at the beginning.
   feature.hasTimestamp = true;
   feature.timestamp    = timestamp + 
         Vamp::RealTime::fromSeconds(0.5 * getBlockSize()/getSrate());

   feature.values.push_back(power);
   // also store the power measurement for later processing in
   // getRemainingFeatures():
   rawpower.push_back(power);

   returnFeatures[0].push_back(feature);

   return returnFeatures;
}



//////////////////////////////
//
// MzPowerCurve::getRemainingFeatures -- This function is called
//    after the last call to process() on the input data stream has 
//    been completed.  Features which are non-causal can be calculated 
//    at this point.  See the comment above the process() function
//    for the format of output Features.
//

MzPowerCurve::FeatureSet MzPowerCurve::getRemainingFeatures(void) {

   int i;
   double filterk = 0.3;
   double oneminusk = 1.0 - filterk;
   int size = rawpower.size();
   vector<double> smoothpower(size, true);

   // Difference equation for smoothing: y[n] = k * x[n] + (1-k) * y[n-1]

   // do forward smoothing:
   list<double>::iterator rawiter = rawpower.begin();
   smoothpower[0] = *rawiter;
   for (i=1; i<size; i++) {
      smoothpower[i] = filterk * *(++rawiter) + oneminusk * smoothpower[i-1];
   }

   // do reverse smoothing: time symmetric filtering to remove filter delays
   for (i=size-2; i>=0; i--) {
      smoothpower[i] = filterk * smoothpower[i] + oneminusk * smoothpower[i+1];
   }

   FeatureSet returnFeatures;
   Feature    feature;
   feature.hasTimestamp = true;

   // process output features #2: smoothed power data

   double timeinsec;
   for (i=0; i<size; i++) {
      timeinsec = (getBlockSize()/2.0 + i * getStepSize())/getSrate();
      feature.timestamp = Vamp::RealTime::fromSeconds(timeinsec);
      feature.values.clear();
      feature.values.push_back(float(smoothpower[i]));
      returnFeatures[1].push_back(feature);
   }

   // process output features #3 here: smooth power slope
   
   vector<double> smoothslope(size-1, true);
   for (i=0; i<size-1; i++) {
      smoothslope[i] = smoothpower[i+1] - smoothpower[i];
      timeinsec = (getBlockSize() + (2*i+1)*getStepSize())/(2.0*getSrate());
      feature.timestamp = Vamp::RealTime::fromSeconds(timeinsec);
      feature.values.clear();
      feature.values.push_back(float(smoothslope[i]));
      returnFeatures[2].push_back(feature);
   }

   // process output features #4 here: smooth power slope product
   
   vector<double> productslope(size-1, true);
   double cutoff = getParameter("cutoffthreshold");
   double width  = getParameter("cutoffwidth")/2.0;
   double scaling;
   for (i=0; i<size-1; i++) {
      scaling = (smoothpower[i] - cutoff)/width;
      scaling = 1.0 / (1.0 + pow(2.718281828, -scaling)); //sigmoid function
      productslope[i] = smoothslope[i] * scaling;
      timeinsec = (getBlockSize() + (2*i+1)*getStepSize())/(2.0*getSrate());
      // check on timinsec value relation to smoothpower and adjust
      feature.timestamp = Vamp::RealTime::fromSeconds(timeinsec);
      feature.values.clear();
      feature.values.push_back(float(productslope[i]));
      returnFeatures[3].push_back(feature);
   }

   return returnFeatures;
}



//////////////////////////////
//
// MzPowerCurve::reset -- This function may be called after data processing
//    has been started with the process() function.  It will be called when
//    processing has been interrupted for some reason and the processing
//    sequence needs to be restarted (and current analysis output thrown out).  
//    After this function is called, process() will start at the beginning
//    of the input selection as if initialise() had just been called.
//    Note, however, that initialise() will NOT be called before processing 
//    is restarted after a reset().
//

void MzPowerCurve::reset(void) {
   rawpower.clear();
}


///////////////////////////////////////////////////////////////////////////
//
// Non-Interface Functions 
//


//////////////////////////////
//
// MzPowerCurve::getPower -- find the average power of an audio block.
//

#define ZEROLOG -120.0

double MzPowerCurve::getPower(float* data, int blocksize) {

   double sum = 0.0; // sample energy accumulator

   for (int i=0; i<blocksize; i++) {
      sum += data[i] * data[i];
   }

   if (sum <= 0.0) {
      return ZEROLOG;
   } else {
      return 10.0 * log10(sum/blocksize);
   }
}



