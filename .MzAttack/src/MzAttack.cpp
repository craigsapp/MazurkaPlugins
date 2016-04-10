//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun Jun 18 00:23:13 PDT 2006
// Last Modified: Fri Jul 14 23:00:01 PDT 2006
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzAttack.cpp
// URL:           http://sv.mazurka.org.uk/src/MzAttack.cpp
// Documentation: http://sv.mazurka.org.uk/MzAttack
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Identify note attacks in an audio signal.
//

// Defines used in getPluginVersion():
#define P_VER    "200607150"
#define P_NAME   "MzAttack"

#include "MzAttack.h" 

#include <stdio.h>
#include <string>
#include <math.h>

#include <algorithm>

#define METHOD_MAGNITUDE_PRODUCT   1
#define METHOD_MAGNITUDE_SUMMATION 2
#define METHOD_COMPLEX_SUMMATION   3

#define influence (0.1)
#define OFFSET    (-0.5)


///////////////////////////////////////////////////////////////////////////
//
// Vamp Interface Functions
//

///////////////////////////////
//
// MzAttack::MzAttack -- class constructor.
//

MzAttack::MzAttack(float samplerate) : 
      MazurkaPlugin(samplerate) {
   mz_harmonics     = 5;
   mz_transformsize = 16384;
   mz_minbin        = 0;
   mz_maxbin        = 511;
   mz_compress      = 0;
}



///////////////////////////////
//
// MzAttack::~MzAttack -- class destructor.
//

MzAttack::~MzAttack() {
   // do nothing
}


////////////////////////////////////////////////////////////
//
// parameter functions --
//

//////////////////////////////
//
// MzAttack::getParameterDescriptors -- return a list of
//      the parameters which can control the plugin.
//

MzAttack::ParameterList 
MzAttack::getParameterDescriptors(void) const {

   ParameterList       pdlist;
   ParameterDescriptor pd;

/*

   // first parameter: The number of samples in the audio window
   pd.identifier   = "windowsamples";
   pd.name         = "Window size";
   pd.unit         = "samples";
   pd.minValue     = 2.0;
   pd.maxValue     = 10000;
   pd.defaultValue = 1500.0;
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);

   // second parameter: The step size between analysis windows.
   pd.identifier   = "stepsamples";
   pd.name         = "Step size";
   pd.unit         = "samples";
   pd.minValue     = 2.0;
   pd.maxValue     = 30000.0;
   pd.defaultValue = int(getSrate() / 100.0 + 0.5);
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);

   // third parameter: The number of harmonics to consider
   pd.identifier   = "harmonics";
   pd.name         = "Harmonics";
   pd.unit         = "";
   pd.minValue     = 2.0;
   pd.maxValue     = 20.0;
   pd.defaultValue = 5.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);

   // fourth parameter: The minimum pitch to consider
   pd.identifier   = "minpitch";
   pd.name         = "Min pitch";
   pd.unit         = "MIDI data";
   pd.minValue     = 0.0;
   pd.maxValue     = 127.0;
   generateMidiNoteList(pd.valueNames, 0, 127);
   pd.defaultValue = 36.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // fifth parameter: The maximum pitch to consider
   pd.identifier   = "maxpitch";
   pd.name         = "Max pitch";
   pd.unit         = "MIDI data";
   pd.minValue     = 0.0;
   pd.maxValue     = 127.0;
   generateMidiNoteList(pd.valueNames, 0, 127);
   pd.defaultValue = 72.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // sixth parameter: The method for harmonic correlation
   pd.identifier   = "method";
   pd.name         = "Method";
   pd.unit         = "";
   pd.minValue     = 1.0;
   pd.maxValue     = 3.0;
   pd.valueNames.push_back("Magnitude Product");
   pd.valueNames.push_back("Magnitude Summation");
   pd.valueNames.push_back("Complex Summation");
   pd.defaultValue = 1.0; 
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // seventh parameter: Magnitude range compression.
   pd.identifier   = "compress";
   pd.name         = "Compress range";
   pd.unit         = "";
   pd.minValue     = 0.0;
   pd.maxValue     = 1.0;
   pd.defaultValue = 1.0;
   pd.valueNames.push_back("no");
   pd.valueNames.push_back("yes");
   pd.isQuantized  = true;
   pd.quantizeStep = 1.0;
   pdlist.push_back(pd);
   pd.valueNames.clear();

   // seventh parameter: Magnitude range compression.
   pd.identifier   = "smoothing";
   pd.name         = "Smoothing factor";
   pd.unit         = "";
   pd.minValue     = -1.0;
   pd.maxValue     = 1.0;
   pd.defaultValue = 0.1;
   pd.isQuantized  = false;
   pd.quantizeStep = 0.0;
   pdlist.push_back(pd);

*/

   return pdlist;
}


////////////////////////////////////////////////////////////
//
// optional polymorphic functions inherited from PluginBase:
//

//////////////////////////////
//
// MzAttack::getPreferredStepSize -- overrides the 
//     default value of 0 (no preference) returned in the 
//     inherited plugin class.
//

size_t MzAttack::getPreferredStepSize(void) const {
   return int(getSrate() / 100.0 + 0.5);
//   return getParameterInt("stepsamples");
}



//////////////////////////////
//
// MzAttack::getPreferredBlockSize -- overrides the 
//     default value of 0 (no preference) returned in the 
//     inherited plugin class.
//

size_t MzAttack::getPreferredBlockSize(void) const {
 
// int transformsize = getParameterInt("transformsamples");
// int blocksize     = getParameterInt("windowsamples");
   int transformsize = 8192;
   int blocksize     = 1500;

   if (blocksize > transformsize) {
      blocksize = transformsize;
   }

   return blocksize;
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzAttack::getIdentifier(void) const
   { return "mzattack"; }

std::string MzAttack::getName(void) const
   { return "Attack detector"; }

std::string MzAttack::getDescription(void) const
   { return "Attack detector"; }

std::string MzAttack::getMaker(void) const
   { return "The Mazurka Project"; }

std::string MzAttack::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

int MzAttack::getPluginVersion(void) const {
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
// MzAttack::getInputDomain -- the host application needs
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

MzAttack::InputDomain 
MzAttack::getInputDomain(void) const { 
   return TimeDomain; 
}



//////////////////////////////
//
// MzAttack::getOutputDescriptors -- return a list describing
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

MzAttack::OutputList 
MzAttack::getOutputDescriptors(void) const {

   OutputList       odlist;
   OutputDescriptor od;

   std::string s;

/*
   // First output channel: Attack Detection Spectrum
#define O_SPECTRUM X
   od.identifier       = "attackdetectionspectrum";
   od.name             = "Attack Detection Spectrum";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = mz_maxbin - mz_minbin + 1 - 2;
   od.hasKnownExtents  = false;   
   // od.minValue      = 0.0;
   // od.maxValue      = 1.0;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::OneSamplePerStep;
   // od.sampleRate    = 0.0;
   odlist.push_back(od);
*/

   // Second output channel: Attack Detection Function
#define O_FUNCTION 0
   od.identifier       = "attackdetectionfunction";
   od.name             = "Attack Detection Function";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = 1;
   od.hasKnownExtents  = false;   // could set to true.
   // od.minValue      = 0.0;
   // od.maxValue      = 1.0;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::VariableSampleRate;
   // od.sampleRate    = 0.0;
   odlist.push_back(od);

   // Third output channel: Detected Attacks
#define O_ATTACK 1
   od.identifier       = "detectedattacks";
   od.name             = "Detected Attacks";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = 0;
   od.hasKnownExtents  = false;   // could set to true.
   // od.minValue      = 0.0;
   // od.maxValue      = 1.0;
   od.isQuantized      = false;
   // od.quantizeStep  = 1.0;
   od.sampleType       = OutputDescriptor::VariableSampleRate;
   // od.sampleRate    = 0.0;
   odlist.push_back(od);

   return odlist; 
}



//////////////////////////////
//
// MzAttack::initialise -- this function is called once
//     before the first call to process().
//

bool MzAttack::initialise(size_t channels, size_t stepsize, 
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

   if (getBlockSize() > mz_transformsize) {
      setBlockSize(mz_transformsize);
   }

//   mz_method       = getParameterInt("method");
//   mz_harmonics    = getParameterInt("harmonics");
//   mz_compress     = getParameterInt("compress");
   mz_method       = 1;
   mz_harmonics    = 5;
   mz_compress     = 1;

   double minfreq, maxfreq, a440interval;

//   a440interval = getParameter("minpitch") - 67.0;
   a440interval = 36 - 67.0;
   minfreq = 440.0 * pow(2.0, a440interval / 12.0);
   mz_minbin = int(minfreq * mz_transformsize / getSrate());

//   a440interval = getParameter("maxpitch") - 67.0;
   a440interval = 72 - 67.0;
   maxfreq = 440.0 * pow(2.0, a440interval / 12.0);
   mz_maxbin = int(maxfreq * mz_transformsize / getSrate() + 0.999);

   if (mz_minbin > mz_maxbin) {
      std::swap(mz_minbin, mz_maxbin);
   }

   if (mz_maxbin >= mz_transformsize) {
      std::cerr << "MzAttack::initialize: maxbin size problem" 
                << std::endl;
      std::cerr << "MzAttack::initialize: maxbin = " 
                << mz_maxbin << std::endl;
      std::cerr << "MzAttack::initialize: transformsize = " 
                << mz_transformsize << std::endl;
      return false;
   }

   if (mz_minbin < 0) {
      std::cerr << "MzAttack::initialize: minbin size problem" 
                << std::endl;
      std::cerr << "MzAttack::initialize: minbin = " 
                << mz_minbin << std::endl;
      return false;
   }

   mz_transformer.setSize(mz_transformsize);
   mz_transformer.zeroSignal();
   mz_windower.setSize(getBlockSize());
   mz_windower.makeWindow("Hann");

   rawframes.setSize(20);      // must be larger than 5 for Gaussian filtering
   rawframes.reset();
   smoothedframes.setSize(20); // must be larger than 3 for Edge detecting
   smoothedframes.reset();

   return true;
}



//////////////////////////////
//
// MzAttack::process -- This function is called sequentially on the 
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

#define sigmoidscale(x,c,w)  (1.0/(1.0+exp(-((x)-(c))/((w)/8.0))))
#define NONPEAKFACTOR 0.2

MzAttack::FeatureSet MzAttack::process(AUDIODATA inputbufs, 
      Vamp::RealTime timestamp) {

   if (getStepSize() <= 0) {
      std::cerr << "ERROR: MzAttack::process: "
                << "MzAttack has not been initialized"
                << std::endl;
      return FeatureSet();
   }

   FeatureSet returnFeatures;
   Feature feature;

   feature.hasTimestamp = false;

   mz_windower.windowNonCausal(mz_transformer, inputbufs[0], getBlockSize());

   mz_transformer.doTransform();

   int bincount = mz_maxbin - mz_minbin + 1;

   std::vector<double> currentrawframe;
   currentrawframe.resize(bincount);

   std::vector<double> currentsmoothedframe;
   currentsmoothedframe.resize(bincount);

   std::vector<double> currentedgeframe;
   currentedgeframe.resize(bincount);

   int spectrumsize = mz_transformsize / 2;
   std::vector<double> magnitudespectrum(spectrumsize);
   std::vector<mz_complex> complexspectrum(spectrumsize);
   std::vector<int>    harmoniccount(bincount);

   int i, j;
   for (i=0; i<bincount; i++) {
      harmoniccount[i] = 0;
   }

   int topbin = mz_maxbin * mz_harmonics;
   if (topbin >= spectrumsize) {
      topbin = spectrumsize - 1;
   }

   int index;
   std::vector<int> maxpeak(spectrumsize);
   mz_complex complexsum;
   mz_complex&cs = complexsum;
   int maxvaluebin = 0;    
   double spectralpower = 0.0;

   switch (mz_method) {

      case METHOD_MAGNITUDE_SUMMATION:

         for (i=0; i<spectrumsize; i++) {
            magnitudespectrum[i] = mz_transformer.getSpectrumMagnitude(i);
            if (i > topbin) {
               // won't need the rest of the magnitude spectrum
               break;
            }
         }

         for (i=mz_minbin; i<=mz_maxbin; i++) {
            currentrawframe[i - mz_minbin] = 0.0;
            for (j=1; j<=mz_harmonics; j++) {
               index = i*j;
               if (index > spectrumsize) {
                  break;
               }
               currentrawframe[i - mz_minbin] += magnitudespectrum[index];
               harmoniccount[i - mz_minbin]++;
            }
         }

         // convert the harmonic spectrum to db
         for (i=0; i<bincount; i++) {
            if (currentrawframe[i] <= 0.0) {
               currentrawframe[i] = -120.0;
            } else {
               spectralpower += currentrawframe[i] / harmoniccount[i];
               currentrawframe[i] = 20.0 
		     * log10(currentrawframe[i] / harmoniccount[i]);
            }
            if (currentrawframe[i] > currentrawframe[maxvaluebin]) {
               maxvaluebin = i;
            }
         }

         break;

      case METHOD_COMPLEX_SUMMATION:

         for (i=0; i<spectrumsize; i++) {
            complexspectrum[i] = mz_transformer.getSpectrum(i);
            if (i > topbin) {
               // won't need the rest of the magnitude spectrum
               break;
            }
         }

         for (i=mz_minbin; i<=mz_maxbin; i++) {
	    complexsum.re = 0.0;
	    complexsum.im = 0.0;
            for (j=1; j<=mz_harmonics; j++) {
               index = i*j;
               if (index > spectrumsize) {
                  break;
               }
               complexsum.re +=  complexspectrum[index].re;
               complexsum.im +=  complexspectrum[index].im;
               harmoniccount[i - mz_minbin]++;
            }
            currentrawframe[i - mz_minbin] = sqrt(cs.re*cs.re + cs.im*cs.im);
         }

         // convert the harmonic spectrum to db
         for (i=0; i<bincount; i++) {
            if (currentrawframe[i] <= 0.0) {
               currentrawframe[i] = -120.0;
            } else {
               spectralpower += currentrawframe[i] / harmoniccount[i];
               currentrawframe[i] = 20.0 
		     * log10(currentrawframe[i] / harmoniccount[i]);
            }
            if (currentrawframe[i] > currentrawframe[maxvaluebin]) {
               maxvaluebin = i;
            }
         }

         break;

      case METHOD_MAGNITUDE_PRODUCT:
      default:

         for (i=0; i<spectrumsize; i++) {
            magnitudespectrum[i] = mz_transformer.getSpectrumMagnitude(i);
            if (i > topbin) {
               // won't need the rest of the magnitude spectrum
               break;
            }
         }

         for (i=mz_minbin; i<=mz_maxbin; i++) {
            currentrawframe[i - mz_minbin] = 1.0;
            for (j=1; j<=mz_harmonics; j++) {
               index = i*j;
               if (index > spectrumsize) {
                  break;
               }
               currentrawframe[i - mz_minbin] *= magnitudespectrum[index];
               harmoniccount[i - mz_minbin]++;
            }
         }

         // convert the harmonic spectrum to db
         for (i=0; i<bincount; i++) {
            if (currentrawframe[i] <= 0.0) {
               currentrawframe[i] = -120.0;
            } else {
               spectralpower += pow(currentrawframe[i], 1.0/harmoniccount[i]);
               currentrawframe[i] = 20.0 / harmoniccount[i] 
		                   * log10(currentrawframe[i]);
            }
            if (currentrawframe[i] > currentrawframe[maxvaluebin]) {
               maxvaluebin = i;
            }
         }

   }

   double cen;
   if (mz_compress) {
      for (i=0; i<bincount; i++) {
	 cen = -40.0 * i / bincount;
         currentrawframe[i] = 
               sigmoidscale(currentrawframe[i], cen, 60);
      }
   }



   // process the first, second, third (and fourth) outputs from the plugin:
   //

   rawframes.insert(currentrawframe);
   if (rawframes.getCount() >= 5) {
      calculateSmoothedFrame(currentsmoothedframe, rawframes);
   } else{
      currentsmoothedframe = currentrawframe;
   }
   smoothedframes.insert(currentsmoothedframe);

   if (smoothedframes.getCount() >= 3) {
      calculateEdgeFrame(currentedgeframe, smoothedframes);
   } else {
      currentedgeframe = currentsmoothedframe;
   }

   std::sort(currentedgeframe.begin(), currentedgeframe.end());

   feature.values.resize(bincount);
   for (i=0; i<(int)currentedgeframe.size(); i++) {
      feature.values[i] = currentedgeframe[i];
   }

//   returnFeatures[O_SPECTRUM].push_back(feature);

   if (rawframes.getCount() < 5) {
      attackfunction.push_back(0.0);
   } else {
      std::sort(currentedgeframe.begin(), currentedgeframe.end());

      int size = currentedgeframe.size();
      double value = 1.0;
      int counter  = 0;
      for (i=size-1; i>size - size/2 ; i--) {
         if (currentedgeframe[i] > 0.0) {
            value *= currentedgeframe[i];
            counter++;
         }
      }

      if (counter > size/8) {
         value = pow(value, 1.0/counter);
      } else {
         value = 0.0;
      }

      double value2 = calculateSlope(currentedgeframe);
      attackfunction.push_back(value * (value2 + value * influence));

      // using only the loudest peak in spectral frame
      // works well, given few false negatives, but does
      // give a lot of false positives
      // attackfunction.push_back(currentedgeframe[size-1]);
   }

   return returnFeatures;
}



//////////////////////////////
//
// MzAttack::getRemainingFeatures -- This function is called
//    after the last call to process() on the input data stream has 
//    been completed.  Features which are non-causal can be calculated 
//    at this point.  See the comment above the process() function
//    for the format of output Features.
//

MzAttack::FeatureSet MzAttack::getRemainingFeatures(void) {

   std::vector<double>& af = attackfunction;
//   double filk = getParameter("smoothing");
   double filk = 0.1;
   double rfilk = 1.0 - filk;
   int size = af.size();
   int i;

   std::vector<double> smoothed(size,true);

   // Difference equation for smoothing: y[n] = k x[n] + (1-k) y[n-1]


   // reverse filtering
   smoothed[size-1] = af[size-1];
   for (i=size-2; i>=0; i--) {
      smoothed[i] = filk * af[i] + rfilk * smoothed[i+1];
   }

   // forward filtering
   for (i=1; i<size; i++) {
      smoothed[i] = filk * smoothed[i] + rfilk * smoothed[i-1];
   }

   FeatureSet returnFeatures;
   Feature feature;
   Vamp::RealTime timestamp = Vamp::RealTime::fromSeconds(0.0);

   // process output #2: attack function
   
   feature.hasTimestamp = true;

   double max = 0.0;
   
   for (i=0; i<size; i++) {
      if (af[i] > max) {
         max = af[i];
      }
   }


   double value = 0.0;
   size = af.size();
   std::vector<double> naf;
   naf.resize(size);

   naf[0] = 0.0;
   naf[1] = 0.0;

   for (i=2; i<size; i++) {
      //value = af[i] - smoothed[i];
      value = af[i];
      if (value < 0.0) {
        value = 0.0;
      }
      value = log10(value/max * 90.0 + 10.0) - 1.0;
      naf[i] = value;
      feature.values.clear();
      feature.values.push_back(value);
      feature.timestamp = positionTimestampInStep(timestamp, i+OFFSET);
      returnFeatures[O_FUNCTION].push_back(feature);
   }

   // output #3: identify the attacks 

   int histeresis = 5;
   int lasttick = 0;
   feature.values.clear();
   for (i=2; i<size-1; i++) {
      if (naf[i] > naf[i-1] && naf[i] > naf[i+1] && naf[i] > 0.045) {
         if (i - lasttick  > histeresis) {
            feature.timestamp = positionTimestampInStep(timestamp, i+OFFSET);
            returnFeatures[O_ATTACK].push_back(feature);
            lasttick = i;
         }
      }
   }

   return returnFeatures;
}



//////////////////////////////
//
// MzAttack::reset -- This function may be called after data 
//    processing has been started with the process() function.  It will 
//    be called when processing has been interrupted for some reason and 
//    the processing sequence needs to be restarted (and current analysis 
//    output thrown out).  After this function is called, process() will 
//    start at the beginning of the input selection as if initialise() 
//    had just been called.  Note, however, that initialise() will NOT 
//    be called before processing is restarted after a reset().
//

void MzAttack::reset(void) {
   attackfunction.clear();

   rawframes.reset();
   smoothedframes.reset();
}


///////////////////////////////////////////////////////////////////////////
//
// Non-Interface Functions 
//


//////////////////////////////
//
// generateMidiNoteList -- Create a list of pitch names for the 
//   specified MIDI key number range.
//

void MzAttack::generateMidiNoteList(std::vector<std::string>& alist,
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
// MzAttack::calculateSmoothedFrame --
//
// http://www.pages.drexel.edu/~weg22/can_tut.html
// 
// A Gaussian mask.  The wider the mask , the greater the localization
// error, but the lower the noise.
// 
// 2  4  5  4  2
// 4  9 12  9  4
// 5 12 15 12  5
// 4  9 12  9  4
// 2  4  5  4  2
// 
// then divide by 115
// 

void MzAttack::calculateSmoothedFrame(std::vector<double>& output, 
      CircularBuffer<std::vector<double> >& input) {

   int i;
   int size = output.size();

   std::vector<double>& p0 = input[0];
   std::vector<double>& p1 = input[-1];
   std::vector<double>& p2 = input[-2];
   std::vector<double>& p3 = input[-3];
   std::vector<double>& p4 = input[-4];

   output[0]      = p2[0];
   output[1]      = p2[1];
   output[size-1] = p2[size-1];
   output[size-2] = p2[size-2];

   for (i=2; i<size-2; i++) {
      output[i] = 0.0;

      output[i] += p4[i+2] *  2.0;
      output[i] += p4[i+1] *  4.0;
      output[i] += p4[i+0] *  5.0;
      output[i] += p4[i-1] *  4.0;
      output[i] += p4[i-2] *  2.0;

      output[i] += p3[i+2] *  4.0;
      output[i] += p3[i+1] *  9.0;
      output[i] += p3[i+0] * 12.0;
      output[i] += p3[i-1] *  9.0;
      output[i] += p3[i-2] *  4.0;

      output[i] += p2[i+2] *  5.0;
      output[i] += p2[i+1] * 12.0;
      output[i] += p2[i+0] * 15.0;
      output[i] += p2[i-1] * 12.0;
      output[i] += p2[i-2] *  5.0;

      output[i] += p1[i+2] *  4.0;
      output[i] += p1[i+1] *  9.0;
      output[i] += p1[i+0] * 12.0;
      output[i] += p1[i-1] *  9.0;
      output[i] += p1[i-2] *  4.0;

      output[i] += p0[i+2] *  2.0;
      output[i] += p0[i+1] *  4.0;
      output[i] += p0[i+0] *  5.0;
      output[i] += p0[i-1] *  4.0;
      output[i] += p0[i-2] *  2.0;

      output[i] = output[i] / 115.0;

   }

}



//////////////////////////////
//
// MzAttack::calculateEdgeFrame --
//
// vertical edges:
//   -1  0 +1
//   -2  0 +2
//   -1  0 +1
//

void MzAttack::calculateEdgeFrame(std::vector<double>& output, 
      CircularBuffer<std::vector<double> >& input) {

   int i;
   int size = output.size();
   
   std::vector<double>& p0 = input[0];
   // std::vector<double>& p1 = input[-1];
   std::vector<double>& p2 = input[-2];

   output[0]      = 0.0;
   output[size-1] = 0.0;

   for (i=1; i<size-1; i++) {
      output[i] = 0.0;

      output[i] += p2[i+1] * -1.0;
      output[i] += p2[i+0] * -2.0;
      output[i] += p2[i-1] * -1.0;

      output[i] += p0[i+1] * +1.0;
      output[i] += p0[i+0] * +2.0;
      output[i] += p0[i-1] * +1.0;
   }

}



//////////////////////////////
//
// MzAttack::calculateSlope -- match a line to part of the data.
//

double MzAttack::calculateSlope(std::vector<double>& data) {
   int length = data.size();
   int N  = length / 8;

   double x;
   double y;
   double sumx  = 0.0;
   double sumxx = 0.0;
   double sumy  = 0.0;
   double sumxy = 0.0;

   int i;
   for (i=0; i<N; i++) {
      x = i;
      y = data[length - i];

      sumx += x;
      sumy += y;
      sumxx += x*x;
      sumxy += x*y;
   }

   double denom = sumxx - sumx * sumx / N;
   if (denom != 0.0) {
      denom += 0.0001;
   }
   
   double slope = (sumxy - sumx * sumy / N) / denom;

   return -data[length-1] * (slope - data[length-1] * influence);
}









