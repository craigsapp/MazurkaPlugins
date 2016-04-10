//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun May  7 18:53:16 PDT 2006
// Last Modified: Sun May  7 18:53:19 PDT 2006
// Filename:      MzSpectrogramHost.cpp
// Syntax:        C++; vamp plugin
//


#include "MzSpectrogramHost.h" 
#include "math.h"



///////////////////////////////
//
// MzSpectrogramHost::MzSpectrum -- class constructor.
//

MzSpectrogramHost::MzSpectrogramHost(float inputSampleRate) :
      Plugin(inputSampleRate) {
   mz_workbuffer = NULL;
   reset();
}



///////////////////////////////
//
// MzSpectrogramHost::~MzSpectrum -- class destructor.
//

MzSpectrogramHost::~MzSpectrogramHost() {
   if (mz_workbuffer != NULL) {
      delete [] mz_workbuffer;
      mz_workbuffer = NULL;
   }
}


////////////////////////////////////////////////////////////
//
// polymorphic functions inherited from Plugin:
//

void MzSpectrogramHost::reset(void) {
   if (mz_workbuffer != NULL) {
      delete [] mz_workbuffer;
      mz_workbuffer = NULL;
   }
   initialise(0, 0, 0);
}

bool MzSpectrogramHost::initialise(size_t channels, size_t stepSize, 
      size_t blockSize) {
   if (channels < getMinChannelCount() ||
       channels > getMaxChannelCount()) {
      return false;
   }

   if (blockSize <= 0) {
      return false;
   }

   mz_channels  = channels;
   mz_stepsize  = stepSize;
   mz_blocksize = blockSize;

   if (mz_workbuffer != NULL) {
      delete [] mz_workbuffer;
      mz_workbuffer = NULL;
   }
   mz_workbuffer = new double[mz_blocksize/4];
   return true;
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzSpectrogramHost::getName(void) const
   { return "mzspectrumhost"; }

std::string MzSpectrogramHost::getDescription(void) const
   { return "Spectrum Host"; }

std::string MzSpectrogramHost::getMaker(void) const
   { return "Craig Stuart Sapp <craig@ccrma.stanford.edu>"; }

std::string MzSpectrogramHost::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

int MzSpectrogramHost::getPluginVersion(void) const
   { return 200605070; }



////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from Plugin:
//

MzSpectrogramHost::InputDomain MzSpectrogramHost::getInputDomain(void) const 
   { return FrequencyDomain; }


MzSpectrogramHost::OutputList 
MzSpectrogramHost::getOutputDescriptors(void) const {
   OutputList list;

   OutputDescriptor od;
   od.name             = "magintude";
   od.description      = "Magnitude Spectrum";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = getBinCount();
   od.hasKnownExtents  = false;
   od.isQuantized      = false;
   od.sampleType       = OutputDescriptor::OneSamplePerStep;
   list.push_back(od);

   return list; 
}


MzSpectrogramHost::FeatureSet MzSpectrogramHost::getRemainingFeatures(void) {
   // no remaining features
   return FeatureSet();
}


MzSpectrogramHost::FeatureSet MzSpectrogramHost::process(float **inputBuffers, 
      Vamp::RealTime timestamp) {

   if (mz_stepsize <= 0 || mz_workbuffer == NULL) {
      std::cerr << "ERROR: MzSpectrogramHost::process: "
                << "MzSpectrogramHost has not been initialized"
                << std::endl;
      return FeatureSet();
   }

   makeMagnitudeSpectrum(inputBuffers, mz_channels, mz_blocksize);

   FeatureSet returnFeatures;
   Feature feature;

   feature.hasTimestamp = false;
   for (int i=0; i<getBinCount(); i++) {
      feature.values.push_back(mz_workbuffer[i]);
   }
   returnFeatures[0].push_back(feature);

   return returnFeatures;
}



////////////////////////////////////////////////////////////
//
// non-interface functions 
//


//////////////////////////////
//
// makeMagnitudeSpectrum -- 
//

#define MAGNITUDE(x, y) sqrt((x)*(x)+(y)*(y))
#define ZEROLOG         120.0

void MzSpectrogramHost::makeMagnitudeSpectrum(float** data, 
      size_t channels, size_t blocksize) {
   int i;
   int count = getBinCount();
   double real, imag;
   for (i=0; i<count; i++) {
      real = data[0][i<<i];
      imag = data[0][(i<<i) + 1];
      mz_workbuffer[i] = MAGNITUDE(real, imag);

      // convert the amplitude magnitude into decibels
      if (mz_workbuffer[i] <= 0.0) {
         mz_workbuffer[i] = ZEROLOG;
      }  else {
         mz_workbuffer[i] = 20.0 * log10(mz_workbuffer[i]);
      }
   }
}



////////////////////////////////
///
// getBinCount -- return the number of channels in each output sample
//

int MzSpectrogramHost::getBinCount(void) const {
   return mz_blocksize / 4;
}



