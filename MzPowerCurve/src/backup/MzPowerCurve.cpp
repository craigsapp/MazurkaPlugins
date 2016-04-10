//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Fri May  5 21:48:47 PDT 2006
// Last Modified: Sat May  6 08:29:17 PDT 2006
// Filename:      MzPowerCurve.cpp
// Syntax:        C++; vamp plugin
//


#include "MzPowerCurve.h" 
#include <math.h>


///////////////////////////////
//
// MzPowerCurve::MzPowerCurve -- class constructor.
//

MzPowerCurve::MzPowerCurve(float inputSampleRate) :
      Plugin(inputSampleRate) {
   reset();
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
// polymorphic functions inherited from Plugin:
//

void MzPowerCurve::reset(void) {
   initialise(0, 0, 0);
}

bool MzPowerCurve::initialise(size_t channels, size_t stepSize, 
      size_t blockSize) {
   if (channels < getMinChannelCount() ||
       channels > getMaxChannelCount()) {
      return false;
   }

   mz_channels  = channels;
   mz_stepsize  = stepSize;
   mz_blocksize = blockSize;

   return true;
}


////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzPowerCurve::getName(void) const
   { return "mzpowercurve"; }

std::string MzPowerCurve::getDescription(void) const
   { return "Power Curve"; }

std::string MzPowerCurve::getMaker(void) const
   { return "Craig Stuart Sapp <craig@ccrma.stanford.edu>"; }

std::string MzPowerCurve::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

int MzPowerCurve::getPluginVersion(void) const
   { return 200605060; }



////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from Plugin:
//

MzPowerCurve::InputDomain MzPowerCurve::getInputDomain(void) const 
   { return TimeDomain; }


MzPowerCurve::OutputList MzPowerCurve::getOutputDescriptors(void) const {
   OutputList list;

   OutputDescriptor od;
   od.name             = "power";
   od.description      = "Raw Power";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = 1;
   od.hasKnownExtents  = false;
   od.isQuantized      = false;
   od.sampleType       = OutputDescriptor::VariableSampleRate;
   list.push_back(od);

   return list; 
}


MzPowerCurve::FeatureSet MzPowerCurve::getRemainingFeatures(void) {
   // no remaining features
   return FeatureSet();
}


MzPowerCurve::FeatureSet MzPowerCurve::process(float **inputBuffers, 
      Vamp::RealTime timestamp) {

   if (mz_stepsize <= 0) {
      std::cerr << "ERROR: MzPowerCurve::process: "
                << "MzPowerCurve has not been initialized"
                << std::endl;
      return FeatureSet();
   }

   double avgpower = getPower(inputBuffers, mz_channels, mz_blocksize);

   FeatureSet returnFeatures;
   Feature feature;

   feature.hasTimestamp = true;
   feature.timestamp    = timestamp;
   feature.values.push_back(avgpower);
   returnFeatures[0].push_back(feature);

   return returnFeatures;
}



////////////////////////////////////////////////////////////
//
// non-interface functions 
//

