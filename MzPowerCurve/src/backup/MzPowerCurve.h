//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Fri May  5 21:48:47 PDT 2006
// Last Modified: Sat May  6 08:28:57 PDT 2006
// Filename:      MzPowerCurve.h
// Syntax:        C++; vamp plugin
//

#ifndef _MZPOWERCURVE_H_INCLUDED
#define _MZPOWERCURVE_H_INCLUDED

#include "Plugin.h"   // Vamp plugin interface for Sonic Visualiser

class MzPowerCurve;

class MzPowerCurve : public Vamp::Plugin {

   public: // plugin interface functions:

                    MzPowerCurve         (float inputSampleRate);
      virtual      ~MzPowerCurve         ();

      // required polymorphic functions inherited from PluginBase:
      std::string   getName              (void) const;
      std::string   getDescription       (void) const;
      std::string   getMaker             (void) const;
      std::string   getCopyright         (void) const;
      int           getPluginVersion     (void) const;

      // required polymorphic functions inherited from Plugin:
      InputDomain   getInputDomain       (void) const;
      OutputList    getOutputDescriptors (void) const;
      FeatureSet    getRemainingFeatures (void);
      FeatureSet    process              (float **inputBuffers, 
                                          Vamp::RealTime timestamp);
      
      // optional polymorphic functions inherited from Plugin:
      void          reset                (void);
      bool          initialise           (size_t channels, size_t stepSize, 
                                          size_t blockSize);

      // non-interface functions and variables:

      static double getPower             (float** data, size_t channels, 
                                          size_t blocksize, int analchan = 0);

   protected: // non-interface functions and variables:


      size_t mz_blocksize;    // number of samples in each input
      size_t mz_stepsize;     // sample hop between blocks in input
      size_t mz_channels;     // number of samples in the input

};



#endif // _MZPOWERCURVE_H_INCLUDED

