//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun May  7 18:19:57 PDT 2006
// Last Modified: Sun May  7 21:33:33 PDT 2006
// Filename:      MzSpectrogramHost.h
// Syntax:        C++; vamp plugin
//

#ifndef _MZSPECTRUMHOST_H_INCLUDED
#define _MZSPECTRUMHOST_H_INCLUDED

#include "Plugin.h"   // Vamp plugin interface for Sonic Visualiser


class MzSpectrogramHost : public Vamp::Plugin {

   public: // plugin interface functions:

                    MzSpectrogramHost    (float inputSampleRate);
      virtual      ~MzSpectrogramHost    ();

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

      void          makeMagnitudeSpectrum(float** data, size_t channel,
                                          size_t blocksize);
      int           getBinCount          (void) const;

   protected: // non-interface functions and variables:


      size_t mz_blocksize;    // number of samples in each input
      size_t mz_stepsize;     // sample hop between blocks in input
      size_t mz_channels;     // number of samples in the input

      double* mz_workbuffer;  // storage space for transform calculations
};



#endif // _MZSPECTRUMHOST_H_INCLUDED

