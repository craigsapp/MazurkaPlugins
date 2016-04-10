//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun May  7 18:19:57 PDT 2006
// Last Modified: Sun May  7 21:33:33 PDT 2006
// Filename:      MzSpectrogramClient.h
// Syntax:        C++; vamp plugin
//

#ifndef _MZSPECTRUMCLIENT_H_INCLUDED
#define _MZSPECTRUMCLIENT_H_INCLUDED

#include "Plugin.h"   // Vamp plugin interface for Sonic Visualiser


class MzSpectrogramClient : public Vamp::Plugin {

   public: // plugin interface functions:

                    MzSpectrogramClient  (float inputSampleRate);
      virtual      ~MzSpectrogramClient  ();

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
      FeatureSet    process              (float **input, 
                                          Vamp::RealTime timestamp);
      
      // optional polymorphic functions inherited from Plugin:
      void          reset                (void);
      bool          initialise           (size_t channels, size_t stepSize, 
                                          size_t blockSize);

      // non-interface functions and variables:

      void          makeMagnitudeSpectrum(double* output, double* realfreq, 
                                          double* imagfreq, size_t blocksize);
      int           getBinCount          (void) const;

      // non-interface functions and variables:

      static void   makeHannWindow       (double* output, int blocksize);
      static void   windowSignal         (double* output, double* window, 
                                          float* input, int blocksize);
      static void   fft                  (int n, bool inverse, 
                                          double *ri, double*ii, 
                                          double *ro, double*io);

   protected: 

      size_t mz_blocksize;     // number of samples in each input
      size_t mz_stepsize;      // sample hop between blocks in input
      size_t mz_channels;      // number of samples in the input

      double*  mz_rtimebuffer; // storage space for transform calculations
      double* mz_rfreqbuffer;  // storage space for transform calculations
      double* mz_ifreqbuffer;  // storage space for transform calculations
      double* mz_window;       // storage for the analysis window
};



#endif // _MZSPECTRUMCLIENT_H_INCLUDED

