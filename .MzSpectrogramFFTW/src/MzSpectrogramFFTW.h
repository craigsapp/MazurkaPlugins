//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun Jun 11 21:02:20 PDT 2006
// Last Modified: Fri Jun 23 01:39:57 PDT 2006 (subclassed to MazurkaPlugin)
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzSpectrogramFFTW.h
// URL:           http://sv.mazurka.org.uk/include/MzSpectrogramFFTW.h
// Documentation: http://sv.mazurka.org.uk/MzSpectrogramFFTW
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Demonstration of how to create spectral data from time data
//                supplied by the host application using the FFTW library
//                for Fourier Transforms.
// 

#ifndef _MZSPECTROGRAMFFTW_H_INCLUDED
#define _MZSPECTROGRAMFFTW_H_INCLUDED

#include "MazurkaPlugin.h"      // Mazurka plugin interface for Sonic Visualiser
#include "MazurkaTransformer.h" // Mazurka interface to FFTW 


class MzSpectrogramFFTW : public MazurkaPlugin {

   public: 

   // plugin interface functions:

                    MzSpectrogramFFTW     (float samplerate);
      virtual      ~MzSpectrogramFFTW     ();

      // required polymorphic functions inherited from PluginBase:
      std::string   getIdentifier         (void) const;
      std::string   getName               (void) const;
      std::string   getDescription        (void) const;
      std::string   getMaker              (void) const;
      std::string   getCopyright          (void) const;
      int           getPluginVersion      (void) const;

      // optional parameter interface functions:
      ParameterList getParameterDescriptors (void) const;


      // required polymorphic functions inherited from Plugin:
      InputDomain   getInputDomain        (void) const;
      OutputList    getOutputDescriptors  (void) const;
      bool          initialise            (size_t channels, 
                                           size_t stepsize, 
                                           size_t blocksize);
      FeatureSet    process               (AUDIODATA inputbufs, 
                                           Vamp::RealTime timestamp);
      FeatureSet    getRemainingFeatures  (void);
      void          reset                 (void);

      // optional polymorphic functions from Plugin:
      // size_t     getPreferredStepSize  (void) const { return 0; }
      // size_t     getPreferredBlockSize (void) const { return 0; }
      // size_t     getMinChannelCount    (void) const { return 1; }
      // size_t     getMaxChannelCount    (void) const { return 1; }


   // non-interface functions and variables:
   
      static void   makeHannWindow        (double* output, int blocksize);
      static void   windowSignal          (MazurkaTransformer& transformer,
		                           double* window, const float* input);
   private: 

      int                mz_minbin;       // minimum spectral bin to display
      int                mz_maxbin;       // maximum spectral bin to display
      double            *mz_wind_buff;    // storage for the window
      MazurkaTransformer mz_transformer;  // FFTW interface
      
};


#endif // _MZSPECTROGRAMFFTW_H_INCLUDED

