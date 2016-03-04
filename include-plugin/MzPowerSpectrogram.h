//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat Oct  9 22:36:57 PDT 2010 (copied from MzNevermore)
// Last Modified: Sat Oct  9 22:37:09 PDT 2010
// Filename:      MzPowerSpectrogram.h
// URL:           http://sv.mazurka.org.uk/include/MzPowerSpectrogram.h
// Documentation: http://sv.mazurka.org.uk/MzPowerSpectrogram
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   DFT spectrogram with independent window and transform size.
// 

#ifndef _MZNEVERMORE_H_INCLUDED
#define _MZNEVERMORE_H_INCLUDED

#include "MazurkaPlugin.h"  // Mazurka plugin interface for Sonic Visualiser
#include "MazurkaTransformer.h"
#include "MazurkaWindower.h"


class MzPowerSpectrogram : public MazurkaPlugin {

   public: 

   // plugin interface functions:

                    MzPowerSpectrogram      (float samplerate);
      virtual      ~MzPowerSpectrogram      ();

      // required polymorphic functions inherited from PluginBase:
      std::string   getIdentifier           (void) const;
      std::string   getName                 (void) const;
      std::string   getDescription          (void) const;
      std::string   getMaker                (void) const;
      std::string   getCopyright            (void) const;
      int           getPluginVersion        (void) const;

      // optional parameter interface functions
      ParameterList getParameterDescriptors (void) const;

      // required polymorphic functions inherited from Plugin:
      InputDomain   getInputDomain          (void) const;
      OutputList    getOutputDescriptors    (void) const;
      bool          initialise              (size_t channels, 
                                             size_t stepsize, 
                                             size_t blocksize);
      FeatureSet    process                 (AUDIODATA inputbufs, 
                                             Vamp::RealTime timestamp);
      FeatureSet    getRemainingFeatures    (void);
      void          reset                   (void);

      // optional polymorphic functions from Plugin:
      size_t        getPreferredStepSize    (void) const;
      size_t        getPreferredBlockSize   (void) const;
      size_t        getMinChannelCount      (void) const { return 1; }
      size_t        getMaxChannelCount      (void) const { return 1; }

   // non-interface functions and variables:

   private: 

      int    mz_transformsize; // DFT transform size
      int    mz_minbin;        // minimum bin to display
      int    mz_maxbin;        // maximum bin to display
      int    mz_compress;      // for compressing the magnigude range
      int    mz_scale;         // for the vertical scale of freq. axis
       
      // long window transformer
      MazurkaTransformer mz_transformerB;  // interface FFTW Fourier transforms
      MazurkaWindower    mz_windowerB;     // interface for windowsing signals

      // short window transformer
      MazurkaTransformer mz_transformerS;  // interface FFTW Fourier transforms
      MazurkaWindower    mz_windowerS;     // interface for windowsing signals

      // input parameters:
      // 
      //    "windowsamples"    -- number of samples in long audio window
      //    "transformsamples" -- number of samples in transform
      //    "stepsamples"      -- number of samples between analysis windows
      //    "minbin"           -- lowest transform bin to display
      //    "maxbin"           -- highest transform bin to display
      //    "scale"            -- linear or logarithmic scaling of the freqs.

};


#endif // _MZNEVERMORE_H_INCLUDED

