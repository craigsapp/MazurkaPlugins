//
// Programmer:    Madeline Huberth <mhuberth@ccrma.stanford.edu>
// Based on code from: Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Thu Jan  4 09:35:05 PST 2007
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzPitchPower2.h
// URL:           http://sv.mazurka.org.uk/include/MzPitchPower.h
// Documentation: http://sv.mazurka.org.uk/MzPitchPower
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Calculate changes in spectral energy for onset detection.
// 

#ifndef _MZPITCHPOWER2_H_INCLUDED
#define _MZPITCHPOWER2_H_INCLUDED

#include "MazurkaPlugin.h"  // Mazurka plugin interface for Sonic Visualiser
#include "MazurkaTransformer.h"
#include "MazurkaWindower.h"

#include <vector>


class MzPitchPower2 : public MazurkaPlugin {

   public: 

   // plugin interface functions:

                    MzPitchPower2            (float samplerate);
      virtual      ~MzPitchPower2            ();

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
      // size_t     getMinChannelCount      (void) const { return 1; }
      // size_t     getMaxChannelCount      (void) const { return 1; }
     
   // non-interface functions and variables:

      static void   generateMidiNoteList    (std::vector<std::string>& alist,
	                                     int minval = 0, int maxval = 127);
      static void   extractHarmonicPowers   (std::vector<double>& 
                                             harmonic_power, 
                                             std::vector<int>& bins, 
                                             MazurkaTransformer& transformer);
   private: 

      int    mz_harmonics;       // number of harmonics
      float*    mz_realsines;       // the real component of the complex sinusoid up to number of harmonics specified
      float*    mz_imagcosines;
      float mz_samplerate;
    
      // input parameters:
      //
      //    "harmonics";       -- number of harmonics to measure power of
      //    "windowsamples";   -- number of samples in audio window
      //    "stepsamples";     -- number of samples between window starts

};


#endif // _MZPITCHPOWER2_H_INCLUDED

