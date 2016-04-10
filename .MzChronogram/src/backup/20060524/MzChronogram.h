//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Tue May  9 05:24:33 PDT 2006
// Last Modified: Sun May 21 00:03:39 PDT 2006 (parameter control added)
// Filename:      MzChronogram.h
// URL:           http://sv.mazurka.org.uk/include/MzChronogram.h
// Documentation: http://sv.mazurka.org.uk/MzChronogram
// Syntax:        ANSI99 C++; vamp 0.9 plugin
//
// Description:   Display audio signal in two dimensions.
// 

#ifndef _MZCHRONOGRAM_H_INCLUDED
#define _MZCHRONOGRAM_H_INCLUDED

#include "MazurkaPlugin.h"  // Mazurka plugin interface for Sonic Visualiser


class MzChronogram : public MazurkaPlugin {

   public: 

   // plugin interface functions:

                    MzChronogram            (float samplerate);
      virtual      ~MzChronogram            ();

      // required polymorphic functions inherited from PluginBase:
      std::string   getName                 (void) const;
      std::string   getMaker                (void) const;
      std::string   getCopyright            (void) const;
      std::string   getDescription          (void) const;
      int           getPluginVersion        (void) const;

      // optional parameter interface functions
      ParameterList getParameterDescriptors (void) const;

      // required polymorphic functions inherited from Plugin:
      InputDomain   getInputDomain          (void) const;
      OutputList    getOutputDescriptors    (void) const;
      bool          initialise              (size_t channels, 
                                             size_t stepsize, 
                                             size_t blocksize);
      FeatureSet    process                 (float **inputbufs, 
                                             Vamp::RealTime timestamp);
      FeatureSet    getRemainingFeatures    (void);
      void          reset                   (void);

      // optional polymorphic functions from Plugin:
      size_t        getPreferredStepSize    (void) const;
      size_t        getPreferredBlockSize   (void) const;
      // size_t     getMinChannelCount      (void) const { return 1; }
      size_t        getMaxChannelCount      (void) const { return 999; }


   // non-interface functions and variables:
   
      inline int    getChannels      (void) const { return mz_channels;  } 
      inline int    getStepSize      (void) const { return mz_stepsize;  } 
      inline int    getBlockSize     (void) const { return mz_blocksize; } 
      inline float  getSrate         (void) const { return m_inputSampleRate; } 

   private: 

      int mz_blocksize;        // number of samples in each input frame
      int mz_stepsize;         // sample hop between blocks in input
      int mz_channels;         // number of audio channels in the input

      // input parameters:
      // 
      //    "verticalperiod"   -- number of samples on vertical axis

};


#endif // _MZCHRONOGRAM_H_INCLUDED

