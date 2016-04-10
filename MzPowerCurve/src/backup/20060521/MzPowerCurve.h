//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat May 13 12:16:45 PDT 2006
// Last Modified: Sat May 20 15:50:06 PDT 2006 (parameters control added)
// Filename:      MzPowerCurve.h
// URL:           http://sv.mazurka.org.uk/MzPowerCurve/src/MzPowerCurve.h
// Documentation: http://sv.mazurka.org.uk/MzPowerCurve
// Syntax:        ANSI99 C++; vamp 0.9 plugin
//
// Description:   Calculate the power of an audio signal.
// 

#ifndef _MZPOWERCURVE_H_INCLUDED
#define _MZPOWERCURVE_H_INCLUDED

#include "MazurkaPlugin.h"  // Mazurka plugin interface for Sonic Visualiser

#include <list>
using std::list;


class MzPowerCurve : public MazurkaPlugin {

   public: 

   // plugin interface functions:

                    MzPowerCurve            (float samplerate);
      virtual      ~MzPowerCurve            ();

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
      bool          initialise              (size_t channels, size_t stepsize, 
                                             size_t blocksize);
      FeatureSet    process                 (float **inputbufs, 
                                             Vamp::RealTime timestamp);
      FeatureSet    getRemainingFeatures    (void);
      void          reset                   (void);

      // optional polymorphic functions from Plugin:
      size_t        getPreferredStepSize    (void) const;
      size_t        getPreferredBlockSize   (void) const;
      // size_t     getMinChannelCount      (void) const { return 1; }
      // size_t     getMaxChannelCount      (void) const { return 1; }

   // non-interface functions and variables:
   
      inline int    getChannels      (void) const { return mz_channels;       }
      inline int    getStepSize      (void) const { return mz_stepsize;       }
      inline int    getBlockSize     (void) const { return mz_blocksize;      }
      inline float  getSrate         (void) const { return m_inputSampleRate; }
       
      static double getPower         (float* data, int blocksize);

   private: 

      int mz_blocksize;             // number of samples in each input frame
      int mz_stepsize;              // sample hop between blocks in input
      int mz_channels;              // number of audio channels in the input
      std::list<double> rawpower;   // power data for non-causal calculations

   /* plugin parameters:
    *    "windowsize"      -- size of the analysis window in milliseconds.
    *    "hopsize"         -- distance between window start times in ms.
    *    "cutoffthreshold" -- noise floor in dB.
    *    "cutoffwidth"     -- transition region around threshold in dB.
    */

};


#endif // _MZPOWERCURVE_H_INCLUDED

