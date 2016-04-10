//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat May 13 12:16:45 PDT 2006
// Last Modified: Tue May 16 18:58:27 PDT 2006
// Filename:      MzPowerCurve.h
// URL:           http://sv.mazurka.org.uk/MzPowerCurve/src/MzPowerCurve.h
// Documentation: http://sv.mazurka.org.uk/MzPowerCurve
// Syntax:        ANSI C++; vamp plugin
//
// Description:   Calculate the power of an audio signal.
// 

#ifndef _MZPOWERCURVE_H_INCLUDED
#define _MZPOWERCURVE_H_INCLUDED

#include "Plugin.h"  // Vamp plugin interface for Sonic Visualiser

#include <list>
using std::list;


class MzPowerCurve : public Vamp::Plugin {

   public: 

   // plugin interface functions:

                    MzPowerCurve          (float samplerate);
      virtual      ~MzPowerCurve          ();

      // required polymorphic functions inherited from PluginBase:
      std::string   getName               (void) const;
      std::string   getMaker              (void) const;
      std::string   getCopyright          (void) const;
      std::string   getDescription        (void) const;
      int           getPluginVersion      (void) const;

      // required polymorphic functions inherited from Plugin:
      InputDomain   getInputDomain        (void) const;
      OutputList    getOutputDescriptors  (void) const;
      bool          initialise            (size_t channels, size_t stepsize, 
                                           size_t blocksize);
      FeatureSet    process               (float **inputbufs, 
                                           Vamp::RealTime timestamp);
      FeatureSet    getRemainingFeatures  (void);
      void          reset                 (void);

      // optional polymorphic functions from Plugin:
      
      // hardwired to 10 ms at 44100 sampling rate for now...
      size_t        getPreferredStepSize  (void) const { return mzp_windowstep;}
      size_t        getPreferredBlockSize (void) const { return mzp_windowsize;}

      // size_t     getMinChannelCount    (void) const { return 1; }
      // size_t     getMaxChannelCount    (void) const { return 1; }

   // non-interface functions and variables:
   
      inline int    getChannels           (void) const { return mz_channels;  } 
      inline int    getStepSize           (void) const { return mz_stepsize;  } 
      inline int    getBlockSize          (void) const { return mz_blocksize; } 
      inline float  getSrate              (void) const 
                                                 { return m_inputSampleRate; }
       
      static double getPower              (float* data, int blocksize);

   private: 

      int mz_blocksize;             // number of samples in each input frame
      int mz_stepsize;              // sample hop between blocks in input
      int mz_channels;              // number of audio channels in the input
      std::list<double> rawpower;   // power data for non-causal calculations

      // Vamp controllable parameters:
      int    mzp_windowsize;        // block size of the power window
      int    mzp_windowstep;        // step size of the power window
      double mzp_cutoff;            // cutoff threshold for scaled power slope
      double mzp_width;             // threshold width for scaled power slope

};



#endif // _MZPOWERCURVE_H_INCLUDED

