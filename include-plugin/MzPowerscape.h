//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat Jul  8 00:19:11 PDT 2006 (derived from MzPowerscape)
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzPowerscape.h
// URL:           http://sv.mazurka.org.uk/include/MzPowerscape.h
// Documentation: http://sv.mazurka.org.uk/MzPowerscape
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Display a hierarchical view of the power in a signal.
// 

#ifndef _MZPOWERSCAPE_H_INCLUDED
#define _MZPOWERSCAPE_H_INCLUDED

#include "MazurkaPlugin.h"  // Mazurka plugin interface for Sonic Visualiser
//#include "MazurkaWindower.h"

#include <list>


class MzPowerscape : public MazurkaPlugin {

   public: 

   // plugin interface functions:

                    MzPowerscape            (float samplerate);
      virtual      ~MzPowerscape            ();

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
   
      static double getPowerLevel    (int i, int j, std::vector<double>& power);

   private: 

      int mz_filterforward;             // true if forward filtering
      int mz_filterbackward;            // true if reverse filtering
      int mz_levels;                    // number of rows in display

      // MazurkaWindower mz_window;        // used for weighted averaging
      std::vector<double> surfacepower; // power data for scape display

   // plugin parameters:
   //    "windowsize"         -- size of the analysis window in milliseconds
   //    "hopsize"            -- distance between window start times in ms
   //    "smoothingfactor"    -- gain value for exponential smoothing filter
   //    "filtermethod"       -- which way to filter raw power
   //    "cutoffthreshold"    -- noise floor in dB
   //    "cutoffwidth"        -- transition region around threshold in dB

};


#endif // _MZPOWERSCAPE_H_INCLUDED

