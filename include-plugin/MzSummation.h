//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Wed Mar  7 23:54:24 PST 2007
// Last Modified: Wed Mar  7 23:54:30 PST 2007
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MzSummation.h
// URL:           http://sv.mazurka.org.uk/include/MzSummation.h
// Documentation: http://sv.mazurka.org.uk/MzSummation
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Simple plugin to test SV behavior in multiple OSes.
// 

#ifndef _MZSUMMATION_H_INCLUDED
#define _MZSUMMATION_H_INCLUDED

#include "MazurkaPlugin.h"  // Mazurka plugin interface for Sonic Visualiser
#include "MazurkaWindower.h"

#include <list>


class MzSummation : public MazurkaPlugin {

   public: 

   // plugin interface functions:

                    MzSummation             (float samplerate);
      virtual      ~MzSummation             ();

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
   
   private: 
      MazurkaTransformer mz_transformer;  // used for transformations

};


#endif // _MZSUMMATION_H_INCLUDED

