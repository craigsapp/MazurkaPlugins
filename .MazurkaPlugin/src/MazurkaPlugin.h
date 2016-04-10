//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat May 20 20:33:36 PDT 2006
// Last Modified: Mon May 22 20:31:53 PDT 2006 (automatic database building)
// Last Modified: Sun Jun 18 07:51:22 PDT 2006 (added getParameterString)
// Last Modified: Tue Jul 11 20:29:11 PDT 2006 (added centerTSInBlock)
// Last Modified: Tue Dec 26 22:11:24 PST 2006 (added getParameterDouble)
// Last Modified: Thu Jan  4 23:46:31 PST 2007 (added AUDIODATA typedef)
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MazurkaPlugin.h
// URL:           http://sv.mazurka.org.uk/include/MazurkaPlugin.h
// Documentation: http://sv.mazurka.org.uk/MazurkaPlugin
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Base class for Mazurka Project vamp plugins.
//                Provides an interface for parameter values via
//                the polymorphic functions setParameter() and getParameter().
//                For example usage, see the plugin MzChronogram:
//                   http://sv.mazurka.org.uk/include/MzChronogram.h
//                   http://sv.mazurka.org.uk/src/MzChronogram.cpp
// 
// Usage:         To use this class, inherit the class as in this example:
//
//                  class MyPlugin : public MazurkaPlugin {...};
//    
//                The constructor for MyPlugin should send the sampling rate
//		  to the MazurkaPlugin constructor:
//
//                  MyPlugin::MyPlugin(float srate) : MazurkaPlugin(srate) {
//                     // whatever you need for constructing
//                  }
//
//                Then any parameter you defined in 
//                MyPlugin::getParameterDescriptors() can be read and changed
//                using the getParameter("name") and setParameter("name", value)
//                functions.  The first time you call getParameter() or
//		  setParameter(), the buildParameterDatabase() function will
//		  be called.
//
//	 	  Note that setParameter() will do bounds checking
//                according to the minValue and maxValue for the parameters
//                defined in MyPlugin::getParameterDescriptors().
// 

#ifndef _MAZURKAPLUGIN_H_INCLUDED
#define _MAZURKAPLUGIN_H_INCLUDED

#include "Plugin.h"  // Vamp plugin interface for Sonic Visualiser

#include <string>
#include <vector>
#include <map>

// input audio data type form the Vamp::process function:
// In version 0.9 of vamp:
//typedef float** AUDIODATA;
// In version 1.0 of vamp:
typedef const float* const* AUDIODATA;


typedef struct {
   bool initialized;
   Vamp::PluginBase::ParameterList pdlist; // parameter descriptions
   std::vector<double>       currentValue; // current parameter value
   std::map<std::string,int> indexMap;     // parameter to index number
   std::vector<bool>         hasChanged;   // true if parameter has been
                                           // set after initial assignment
                                           // to the default value.
   std::vector<bool>         isFrozen;     // true if the parameter 
                                           // has been marked as read only.
} ParameterDatabase;


class MazurkaPlugin : public Vamp::Plugin {

   public: 

   // Vamp plugin polymorphic interface functions:
      void          setParameter              (std::string name, 
                                               float value);
      float         getParameter              (std::string name) const;
      float         getParameterFloat         (std::string name) const;
      double        getParameterDouble        (std::string name) const;


   // Mazurka plugin enhanced interface functions:

      int           getParameterInt           (std::string name) const;
      bool          getParameterBool          (std::string name) const;
      std::string   getParameterString        (std::string name) const;
      float         getParameterMin           (std::string name) const;
      float         getParameterMax           (std::string name) const;

      float         getParameterDefault       (std::string name) const;
      int           getParameterDefaultInt    (std::string name) const;
      bool          isParameterAtDefault      (std::string name) const;


      bool          isParameterFrozen         (std::string name) const;
      bool          isValid                   (std::string name) const;

      inline float  getSrate                  (void) const
                                              { return m_inputSampleRate; }
      inline int    getBlockSize              (void) const 
                                              { return mz_blocksize; }
      inline int    getStepSize               (void) const
                                              { return mz_stepsize; }
      inline int    getChannelCount           (void) const
                                              { return mz_channels; }

      Vamp::RealTime centerTimestampInBlock   (Vamp::RealTime timestamp);
      Vamp::RealTime centerTimestampInStep    (Vamp::RealTime timestamp);
      Vamp::RealTime positionTimestampInBlock (Vamp::RealTime timestamp,
		                               double fraction);
      Vamp::RealTime positionTimestampInStep  (Vamp::RealTime timestamp,
		                               double fraction);

   protected:
                    MazurkaPlugin             (float samplerate);

      void          freezeParameter           (std::string name) const;
      void          unfreezeParameter         (std::string name) const;
      void          freezeAllParameters       (void) const;
      void          unfreezeAllParameters     (void) const;

      void          freezeInitialiseData      (void);
      void          unfreezeInitialiseData    (void);

      int           setBlockSize              (int value);
      int           setStepSize               (int value);
      int           setChannelCount           (int value);


      void          buildParameterDatabase    (const ParameterList& list) const;
   private: 

      int           getIndex                  (std::string name) const;

      int objectnumber; // for building database even in const funcs.
      static std::vector<ParameterDatabase> paramdata;

      // variables which will be set in inherited initialise() function:
      int mz_blocksize;             // number of samples in each input frame
      int mz_stepsize;              // sample hop between blocks in input
      int mz_channels;              // number of audio channels in the input
      int mz_initfreeze;            // true if writing to following no allowed:

};


#endif // _MAZURKAPLUGIN_H_INCLUDED

