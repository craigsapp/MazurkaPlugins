//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Mon May  8 18:10:56 PDT 2006
// Last Modified: Fri Mar  4 04:49:58 PST 2016
// Filename:      MazurkaPlugins/src/plugins.cpp
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Automatically generated file based on the 
//                vamp sdk file examples/plugins.cpp
//

#include "vamp.h"
#include "PluginAdapter.h"

#include "MzPowerCurve.h"

static Vamp::PluginAdapter<MzPowerCurve> mzPowerCurveAdapter;

// Pre vamp 1.0 interface:
//const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int index) {

const VampPluginDescriptor *vampGetPluginDescriptor(
      unsigned int vampApiVersion, unsigned int index) {
   if (vampApiVersion < 1) {
       return 0;
   }

   const char* setinfo = "@@VampPluginSet@" __DATE__ "@MzPowerCurve@@";
   if (setinfo[0] != '@') { 
      std::cerr << "This is a dummy statment: " << setinfo << std::endl;
   }

    switch (index) {
       case  0: return mzPowerCurveAdapter.getDescriptor();
       default: return 0;
    }
}
