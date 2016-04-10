//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sat May 20 20:33:36 PDT 2006
// Last Modified: Mon May 22 20:31:53 PDT 2006 (automatic database building)
// Last Modified: Sun Jun 18 07:51:22 PDT 2006 (added getParameterString)
// Last Modified: Wed Jun 21 07:18:35 PDT 2006 (added some variables)
// Last Modified: Tue Dec 26 22:11:24 PST 2006 (added getParameterDouble)
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      MazurkaPlugin.cpp
// URL:           http://sv.mazurka.org.uk/src/MazurkaPlugin.cpp
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

#include "MazurkaPlugin.h"  


// paramdata is storage space for parameter data.  It is declared static
// to allow for automatic parameter database building inside of virtual
// const functions inherited from Vamp::PluginBase.  If the parameter
// database functionality is implemented in the class's parent class,
// then the automatic building of the database can be stored inside of
// the class rather than in this static structure (and then checking to
// see if the database has been built would not be needed in all parameter
// access functions).

std::vector<ParameterDatabase> MazurkaPlugin::paramdata;




//////////////////////////////////////////////////////////////////////////
//
// vamp polymorphic functions inherited from the PluginBase class
//

//////////////////////////////
//
// MazurkaPlugin::setParameter -- set the value of a parameter.
//    If the parameter does not exist, then ignore the request.
//    Also, this function enforces the min..max range as defined
//    in ParameterDescriptors() for each parameter.  The 
//    quantization settings of the Parameter Descriptors are not 
//    enforced.  If the parameter database has not yet been built
//    (such as the first time setParameter() is called) it will
//    be built.  Subsequent calls to setParameter() will not rebuild
//    the parameter database.
//

void MazurkaPlugin::setParameter(std::string name, float value) {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return;
   }

   // don't do anything if the parameter has been frozen
   if (pd.isFrozen[index] == true) {
      return;
   }

   // check to make sure that the value is in range and then store
   if (value < pd.pdlist[index].minValue) { 
      value = pd.pdlist[index].minValue; 
   } else if (value > pd.pdlist[index].maxValue) { 
      value = pd.pdlist[index].maxValue; 
   }
   pd.currentValue[index] = value;

std::cerr << "Setting " << name << " to " << value << std::endl;

   if (value != pd.pdlist[index].defaultValue) {
      // hasChanged remains false until the first time that
      // the parameter is set to anything other than the default value.
      // This is to prevent fake changes, where the host application
      // gratutitously sets all parameters to their defaults (which
      // the setParameter() function already knows about).
      pd.hasChanged[index]   = true;
   }
}



//////////////////////////////
//
// MazurkaPlugin::getParameter -- return a parameter's current value.
//    If the parameter database has not yet been built (such as the 
//    first time getParameter() is called) it will be built.  Subsequent 
//    calls to getParameter() will not rebuild the parameter database.

float MazurkaPlugin::getParameter(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }
   int index = getIndex(name);
   if (index < 0) {
      // unknown parameter, just return 0.0
      return 0.0;
   } else {
      return float(pd.currentValue[index]);
   }
}

float MazurkaPlugin::getParameterFloat(std::string name) const {
   return getParameter(name);
}

double MazurkaPlugin::getParameterDouble(std::string name) const {
   return (double)getParameter(name);
}

//////////////////////////////////////////////////////////////////////////
//
// non-vamp functions 
//

//////////////////////////////
//
// MazurkaPlugin::freezeInitialiseData -- prevent changing
//   data with setChannelCount(), setBlockSize(), and
//   setStepSize().
//

void MazurkaPlugin::freezeInitialiseData(void) { 
   mz_initfreeze = 0;
}



//////////////////////////////
//
// MazurkaPlugin::unfreezeInitialiseData -- allow changing
//   data with setChannelCount(), setBlockSize(), and
//   setStepSize().
//

void MazurkaPlugin::unfreezeInitialiseData(void) { 
   mz_initfreeze = 0;
}



//////////////////////////////
//
// MazurkaPlugin::setBlockSize -- mz_blocksize stores the number
//    of samples in each frame of data sent to the plugin with
//    the process() function.
//

int MazurkaPlugin::setBlockSize(int value) { 
   if (mz_initfreeze == 0) {
      mz_blocksize = value;
      return 1;
   } else {
      return 0;
   }
}



//////////////////////////////
//
// MazurkaPlugin::centerTimestampInBlock --
//

Vamp::RealTime MazurkaPlugin::centerTimestampInBlock(Vamp::RealTime timestamp) {
   return positionTimestampInBlock(timestamp, 0.5);
}



//////////////////////////////
//
// MazurkaPlugin::centerTimestampInStep --
//

Vamp::RealTime MazurkaPlugin::centerTimestampInStep(Vamp::RealTime timestamp) {
   return positionTimestampInStep(timestamp, 0.5);
}



//////////////////////////////
//
// MazurkaPlugin::positionTimestampInBlock --  move the timestamp
//   to a relative position in the block, with 0.0 being the start,
//   and 1.0 being the end of the block.
//

Vamp::RealTime 
MazurkaPlugin::positionTimestampInBlock(Vamp::RealTime timestamp, 
      double fraction) {
   return timestamp + 
          Vamp::RealTime::fromSeconds(fraction * getBlockSize() / getSrate());
}



//////////////////////////////
//
// MazurkaPlugin::positionTimestampInStep --  move the timestamp
//   to a relative position in the step size, with 0.0 being the start,
//   and 1.0 being the end of the step duration.
//

Vamp::RealTime 
MazurkaPlugin::positionTimestampInStep(Vamp::RealTime timestamp,
      double fraction) {
   return timestamp + 
          Vamp::RealTime::fromSeconds(fraction * getStepSize() / getSrate());
}



//////////////////////////////
//
// MazurkaPlugin::setStepSize -- mz_stepsize is the number of
//   samples between frame sample starting points.
//

int MazurkaPlugin::setStepSize(int value) { 
   if (mz_initfreeze == 0) {
      mz_stepsize = value;
      return 1;
   } else {
      return 0;
   }
}



//////////////////////////////
//
// MazurkaPlugin::setChannelCount -- mz_channels is the
//   number of channels of data sent to the plugin using
//   the process() function.
//

int MazurkaPlugin::setChannelCount(int value) { 
   if (mz_initfreeze == 0) {
      mz_channels = value;
      return 1;
   } else {
      return 0;
   }
}



//////////////////////////////
//
// MazurkaPlugin::getParameterInt -- return a parameter's current value
//    rounded to the nearest integer.  If the parameter database has not 
//    yet been built (such as the first time getParameter*() is called) 
//    it will be built.  Subsequent calls to getParameterInt() will not 
//    rebuild the parameter database.
//

int MazurkaPlugin::getParameterInt(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return 0;
   } 

   double value = pd.currentValue[index];

   // ceiling for negative numbers; floor for positive numbers:
   if (value < 0.0) { return int(value - 0.5); }
   else             { return int(value + 0.5); }
}



//////////////////////////////
//
// MazurkaPlugin::getParameterBool -- return true if a parameter has
//    been set after it initial assignment to the default value.
//    This function is useful for choosing between two inputs 
//    according to the one which has been changed.  For example, 
//    one parameter could set a window size in samples, another
//    parameter could set the window size in seconds.  The plugin
//    can choose which one to use depending if the user specified
//    one or the other parameter.
//
//    If the parameter database has not yet been built (such as the first
//    time getParameter*() is called) it will be built.  Subsequent
//    calls to getParameterBool() will not rebuild the parameter database.
//

bool MazurkaPlugin::getParameterBool(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return false;
   } 

   return pd.hasChanged[index];
}



//////////////////////////////
//
// MazurkaPlugin::getParameterString -- return the string value
//   which is associated with the current (quantized) value of 
//   of the given parameter.
//

std::string MazurkaPlugin::getParameterString(std::string name) const {
	   
   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return "";
   } 

   if (pd.pdlist[index].valueNames.size() <= 0) {
      // there are no strings associated with the parameter
      return "";
   }

   if (pd.pdlist[index].isQuantized == false) {
       // strange, there should only be strings when the parameter
       // is quantized.
       return "";
   }

   if (pd.pdlist[index].quantizeStep <= 0.0) {
       // strange, the quantization step should be positive number.
       return "";
   }

   float& minval = pd.pdlist[index].minValue;
   float& maxval = pd.pdlist[index].maxValue;
   double& curval = pd.currentValue[index];
   int stringcount = pd.pdlist[index].valueNames.size();

   int stringindex = int(stringcount * (curval-minval)/(maxval-minval+1) + 0.5);

   if (stringindex > 0 && stringindex < stringcount) {
      return pd.pdlist[index].valueNames[stringindex];      
   }

   return "";
}



//////////////////////////////
//
// MazurkaPlugin::getParameterMin -- return the minimum allowed
//    value for the specified parameter.   Returns 0.0 if the
//    parameter is not defined.  If the parameter database has not 
//    yet been built (such as the first time getParameter*() is called) 
//    it will be built.  Subsequent calls to getParameterMin() will not 
//    rebuild the parameter database.
//

float MazurkaPlugin::getParameterMin(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return 0.0;
   } 

   return pd.pdlist[index].minValue;
}



//////////////////////////////
//
// MazurkaPlugin::getParameterMax -- return the maximum allowed
//    value for the specified parameter.   Returns 0.0 if the
//    parameter is not defined.  If the parameter database has not 
//    yet been built (such as the first time getParameter*() is called) 
//    it will be built.  Subsequent calls to getParameterMin() will not 
//    rebuild the parameter database.
//

float MazurkaPlugin::getParameterMax(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return 0.0;
   } 

   return pd.pdlist[index].maxValue;
}



//////////////////////////////
//
// MazurkaPlugin::getParameterDefault -- return the default value
//    for a parameter.  If the parameter database has not yet been built 
//    (such as the first time getParameter() is called) it will be built.  
//    Subsequent calls to getParameterDefault() will not rebuild the 
//    parameter database.
//

float MazurkaPlugin::getParameterDefault(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }
   int index = getIndex(name);
   if (index < 0) {
      // unknown parameter, just return 0.0
      return 0.0;
   } else {
      return float(pd.pdlist[index].defaultValue);
   }
}



//////////////////////////////
//
// MazurkaPlugin::getParameterDefaultInt -- return a parameter's default 
//    value, rounded to the nearest integer.  If the parameter database has 
//    not yet been built (such as the first time getParameter() is called) 
//    it will be built.  Subsequent calls to getParameterDefaultInt() will 
//    not rebuild the parameter database.
//

int MazurkaPlugin::getParameterDefaultInt(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return 0;
   } 

   double value = pd.pdlist[index].defaultValue;

   // ceiling for negative numbers; floor for positive numbers:
   if (value < 0.0) { return int(value - 0.5); }
   else             { return int(value + 0.5); }
}



//////////////////////////////
//
// MazurkaPlugin::isParameterAtDefault -- return true if a 
//    parameter's current value is equal to the default value.
//    If the parameter database has not yet been built (such as the first
//    time getParameter*() is called) it will be built.  Subsequent
//    calls to getParameterBool() will not rebuild the parameter database.
//

bool MazurkaPlugin::isParameterAtDefault(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return 0;
   } 

   return pd.pdlist[index].defaultValue == pd.currentValue[index];
}



//////////////////////////////
//
// MazurkaPlugin::isParameterFrozen -- returns true if the
//     specified parameter is read-only.  Will also return true for
//     invalid parameters, since they cannot be set to anything
//     other than 0.0.  The parameter database will be
//     built if no parameter access has yet been done.
//

bool MazurkaPlugin::isParameterFrozen(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return true;
   } 

   return pd.isFrozen[index];
}



//////////////////////////////
//
// MazurkaPlugin::isValid -- returns true if the specified 
//    parameter name exists in the database.  The parameter
//    database will be built if no parameter access has yet
//    been done.
//

bool MazurkaPlugin::isValid(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return false;
   } 

   return true;
}


////////////////////////////////////////////////////////////
//
// private/protected functions
//

//////////////////////////////
//
// MazurkaPlugin::MazurkaPlugin -- constructor
//

MazurkaPlugin::MazurkaPlugin(float samplerate) : Plugin(samplerate) {
   mz_blocksize  = 0;
   mz_stepsize   = 0;
   mz_channels   = 0;
   mz_initfreeze = 0;

   ParameterDatabase pd;
   MazurkaPlugin::paramdata.push_back(pd);
   objectnumber = MazurkaPlugin::paramdata.size() - 1;

   MazurkaPlugin::paramdata[objectnumber].initialized = false;
   // maybe need to decide how to clean contents of 
   // parameter data on destruction of object so that too much
   // memory isn't wasted...
}



//////////////////////////////
//
// MazurkaPlugin::freezeParameter -- Mark the specified
//     parameter as read-only so that the setParameter()
//     function will ignore any requests to change the
//     parameter value.  The parameter database will be
//     built if no parameter access has yet been done.
//

void MazurkaPlugin::freezeParameter(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return;
   } 

   pd.isFrozen[index] = true;
}



//////////////////////////////
//
// MazurkaPlugin::unfreezeParameter -- Mark the specified
//     parameter as read/write so that the setParameter() 
//     function can change the value of the parameter.
//     The parameter database will be built if no parameter
//     access has yet been done.
//

void MazurkaPlugin::unfreezeParameter(std::string name) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int index = getIndex(name);
   if (index < 0) {
      return;
   } 

   pd.isFrozen[index] = false;
}



//////////////////////////////
//
// MazurkaPlugin::freezeAllParameters -- call freezeParameter()
//     on all parameters in the database.  The parameter database
//     will be built if no parameter access has yet been done.
//

void MazurkaPlugin::freezeAllParameters(void) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int size = pd.isFrozen.size();
   for (int i=0; i<size; i++) {
      pd.isFrozen[i] = true;
   }
}



//////////////////////////////
//
// MazurkaPlugin::unfreezeAllParameters -- call unfreezeParameter()
//     on all parameters in the database.  
//

void MazurkaPlugin::unfreezeAllParameters(void) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (!pd.initialized) {
      buildParameterDatabase(getParameterDescriptors());
   }

   int size = pd.isFrozen.size();
   for (int i=0; i<size; i++) {
      pd.isFrozen[i] = false;
   }
}



//////////////////////////////
//
// MazurkaPlugin::buildParameterDatabase -- Store the list of
//     parameters understood by this plugin.  You should not
//     define more than one parameter with the same name since
//     one of them will disappear (could be any one of the duplicates
//     but most likely the earlier defined ones will be lost).
//

void MazurkaPlugin::buildParameterDatabase(const ParameterList& list) const {


   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   if (pd.initialized) {
      // Sorry, you can only initialize once
      return;
   } else{
      pd.initialized = true;
   }
	  
   pd.pdlist = list;
   pd.indexMap.clear();
   pd.currentValue.clear();
   pd.hasChanged.clear();
   pd.isFrozen.clear();

   typedef std::pair<std::string,int> mapentry;
   
   // store the default values in currentValue array, checking bounds
   // defined in minValue and maxValue for each parameter:
   double value;
   int size = pd.pdlist.size();
   for (int i=0; i<size; i++) {
      value = pd.pdlist[i].defaultValue;
      if      (value < pd.pdlist[i].minValue) { value = pd.pdlist[i].minValue; }
      else if (value > pd.pdlist[i].maxValue) { value = pd.pdlist[i].maxValue; }

      // change the default value if it is out of range:
      pd.pdlist[i].defaultValue = value;

      pd.currentValue.push_back(value);
      pd.hasChanged.push_back(false);
      pd.isFrozen.push_back(false);
      pd.indexMap.insert(mapentry(pd.pdlist[i].identifier, i));
   }

}



//////////////////////////////
//
// MazurkaPlugin::getIndex -- return the index of the associated plugin
//     name.  Returns -1 if there is no parameter by that name.
//

int MazurkaPlugin::getIndex(std::string identifier) const {

   ParameterDatabase& pd = MazurkaPlugin::paramdata[objectnumber];

   std::map<std::string,int>::const_iterator citer;

   citer = pd.indexMap.find(identifier);

   if (citer != pd.indexMap.end()) {
      return citer->second;
   } else {
      return -1;
   }
}


