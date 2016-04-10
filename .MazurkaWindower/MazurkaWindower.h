//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Fri Jun 16 23:19:44 PDT 2006
// Last Modified: Fri Jun 16 23:19:47 PDT 2006
// Last Modified: Sun May  6 01:48:58 PDT 2007 (upgraded to vamp 1.0)
// Filename:      ...plugin/MazurkaWindower/MazurkaWindower.h
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Interface for windowing audio signals.
//

#ifndef _MAZURKAWINDOWER_H_INCLUDED
#define _MAZURKAWINDOWER_H_INCLUDED

#include "MazurkaTransformer.h"
#include <string>
#include <vector>


class MazurkaWindower {

   public:
                        MazurkaWindower       (void);
                        MazurkaWindower       (int size);
                        MazurkaWindower       (int size, std::string type);
                       ~MazurkaWindower       ();

      double            getWindowSum          (void);
      int               getSize               (void);
      void              setSize               (int size);
      std::string       getWindowType         (void);
      double&           operator[]            (int index);
      MazurkaWindower&  operator=             (MazurkaWindower& aWindow);

      void              windowNonCausal       (MazurkaTransformer& transformer, 
                                               const float* buffer, 
                                               int size);

      static int        getWindowList         (std::vector<std::string>& 
		                              windows);
      static std::string getEnumeratedWindow  (int enumeration);
      int               makeWindow            (std::string type);
      int               makeWindow            (int windowEnum);
      int               makeWindow            (std::string type, int size);
      int               makeWindow            (int windowEnum, int size);
      int               makeWindow            (std::string type, double* data, 
                                              int size);

      // static window generation functions:
      static void       makeRectangularWindow (double* data, int size);
      static void       makeSquareWindow      (double* data, int size);
      static void       makeHannWindow        (double* data, int size);
      static void       makeHanningWindow     (double* data, int size);
      static void       makeBlackmanWindow    (double* data, int size,
                                              double p1 = 0.426590713671539,
                                              double p2 = 0.496560619088564,
                                              double p3 = 0.076848667239896,
                                              double p4 = 0.0);
      static void       makeBlackmanHarris4_92Window  
	                                      (double* data, int size);
      static void       makeTriangularWindow  (double* data, int size);
      static void       makeFejerWindow       (double* data, int size);
      static void       makeBartlettWindow    (double* data, int size);

   protected:

      void              initialize            (int size);
      void              deinitialize          (void);

   private:

      int               dataSize;    // size of the window
      double           *data;        // storage for window shape
      std::string       dataType;    // current window type in storage

};



#endif // _MAZURKAWINDOWER_H_INCLUDED



