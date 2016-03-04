//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Tue Dec 19 20:37:48 PST 2006
// Last Modified: Tue Dec 19 20:37:53 PST 2006
// Filename:      binmaptest.cpp
// Syntax:        ANSI99 C++ 
//
// Description:   Calculate DFT to spectral flux spectrum bin mapping.
//
// Reference:     http://www.ofai.at/~simon.dixon/beatbox
// Reference:     http://en.wikipedia.org/wiki/Spectral_flux
//

#include <iostream>
#include <math.h>

using namespace std;

int makeFreqMap(int* &binmap, int fftsize, float srate);


///////////////////////////////////////////////////////////////////////////

int main(void) {
   int *mapping = NULL;
   int  mapsize = makeFreqMap(mapping, 2048, 44100.0);

   cout << "Size of mapping array: " << mapsize << endl;
   cout << "Mapping contents:\n" << endl; 

   for (int i=0; i<mapsize; i++) {
      cout << i << ":" << mapping[i] << endl;
   }

   delete [] mapping;
   mapping = NULL;

   return 0;
}


///////////////////////////////////////////////////////////////////////////

//////////////////////////////
//
// makeFreqMap -- map from DFT bins to musical note-like binspacings.
//

int makeFreqMap(int* &binmap, int fftsize, float srate) {
   double width  = srate / fftsize;
   double a4freq = 440.0;
   int    a4midi = 69;
   int    mapsize= fftsize/2+1;
   int    xbin   = (int)(2.0/(pow(2.0, 1.0/12.0) - 1.0));
   int    xmidi  = (int)(log(xbin*width/a4freq)/log(2.0)*12 + a4midi + 0.5);
   int    midi;
   int    i;

   delete [] binmap;
   binmap = new int[mapsize];

   for (i=0; i<=xbin; i++) {
      binmap[i] = i;
   }
   for (i=xbin+1; i<mapsize; i++) {
      midi = (int)(log(i*width/a4freq)/log(2.0)*12 + a4midi + 0.5);
      midi = midi > 127 ? 127 : midi;
      binmap[i] = xbin + midi - xmidi;
   }

   return mapsize;
}


