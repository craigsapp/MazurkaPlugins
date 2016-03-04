//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun Jun 11 05:13:19 PDT 2006
// Last Modified: Sun Jun 11 05:13:19 PDT 2006
// Last Modified: Tue Apr 15 00:38:51 PDT 2008 (added inverse transform)
// Filename:      ...plugin/MazurkaTransformer/MazurkaTransformer.cpp
// Web:           http://sv.mazurka.org.uk/src/MazurkaTransformer.cpp
// Syntax:        ANSI99 C++ 
//
// Description:   Interface to FFTW for doing Discrete Fourier Transforms 
//                on real signals.
//


#include "MazurkaTransformer.h"
#include <math.h>


//////////////////////////////
//
// MazurkaTransformer::MazurkaTransformer -- constructor
//

MazurkaTransformer::MazurkaTransformer(void) {
   fftSignalN      = 0;
   fftSignalNd2    = 0;
   fftSpectrumN    = 0;
   fftSignal       = NULL;
   fftSpectrum     = NULL;
   fftPlan         = NULL;
   fftPlan_inverse = NULL;
}


MazurkaTransformer::MazurkaTransformer(int size) {
   fftSignalN      = 0;
   fftSignalNd2    = 0;
   fftSpectrumN    = 0;
   fftSignal       = NULL;
   fftSpectrum     = NULL;
   fftPlan         = NULL;
   fftPlan_inverse = NULL;

   initialize(size);
}



//////////////////////////////
//
// MazurkaTransformer::~MazurkaTransformer -- destructor
//

MazurkaTransformer::~MazurkaTransformer() {
   deinitialize();
}



//////////////////////////////
//
// MazurkaTransformer::doTransform -- convert the real signal 
//      into a complex spectrum.
//

int MazurkaTransformer::doTransform(void) {
   if (fftPlan != NULL) {
      fftw_execute(fftPlan);
      return 1;  // fftPlan is initialized
   } else {
      return 0;  // fftPlan is not initialized
   }
}



//////////////////////////////
//
// MazurkaTransformer::doInverseTransform -- convert complex spectrum
//     into a real signal.
//

int MazurkaTransformer::doInverseTransform(void) {
   if (fftPlan_inverse != NULL) {
      fftw_execute(fftPlan_inverse);
      for (int i=0; i<fftSignalN; i++) {
         fftSignal[i] /= fftSignalN; 
      }
      return 1;  // fftPlan_inverse is initialized
   } else {
      return 0;  // fftPlan_inverse is not initialized
   }
}



//////////////////////////////
//
// MazurkaTransformer::getSize --
//

int MazurkaTransformer::getSize(void) {
   return fftSignalN;
}



//////////////////////////////
//
// MazurkaTransformer::setSize --
//

void MazurkaTransformer::setSize(int size) {
   initialize(size);
}



//////////////////////////////
//
// MazurkaTransformer::signalNonCausal -- index the array with
//     time 0 at the middle index.
//

double& MazurkaTransformer::signalNonCausal(int index) {

   index += fftSignalNd2;
   if (index >= fftSignalN) {
      index -= fftSignalN;
   }

   static double dummy;
   if (index < 0 || index >= fftSignalN) {
      return dummy;
   }

   return fftSignal[index];
}



//////////////////////////////
//
// MazurkaTransformer::zeroSignal -- set signal to all zeros.
//

void MazurkaTransformer::zeroSignal(void) {
   for (int i=0; i<fftSignalN; i++) {
      fftSignal[i] = 0.0;
   }
}



//////////////////////////////
//
// MazurkaTransformer::signalCausal -- index the array with
//     time 0 at the first index.
//

double& MazurkaTransformer::signalCausal(int index) {

   static double dummy = 1.0;
   if (index < 0 || index >= fftSignalN) {
      return dummy;
   }

   return fftSignal[index];
}



//////////////////////////////
//
// MazurkaTransformer::operator[] --
//

double& MazurkaTransformer::operator[](int index) {
   return  fftSignal[index];
}


//////////////////////////////
//
// MazurkaTransformer::getSpectrum -- Gives the spectrum (valid after
//     calling doTransformation() on a given signal).
//

mz_complex MazurkaTransformer::getSpectrum(int index) {
   static mz_complex dummy;
   if (index < 0) {
      index += fftSignalN;
      if (index < 0) {
         return dummy;
      }
   }
   if (index < fftSpectrumN) {
      dummy.re = fftSpectrum[index][0];
      dummy.im = fftSpectrum[index][1];
      return dummy;
   }
   if (index < fftSignalN) {
      dummy.re = fftSpectrum[fftSignalN - index][0];
      dummy.im = -fftSpectrum[fftSignalN - index][1];
   }
   return dummy;
}



//////////////////////////////
//
// MazurkaTransformer::getSpectrumMagnitude -- return the absolute
//     value of the complex spectrum value.
//

double MazurkaTransformer::getSpectrumMagnitude(int index) {
   mz_complex num = getSpectrum(index);
   return sqrt(num.re * num.re + num.im * num.im);
}

//////////////////////////////
//
// MazurkaTransformer::getSpectrumSquared -- return the squared
//     value of the complex spectrum 
//

double MazurkaTransformer::getSpectrumSquared(int index) {
   mz_complex num = getSpectrum(index);
   return num.re * num.re + num.im * num.im;
}



//////////////////////////////
//
// MazurkaTransformer::getSpectrumMagnitudeDb -- return the magnitude
//     of th complex spectrum value in decibels.
//
//     default value: reference = 1.0;
//

double MazurkaTransformer::getSpectrumMagnitudeDb(int index, double reference) {
   double num = getSpectrumSquared(index);
   #define ZEROLOG -120.0

   if (num <= 0.0) {
      return ZEROLOG;
   } else if (reference == 1.0) {
      return 10 * log10(num);
   } else {
      return 10 * log10(num/(reference * reference));
   }
}


///////////////////////////////////////////////////////////////////////////
//
// protected functions
//

//////////////////////////////
//
// MazurkaTransformer::initialize -- Prepare for doing transforms
//    at a given signal length.
//

int MazurkaTransformer::initialize(int size) {
   if (fftSignalN == size) {
      return 1;
   }

   deinitialize();

   if (size <= 0) {
      return 0;
   }

   fftSignalN   = size;
   fftSignalNd2 = size/2;
   fftSpectrumN = size/2 + 1;

   fftSignal   = (double*)fftw_malloc(fftSignalN * sizeof(double));
   fftSpectrum = (fftw_complex*)fftw_malloc(fftSpectrumN * sizeof(fftw_complex));

   // r2c_1d = real to complex, 1-dimensional
   // signal size = n, spectrum size = (n/2+1) complex numbers
   fftPlan = fftw_plan_dft_r2c_1d(size, fftSignal, fftSpectrum, FFTW_ESTIMATE);
   fftPlan_inverse = fftw_plan_dft_c2r_1d(size, fftSpectrum, fftSignal,
		     FFTW_ESTIMATE);
                     // Also can use: FFTW_ESTIMATE | FFTW_PRESERVE_INPUT
		     // if you don't want the input to be destroyed
		     // during the transform process.

   // http://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags
   //    FFTW_ESTIMATE
   //    FFTW_MEASURE
   //    FFTW_PATIENT
   //    FFTW_EXHAUSTIVE

   if (fftPlan == NULL || fftPlan_inverse == NULL) {
      deinitialize();
      return 0; // unsuccessful initialization
   }

   return 1; // successful initialization
}




//////////////////////////////
//
// MazurkaTransformer::deinitialize -- Clean up transform data.
//

void MazurkaTransformer::deinitialize(void) {
   if (fftPlan != NULL) {
      fftw_destroy_plan(fftPlan);
      fftPlan = NULL;
   }
   if (fftPlan_inverse != NULL) {
      fftw_destroy_plan(fftPlan_inverse);
      fftPlan = NULL;
   }
   if (fftSignal != NULL) {
      fftw_free(fftSignal);
      fftSignal = NULL;
   }
   if (fftSpectrum != NULL) {
      fftw_free(fftSpectrum);
      fftSpectrum = NULL;
   }
   fftSignalN   = 0;
   fftSignalNd2 = 0;
   fftSpectrumN = 0;
}




