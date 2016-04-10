//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun Jun 11 05:13:19 PDT 2006
// Last Modified: Sun Jun 11 20:43:58 PDT 2006
// Last Modified: Tue Apr 15 00:38:51 PDT 2008 (added inverse transform)
// Filename:      ...plugin/MazurkaTransformer/MazurkaTransformer.h
// Web:           http://sv.mazurka.org.uk/include/MazurkaTransformer.h
// Syntax:        ANSI99 C++; vamp 1.0 plugin
//
// Description:   Interface to FFTW for doing Fourier Transforms
//

#include "fftw3.h"

#ifndef _MAZURKATRANSFORMER_H_INCLUDED
#define _MAZURKATRANSFORMER_H_INCLUDED

typedef struct {
   double re;
   double im;
} mz_complex;

class MazurkaTransformer {

   public:
                    MazurkaTransformer     (void);
                    MazurkaTransformer     (int size);
                   ~MazurkaTransformer     ();

      int           doTransform            (void);
      int           doInverseTransform     (void);

      int           getSize                (void);
      void          setSize                (int size);
      double&       signalCausal           (int index);
      double&       signalNonCausal        (int index);
      void          zeroSignal             (void);

      mz_complex    getSpectrum            (int index);
      double        getSpectrumMagnitude   (int index);
      double        getSpectrumMagnitudeDb (int index, double reference = 1.0);
      double        getSpectrumSquared     (int index);

      double&       operator[]             (int index);

   protected:

      int           initialize             (int size);
      void          deinitialize           (void);

   private:

      fftw_plan     fftPlan;          // maintenance variable for fftw
      fftw_plan     fftPlan_inverse;  // setup for complex to real ifft
      int           fftSignalN;       // length of signal input
      int           fftSignalNd2;     // 1/2 length of signal input
      int           fftSpectrumN;     // internal size of freq. array
      double       *fftSignal;        // storage space for real signal
      fftw_complex *fftSpectrum;      // storage space for complex spectrum
      // typedef double[2] fftw_complex;

};



#endif // _MAZURKATRANSFORMER_H_INCLUDED

