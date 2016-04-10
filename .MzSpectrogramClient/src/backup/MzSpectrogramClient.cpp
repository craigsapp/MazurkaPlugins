//
// Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
// Creation Date: Sun May  7 18:53:16 PDT 2006
// Last Modified: Mon May  8 07:28:30 PDT 2006
// Filename:      MzSpectrogramClient.cpp
// Syntax:        C++; vamp plugin
//


#include "MzSpectrogramClient.h" 
#include "math.h"


///////////////////////////////
//
// MzSpectrogramClient::MzSpectrum -- class constructor.
//

MzSpectrogramClient::MzSpectrogramClient(float inputSampleRate) :
      Plugin(inputSampleRate) {
   mz_rtimebuffer = NULL;
   mz_rfreqbuffer = NULL;
   mz_ifreqbuffer = NULL;
   mz_window      = NULL;
   reset();
}



///////////////////////////////
//
// MzSpectrogramClient::~MzSpectrum -- class destructor.
//

MzSpectrogramClient::~MzSpectrogramClient() {
   if (mz_rtimebuffer != NULL) {
      delete [] mz_rtimebuffer;
      mz_rtimebuffer = NULL;
   }
   if (mz_rfreqbuffer != NULL) {
      delete [] mz_rfreqbuffer;
      mz_rfreqbuffer = NULL;
   }
   if (mz_ifreqbuffer != NULL) {
      delete [] mz_ifreqbuffer;
      mz_ifreqbuffer = NULL;
   }
   if (mz_window != NULL) {
      delete [] mz_window;
      mz_window = NULL;
   }
}



////////////////////////////////////////////////////////////
//
// polymorphic functions inherited from Plugin:
//

void MzSpectrogramClient::reset(void) {
   if (mz_rtimebuffer != NULL) {
      delete [] mz_rtimebuffer;
      mz_rtimebuffer = NULL;
   }
   if (mz_rfreqbuffer != NULL) {
      delete [] mz_rfreqbuffer;
      mz_rfreqbuffer = NULL;
   }
   if (mz_ifreqbuffer != NULL) {
      delete [] mz_ifreqbuffer;
      mz_ifreqbuffer = NULL;
   }
   if (mz_window != NULL) {
      delete [] mz_window;
      mz_window = NULL;
   }
   initialise(0, 0, 0);
}

bool MzSpectrogramClient::initialise(size_t channels, size_t stepSize, 
      size_t blockSize) {
   if (channels < getMinChannelCount() ||
       channels > getMaxChannelCount()) {
      return false;
   }

   if (blockSize <= 0) {
      return false;
   }

   mz_channels  = channels;
   mz_stepsize  = stepSize;
   mz_blocksize = blockSize;

   if (mz_rtimebuffer != NULL) {
      delete [] mz_rtimebuffer;
      mz_rtimebuffer = NULL;
   }
   if (mz_rfreqbuffer != NULL) {
      delete [] mz_rfreqbuffer;
      mz_rfreqbuffer = NULL;
   }
   if (mz_ifreqbuffer != NULL) {
      delete [] mz_ifreqbuffer;
      mz_ifreqbuffer = NULL;
   }
   if (mz_window != NULL) {
      delete [] mz_window;
      mz_window = NULL;
   }

   mz_rtimebuffer = new double[mz_blocksize];
   mz_rfreqbuffer = new double[mz_blocksize];
   mz_ifreqbuffer = new double[mz_blocksize];
   mz_window      = new double[mz_blocksize];

   makeHannWindow(mz_window, mz_blocksize);

   return true;
}



////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from PluginBase:
//

std::string MzSpectrogramClient::getName(void) const
   { return "mzspectrumclient"; }

std::string MzSpectrogramClient::getDescription(void) const
   { return "Spectrum Client"; }

std::string MzSpectrogramClient::getMaker(void) const
   { return "Craig Stuart Sapp <craig@ccrma.stanford.edu>"; }

std::string MzSpectrogramClient::getCopyright(void) const
   { return "2006 Craig Stuart Sapp"; }

int MzSpectrogramClient::getPluginVersion(void) const
   { return 200605080; }



////////////////////////////////////////////////////////////
//
// required polymorphic functions inherited from Plugin:
//

MzSpectrogramClient::InputDomain 
MzSpectrogramClient::getInputDomain(void) const 
   { return TimeDomain; }


MzSpectrogramClient::OutputList 
MzSpectrogramClient::getOutputDescriptors(void) const {
   OutputList list;

   OutputDescriptor od;
   od.name             = "magintude";
   od.description      = "Magnitude Spectrum";
   od.unit             = "";
   od.hasFixedBinCount = true;
   od.binCount         = getBinCount();
   od.hasKnownExtents  = false;
   od.isQuantized      = false;
   od.sampleType       = OutputDescriptor::OneSamplePerStep;
   list.push_back(od);

   return list; 
}


MzSpectrogramClient::FeatureSet 
MzSpectrogramClient::getRemainingFeatures(void) {
   // no remaining features
   return FeatureSet();
}


MzSpectrogramClient::FeatureSet MzSpectrogramClient::process(float **input, 
      Vamp::RealTime timestamp) {

   if (mz_stepsize <= 0 || mz_rtimebuffer == NULL) {
      std::cerr << "ERROR: MzSpectrogramClient::process: "
                << "MzSpectrogramClient has not been initialized"
                << std::endl;
      return FeatureSet();
   }


   windowSignal(mz_rtimebuffer, mz_window, input[0], mz_blocksize);
   // fft(mz_blocksize, false, mz_rtimebuffer, NULL, 
   //       mz_rfreqbuffer, mz_ifreqbuffer);
   
   // makeMagnitudeSpectrum(mz_rtimebuffer, mz_rfreqbuffer, mz_ifreqbuffer, 
   //       mz_blocksize);

   FeatureSet returnFeatures;
   Feature feature;

   feature.hasTimestamp = false;
   int i;
   for (i=0; i<getBinCount(); i++) {
      feature.values.push_back(mz_rtimebuffer[i]);
   }
   returnFeatures[0].push_back(feature);

   return returnFeatures;
}



////////////////////////////////////////////////////////////
//
// non-interface functions 
//


//////////////////////////////
//
// MzSpectrogramClient::getBinCount -- 
//

int MzSpectrogramClient::getBinCount(void) const {
   return mz_blocksize / 4;
}



//////////////////////////////
//
// MzSpectrogramClient::makeHannWindow -- create a raised cosine (Han)
//     window.
//

void MzSpectrogramClient::makeHannWindow(double* output, int blocksize) {
   int i;
   for (i=0; i<blocksize; i++) {
      output[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i/(blocksize-1)));
      // make a square window for now for debugging:
      output[i] = 1.0;
   }
}



//////////////////////////////
//
// MzSpectrogramClient::windowSignal -- window the time-domain signal.
//

void MzSpectrogramClient::windowSignal(double* output, double* window, 
      float* input, int blocksize) {
   int i;
   for (i=0; i<blocksize; i++) {
      //output[i] = window[i] * input[i];
      // for debugging:
      output[i] = input[i];
   }
}



//////////////////////////////
//
// makeMagnitudeSpectrum -- 
//

#define MAGNITUDE(x, y) sqrt((x)*(x)+(y)*(y))
#define ZEROLOG 120.0

void MzSpectrogramClient::makeMagnitudeSpectrum(double* output, 
      double* realfreq, double* imagfreq, size_t blocksize) {
   int i;
   int count = getBinCount();
   for (i=0; i<count; i++) {
      mz_rtimebuffer[i] = MAGNITUDE(realfreq[i], imagfreq[i]);

      // convert the amplitude magnitude into decibels
      if (mz_rtimebuffer[i] <= 0.0) {
         mz_rtimebuffer[i] = ZEROLOG;
      } else {
         mz_rtimebuffer[i] = 20.0 * log10(mz_rtimebuffer[i]);
      }
   }
}



//////////////////////////////
//
//  MzSpectrogramClient::fft -- calculate the Fast Fourier Transform.
//     Code stolen from the vamp plugin sdk.
//

void MzSpectrogramClient::fft(int n, bool inverse, double *ri, 
      double *ii, double *ro, double *io) {

    if (!ri || !ro || !io) return;

    int bits;
    int i, j, k, m;
    int blockSize, blockEnd;
    double tr, ti;

    if (n < 2) return;
    if (n & (n-1)) return;
    double angle = 2.0 * M_PI;
    if (inverse) angle = -angle;

    for (i = 0; ; ++i) {
	if (n & (1 << i)) {
	    bits = i;
	    break;
	}
    }

    static int tableSize = 0;
    static int *table = 0;

    if (tableSize != n) {
	delete[] table;
	table = new int[n];
	for (i = 0; i < n; ++i) {
	    m = i;
	    for (j = k = 0; j < bits; ++j) {
		k = (k << 1) | (m & 1);
		m >>= 1;
	    }
	    table[i] = k;
	}
	tableSize = n;
    }

    if (ii) {
	for (i = 0; i < n; ++i) {
	    ro[table[i]] = ri[i];
	    io[table[i]] = ii[i];
	}
    } else {
	for (i = 0; i < n; ++i) {
	    ro[table[i]] = ri[i];
	    io[table[i]] = 0.0;
	}
    }

    blockEnd = 1;

    for (blockSize = 2; blockSize <= n; blockSize <<= 1) {
	double delta = angle / (double)blockSize;
	double sm2 = -sin(-2 * delta);
	double sm1 = -sin(-delta);
	double cm2 = cos(-2 * delta);
	double cm1 = cos(-delta);
	double w = 2 * cm1;
	double ar[3], ai[3];

	for (i = 0; i < n; i += blockSize) {
	    ar[2] = cm2;
	    ar[1] = cm1;
	    ai[2] = sm2;
	    ai[1] = sm1;
	    for (j = i, m = 0; m < blockEnd; j++, m++) {
		ar[0] = w * ar[1] - ar[2];
		ar[2] = ar[1];
		ar[1] = ar[0];
		ai[0] = w * ai[1] - ai[2];
		ai[2] = ai[1];
		ai[1] = ai[0];
		k = j + blockEnd;
		tr = ar[0] * ro[k] - ai[0] * io[k];
		ti = ar[0] * io[k] + ai[0] * ro[k];
		ro[k] = ro[j] - tr;
		io[k] = io[j] - ti;
		ro[j] += tr;
		io[j] += ti;
	    }
	}

	blockEnd = blockSize;
    }

    if (inverse) {

	double denom = (double)n;

	for (i = 0; i < n; i++) {
	    ro[i] /= denom;
	    io[i] /= denom;
	}
    }
}



