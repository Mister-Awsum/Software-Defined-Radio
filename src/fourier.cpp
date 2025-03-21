/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// Source code for Fourier-family of functions

#include "fourier.h"

// Just DFT function (no FFT)
void DFT(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * PI * (k * m) / x.size());
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// Function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<real>> &Xf, std::vector<real> &Xmag)
{void DFT_1d_twiddle(const std::vector<real> &x, std::vector<std::complex<real>> &Xf);
	Xmag.clear();
	Xmag.resize(Xf.size(), real(0));
	for (int i = 0; i < (int)Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i]) / Xf.size();
	}
}

// Add your own code to estimate the PSD
void estimatePSD(const std::vector<real> &samples, unsigned int Fs, std::vector<real> &freq, std::vector<real> &psd_est) {
    int freq_bins = NFFT;
    real df = 1.0 * Fs / freq_bins;

    // Create frequency vector
    freq.clear();
    for (double f = 0; f < Fs / 2; f += df) {
        freq.push_back(f);
    }

    // Create and apply Hann window
    std::vector<real> hann(freq_bins, 0.0);
    for (unsigned int i = 0; i < hann.size(); i++) {
        hann[i] = 0.5 * (1 - cos(2 * PI * i / (freq_bins - 1)));
    }

    // Compute number of segments
    int no_segments = int(floor(samples.size() / float(freq_bins)));

    // PSD storage
    std::vector<real> psd_list;
    std::vector<real> psd_seg(freq_bins / 2, 0.0);

    for (int k = 0; k < no_segments; k++) {
        // Apply Hann window
        std::vector<real> windowed_samples(freq_bins, 0.0);
        for (int i = 0; i < freq_bins; i++) {
            windowed_samples[i] = samples[k * freq_bins + i] * hann[i];
        }

        // Compute Fourier Transform
        std::vector<std::complex<real>> Xf;
        DFT(windowed_samples, Xf);

        // Keep only positive half of the spectrum
        std::vector<std::complex<real>> Xf_half(Xf.begin(), Xf.begin() + freq_bins / 2);

        // Compute power spectrum
        for (unsigned int i = 0; i < Xf_half.size(); i++) {
            real power = 2 * (1.0 / (Fs * freq_bins / 2)) * std::norm(Xf_half[i]);
            psd_list.push_back(power);
        }
    }

    // Average across segments
    for (int k = 0; k < int(freq_bins / 2); k++) {
        for (int l = 0; l < no_segments; l++) {
            psd_seg[k] += psd_list[k + l * int(freq_bins / 2)];
        }
        psd_seg[k] /= no_segments;
    }

    // Convert to dB
    psd_est.clear();
    psd_est.resize(int(freq_bins / 2));
    for (int k = 0; k < int(freq_bins / 2); k++) {
        psd_est[k] = 10 * log10(std::max(psd_seg[k], 1e-10)); // Avoid log10(0)
    }
}


//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void DFT_reference(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / x.size());
			Xf[m] +=  + x[k] * std::exp(expval);
		}
	}
}

void DFT_init_bins(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / N);
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

void DFT_code_motion(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	std::complex<real> base_exp(0, -2 * M_PI / N);
	
	for (int m = 0; m < N; m++) {
	        std::complex<real> inter_exp = base_exp * std::complex<real> (m, 0.);
		for (int k = 0; k < N; k++) {
        	        std::complex<real> twiddle_exp = inter_exp * std::complex<real> (k, 0.);
			Xf[m] += x[k] * std::exp(twiddle_exp);
		}
	}
}

void generate_DFT_twiddles(const int& N, std::vector<std::complex<real>> &Twiddle1D) {

	Twiddle1D.resize(N);
	for (int k = 0; k < N; k++) {
		std::complex<real> expval(0, -2 * M_PI * k / N);
		Twiddle1D[k] = std::exp(expval);
	}
}

void generate_DFT_matrix(const int& N, std::vector<std::vector<std::complex<real>>> &Twiddle2D) {

	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));
	std::vector<std::complex<real>> Twiddle1D;
	generate_DFT_twiddles(N, Twiddle1D);

	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			Twiddle2D[m][k] = Twiddle1D[(k * m) % N];
		}
	}
}

void DFT_1d_twiddle(const std::vector<real> &x, std::vector<std::complex<real>> &Twiddle1D, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();
	
	// Main Loop
	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
        	Xf[m] += x[k] * Twiddle1D[(k*m)%N];
		}
	}
}

void DFT_2d_twiddle(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();	
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	
	// Main Loop
	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
        	Xf[m] += x[k] * Twiddle2D[m][k];
		}
	}
}

void DFT_loop_m(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf) {
	int N = (int)x.size();
	
	for (int m=0; m<N; m++) {
		for (int k=0; k<N; k++) {
			Xf[m] += x[k] * Twiddle2D[m][k];
		}
	}
}

void DFT_loop_k(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf) {
	int N = (int)x.size();
	
	for (int k=0; k<N; k++) {
		for (int m=0; m<N; m++) {
			Xf[m] += x[k] * Twiddle2D[m][k];
		}
	}
}

void DFT_unrolling(const std::vector<real> &x, std::vector<std::vector<std::complex<real>>> &Twiddle2D, std::vector<std::complex<real>> &Xf) {
	int N = (int)x.size();
	int unrolled_limit = N - (N%4);

	for(int m=0; m<N; m++) {
		for(int k=0; k<unrolled_limit; k += 4) {
			Xf[m] += x[k] * Twiddle2D[m][k] + x[k+1] * Twiddle2D[m][k+1] + x[k+2] * Twiddle2D[m][k+2] + x[k+3] * Twiddle2D[m][k+3];
		}

		for(int k=unrolled_limit; k<N; k++) {
			Xf[m] += x[k] * Twiddle2D[m][k];
		}
	}
}

// Calculate PSD using matrix operations
void psd_matrix(const std::vector<double>& samples, const real Fs, std::vector<double>& freq, std::vector<double>& psdEst) {
	int freq_bins = NFFT;
	real df = 1.0 * Fs / freq_bins;

	freq.clear();

	for (double f = 0; f < Fs / 2; f += df) {
        freq.push_back(f);
    }

    int noSegments = std::floor(samples.size() / freq_bins);

	// generate the DFT matrix for the given size N
    std::vector<std::vector<std::complex<double>>> dftMatrix(freq_bins, std::vector<std::complex<double>>(freq_bins));

    for (int m = 0; m < freq_bins; m++) {
        for (int k = 0; k < freq_bins; k++) {
            double angle = 2 * PI * (-k) * m / freq_bins;
            dftMatrix[m][k] = std::exp(std::complex<double>(0, angle));
        }
    }

    // generate the Hann window for the given size N
	std::vector<real> hannWindow;
	hannWindow.resize(freq_bins);

    for (int i = 0; i < freq_bins; i++) {
        hannWindow[i] = 0.5 * (1 - std::cos(2 * PI * i / (freq_bins - 1)));
    }

    // apply Hann window and perform matrix multiplication using nested loops
	std::vector<std::vector<std::complex<double>>> Xf(noSegments, std::vector<std::complex<double>>(freq_bins));
    for (int seg = 0; seg < noSegments; seg++) {
        for (int m = 0; m < freq_bins; m++) {
            std::complex<double> sum = 0;

            for (int k = 0; k < freq_bins; k++) {
                sum += samples[seg * freq_bins + k] * hannWindow[k] * dftMatrix[m][k];
            }

            Xf[seg][m] = sum;
        }
    }

    // compute power, keep only positive frequencies, average across segments, and convert to dB
	psdEst.resize(freq_bins / 2);

    for (int m = 0; m < freq_bins / 2; m++) {
        double sumPower = 0.0;
        
		for (int seg = 0; seg < noSegments; seg++) {
            sumPower += (2.0 / (Fs * freq_bins / 2)) * std::norm(Xf[seg][m]);
        }

        psdEst[m] = 10 * std::log10(sumPower / noSegments);
    }
}

// Custom Optimizations
void psd_2d_twiddle(const std::vector<double>& samples, const real Fs, std::vector<double>& freq, std::vector<double>& psdEst) {
	int freq_bins = NFFT;
	real df = 1.0 * Fs / freq_bins;

	// Create frequency vector
	freq.clear();
	for (double f = 0; f < Fs / 2; f += df) {
        freq.push_back(f);
    }

	int noSegments = std::floor(samples.size() / freq_bins);

	// Generate the precomputed 2D DFT matrix (only once)
	std::vector<std::vector<std::complex<real>>> Twiddle2D;
	generate_DFT_matrix(freq_bins, Twiddle2D);

	// Generate the Hann window
	std::vector<real> hannWindow(freq_bins);
	for (int i = 0; i < freq_bins; i++) {
		hannWindow[i] = 0.5 * (1 - std::cos(2 * PI * i / (freq_bins - 1)));
	}

	// Allocate space for DFT output
	std::vector<std::vector<std::complex<double>>> Xf_all(noSegments, std::vector<std::complex<double>>(freq_bins));

	for (int seg = 0; seg < noSegments; seg++) {
		// Apply Hann window
		std::vector<real> windowed_samples(freq_bins);
		for (int i = 0; i < freq_bins; i++) {
			windowed_samples[i] = samples[seg * freq_bins + i] * hannWindow[i];
		}

		// Compute DFT using precomputed twiddle matrix (Inlined DFT_2d_twiddle)
		std::vector<std::complex<real>> Xf(freq_bins, std::complex<real>(0.0, 0.0));
		for (int m = 0; m < freq_bins; m++) {
			for (int k = 0; k < freq_bins; k++) {
				Xf[m] += windowed_samples[k] * Twiddle2D[m][k];
			}
		}
		Xf_all[seg] = Xf; // Store the transformed segment
	}

	// Compute power, keep only positive frequencies, average across segments, and convert to dB
	psdEst.resize(freq_bins / 2);

	for (int m = 0; m < freq_bins / 2; m++) {
		double sumPower = 0.0;

		for (int seg = 0; seg < noSegments; seg++) {
			sumPower += (2.0 / (Fs * freq_bins / 2)) * std::norm(Xf_all[seg][m]);
		}

		psdEst[m] = 10 * std::log10(sumPower / noSegments);
	}
}

void psd_unrolling(const std::vector<double>& samples, const real Fs, std::vector<double>& freq, std::vector<double>& psdEst) {
	int freq_bins = NFFT;
	real df = 1.0 * Fs / freq_bins;

	freq.clear();

	for (double f = 0; f < Fs / 2; f += df) {
        freq.push_back(f);
    }

    int noSegments = std::floor(samples.size() / freq_bins);

	// generate the DFT matrix for the given size N
    std::vector<std::vector<std::complex<double>>> dftMatrix(freq_bins, std::vector<std::complex<double>>(freq_bins));

    for (int m = 0; m < freq_bins; m++) {
        for (int k = 0; k < freq_bins; k++) {
            double angle = 2 * PI * (-k) * m / freq_bins;
            dftMatrix[m][k] = std::exp(std::complex<double>(0, angle));
        }
    }

    // generate the Hann window for the given size N
	std::vector<real> hannWindow;
	hannWindow.resize(freq_bins);

    for (int i = 0; i < freq_bins; i++) {
        hannWindow[i] = 0.5 * (1 - std::cos(2 * PI * i / (freq_bins - 1)));
    }

    // apply Hann window and perform matrix multiplication using nested loops
	std::vector<std::vector<std::complex<double>>> Xf(noSegments, std::vector<std::complex<double>>(freq_bins));
    for (int seg = 0; seg < noSegments; seg++) {
		for (int m = 0; m < freq_bins; m++) {
			std::complex<double> sum = 0;
			int k = 0;

			// Loop unrolling for better performance
			for (; k + 3 < freq_bins; k += 4) {
				sum += samples[seg * freq_bins + k] * hannWindow[k] * dftMatrix[m][k] +
				       samples[seg * freq_bins + k + 1] * hannWindow[k + 1] * dftMatrix[m][k + 1] +
				       samples[seg * freq_bins + k + 2] * hannWindow[k + 2] * dftMatrix[m][k + 2] +
				       samples[seg * freq_bins + k + 3] * hannWindow[k + 3] * dftMatrix[m][k + 3];
			}

			// Process remaining elements
			for (; k < freq_bins; k++) {
				sum += samples[seg * freq_bins + k] * hannWindow[k] * dftMatrix[m][k];
			}

			Xf[seg][m] = sum;
		}
	}

    // compute power, keep only positive frequencies, average across segments, and convert to dB
	psdEst.resize(freq_bins / 2);

    for (int m = 0; m < freq_bins / 2; m++) {
        double sumPower = 0.0;
        
		for (int seg = 0; seg < noSegments; seg++) {
            sumPower += (2.0 / (Fs * freq_bins / 2)) * std::norm(Xf[seg][m]);
        }

		psdEst[m] = 10 * std::log10(sumPower / noSegments);
    }
}