/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "filter.h"

// Function to compute the impulse response "h" based on the sinc function
void generateLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h)
{
    h.clear();
    h.resize(num_taps, 0.0);
    
    real NormCutOff = Fc / (Fs / 2);
    real Mid_Index = (num_taps - 1) / 2;
    
    for (int i = 0; i < num_taps; i++) {
        if (i == Mid_Index) {
            h[i] = NormCutOff;
        } else {
            h[i] = NormCutOff * (sin(PI * NormCutOff * (i - Mid_Index))) / (PI * NormCutOff * (i - Mid_Index));
        }
        h[i] *= (1.0 / 2 - (1.0 / 2) * cos((2 * i * PI) / (num_taps - 1)));
    }
}

void generateHPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h) {
    h.clear();
    h.resize(num_taps, 0.0);

    real NormCutOff = Fc / (Fs / 2);
    real Mid_Index = (num_taps - 1) / 2;

    for (int i = 0; i < num_taps; i++) {
        if (i == Mid_Index) {
            h[i] = 1.0 - NormCutOff;
        } else {
            h[i] = -NormCutOff * (sin(PI * NormCutOff * (i - Mid_Index)) / (PI * NormCutOff * (i - Mid_Index)));
        }
        h[i] *= (1.0 / 2 - (1.0 / 2) * cos((2 * i * PI) / (num_taps - 1)));
    }
    
    // Apply spectral inversion to convert LPF to HPF
    for (int i = 0; i < num_taps; i++) {
        h[i] = (i == Mid_Index) ? (1 - h[i]) : (-h[i]);
    }
}


void generateBPF(std::vector<dy4::real> &h, dy4::real f_low, dy4::real f_high, dy4::real fs, int taps) {
    h.clear();
    h.resize(taps, 0.0);
    
    int centre = (taps - 1) / 2;
    dy4::real f_mid = (f_low + f_high) / 2;
    dy4::real scale_factor = 0.0;
    
    // Compute the raw bandpass impulse response
    for (int k = 0; k < taps; k++) {
        int n = k - centre;
        if (n == 0) {
            h[k] = 2.0 / fs * (f_high - f_low);
        } else {
            h[k] = (sin(2 * M_PI * f_high * n / fs) - sin(2 * M_PI * f_low * n / fs)) / (M_PI * n);
        }
        // Apply Hann window:
        h[k] *= 0.5 * (1 - cos(2 * M_PI * k / (taps - 1)));
        
        // Accumulate scale factor using the desired frequency response at f_mid:
        scale_factor += h[k] * cos(2 * M_PI * n * f_mid / fs);
    }
    
    // Normalize the coefficients so that the frequency response at f_mid is 1.
    for (int k = 0; k < taps; k++) {
        h[k] /= scale_factor;
    }
}

// Function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h)
{
	// Bring your own functionality
    // Bring your own functionality
	// Allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size() + h.size() - 1, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

	for(unsigned int c=0; c<y.size(); c++) {
		for(int l = h.size()-1; l >= 0; l--) {
			if (c-l < 0) {
				y[c] += 0;
			} else {
				y[c] += h[l] * x[c-l];
			}
		}
	}
}

void convolveBlock(const std::vector<dy4::real>& signal, const std::vector<dy4::real>& coeffs, std::vector<dy4::real>& state, std::vector<dy4::real>& output) {
    for (size_t i = 0; i < signal.size(); i++) {
        for (size_t j = 0; j < coeffs.size(); j++) {
            if (i >= j) {
                output[i] += signal[i - j] * coeffs[j];
            } else {
                output[i] += state[state.size() - j + i] * coeffs[j];
            }
        }
    }
    state.assign(signal.end() - coeffs.size() + 1, signal.end());
}

void convolveAndDecimate(const std::vector<dy4::real>& x, const std::vector<dy4::real>& h, std::vector<dy4::real>& state, int d, std::vector<dy4::real>& y) {
    int y_n = 0;
    
    for (int n = 0; n < x.size(); n += d) {
        dy4::real sum = 0.0;

        for (int k = 0; k < h.size(); k++) {
            if (n >= k) {
                sum += h[k] * x[n - k];
            } else {
                int state_idx = state.size() + n - k;

                if (state_idx >= 0 && state_idx < state.size()) {
                    sum += h[k] * state[state_idx];
                }
            }
        }

        y[y_n++] = sum;
    }

    int state_len = h.size() - 1;
    if (x.size() >= state_len) {
        for (int i = 0; i < state_len; i++) {
            state[i] = x[x.size() - state_len + i];
        }
    }
}

void ConvolveFast(const std::vector<dy4::real>& x, 
                  const std::vector<dy4::real>& h, 
                  std::vector<dy4::real>& state, 
                  int u, int d, 
                  std::vector<dy4::real>& y) {
    int out_size = x.size() * u / d;
    y.resize(out_size, 0.0);
    
    for (int n = 0; n < out_size; n++) {
        dy4::real sum = 0.0;
        // Effective index (for u = d = 1, effective = n)
        int effective = n * d / u;  
        for (int j = 0; j < h.size(); j++) {
            int idx = effective - j;
            if (idx >= 0) {
                sum += h[j] * x[idx];
            } else {
                sum += h[j] * state[state.size() + idx];
            }
        }
        y[n] = sum * u;
    }
    
    int state_len = h.size() - 1;
    if (x.size() >= (unsigned)state_len) {
        for (int i = 0; i < state_len; i++) {
            state[i] = x[x.size() - state_len + i];
        }
    }
}






// The rest of the code in this function is to be completed by you
// based on your understanding and the Python

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void convolveFIR_inefficient(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);
    
    for (auto n = 0; n < (int)y.size(); n++) {
        for (auto k = 0; k < (int)x.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)h.size())
                y[n] += x[k] * h[n - k];
        }
    }
}

void convolveFIR_reference(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    for (auto n = 0; n < (int)y.size(); n++) {
        for (auto k = 0; k < (int)h.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)x.size())
                y[n] += h[k] * x[n - k];
        }
    }
}

void convolveFIR_loop_x(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    int N = x.size();
    int M = h.size();

    for(int n=0; n < N + M - 1; n++) {
        for(int k=0; k<N; k++) {
            if((n-k >= 0) && (n-k < M)) {
                y[n] += x[k] * h[n-k];
            }
        }
    }
}

void convolveFIR_loop_h(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    int N = x.size();
    int M = h.size();

    for(int n=0; n < N + M - 1; n++) {
        for(int k=0; k<M; k++) {
            if((n-k >= 0) && (n-k < N)) {
                y[n] += h[k] * x[n-k];
            }
        }
    }
}

void convolveFIR_Valid_Range(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    int N = x.size();
    int M = h.size();

    for (auto n = 0; n < N + M - 1; n++) {
        int kmin = std::max(0, n - N + 1);
        int kmax = std::min(M, n + 1);
        
        for (auto k = kmin; k < kmax; k++) {
            y[n] += h[k] * x[n - k];
        }
    }
}

void convolveFIR_Range_and_Unrolling(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    int N = x.size();
    int M = h.size();

    for (auto n = 0; n < N + M - 1; n++) {
        int kmin = std::max(0, n - N + 1);
        int kmax = std::min(M, n+1);
        int unrolled_limit = kmax - ((kmax-kmin)%4);
        
        for (int k = kmin; k < unrolled_limit; k += 4) {
            y[n] += (h[k] * x[n - k]) + (h[k+1]*x[n-(k+1)]) + (h[k+2]*x[n-(k+2)]) + (h[k+3]*x[n-(k+3)]);
        }
        
        for (int k = unrolled_limit; k < kmax; k++){
            y[n] += h[k]*x[n-k];
        }
    }
}

void convolveFIR_il4a(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    int N = x.size();
    int M = h.size();

    for(int k=0; k<M; k++) {
        for(int n=0; n < N + M - 1; n++) {
            if((n-k >= 0) && (n-k < N)) {
                y[n] += h[k] * x[n-k];
            }
        }
    }
}

void convolveFIR_il4b(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    int N = x.size();
    int M = h.size();

    for(int k=0; k < M; k++) {
        int nmin = k;
        int nmax = N + k;

        for(int n = nmin; n < nmax; n++) {
            y[n] += h[k] * x[n - k];
        }
    }
}

void convolveFIR_il4c(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    int N = x.size();
    int M = h.size();

    for (int k=0; k < M; k++) {
        int nmin = k;
        int nmax = N + k;
        int unrolled_limit = nmax - ((nmax-nmin)%4);
        
        for (int n = nmin; n < unrolled_limit; n += 4) {
            y[n+0] += (h[k] * x[(n+0) - k]);
            y[n+1] += (h[k] * x[(n+1) - k]);
            y[n+2] += (h[k] * x[(n+2) - k]);
            y[n+3] += (h[k] * x[(n+3) - k]);
        }
        
        for (int n = unrolled_limit; n < nmax; n++){
            y[n] += h[k]*x[n-k];
        }
    }
}