#include <iostream>
#include <fftw3.h>
#include <complex>
#include <vector>
#include <cmath>
#include <numbers>

constexpr auto REAL = 0;
constexpr auto IMAG = 1;
constexpr auto N = 3;

constexpr auto PI = 3.14159265358979323846;

typedef std::complex<double> comp;
typedef std::vector<comp> complex_data;

int rev (int num, int lg_n) {
	int res = 0;
	for (int i=0; i<lg_n; ++i)
		if (num & (1<<i))
			res |= 1<<(lg_n-1-i);
	return res;
}

void radix2_fft(complex_data& a, bool inverse) {
	int n = (int) a.size();
	int lg_n = 0;
	while ((1 << lg_n) < n)  ++lg_n;
 
	for (int i=0; i<n; ++i)
		if (i < rev(i,lg_n))
			swap (a[i], a[rev(i,lg_n)]);
 
	for (int len=2; len<=n; len<<=1) {
		double ang = 2*PI/len * (inverse ? -1 : 1);
		comp wlen (cos(ang), sin(ang));
		for (int i=0; i<n; i+=len) {
			comp w (1);
			for (int j=0; j<len/2; ++j) {
				comp u = a[i+j],  v = a[i+j+len/2] * w;
				a[i+j] = u + v;
				a[i+j+len/2] = u - v;
				w *= wlen;
			}
		}
	}
	if (inverse)
		for (int i=0; i<n; ++i)
			a[i] /= n;

}

void blustein_fft(complex_data& a, bool inverse) {

};

void direct_fftw(fftw_complex *in, fftw_complex *out) {
	// Createing plan for fft.
	fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	// Execute the plan
	fftw_execute(p);
	// Do some cleaning
	fftw_destroy_plan(p);
	fftw_cleanup();

}

void inverse_fftw(fftw_complex *in, fftw_complex *out) {
	// Createing plan for fft. 
	fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	// Execute the plan
	fftw_execute(p);
	// Do some cleaning
	fftw_destroy_plan(p);
	fftw_cleanup();
	// Scale the plan output to obtain the exact inverse
	for (int i = 0; i < N; ++i) {
		out[i][REAL] /= N;
		out[i][IMAG] /= N;
	}
}
// Display complex numbers in the forma +/- bi.
void displayComplexFFTW(fftw_complex* y) {
	for (int i = 0; i < N; ++i) {
		if (y[i][IMAG] < 0)
			std::cout << y[i][REAL] << " - " << abs(y[i][IMAG]) << 'i' << std::endl;
		else
			std::cout << y[i][REAL] << " + " << abs(y[i][IMAG]) << 'i' << std::endl;
	}
}
// Display the real parts of complex numbers.
void displayReal(fftw_complex* y) {
	for (int i = 0; i < N; ++i) {
		std::cout << y[i][REAL] << std::endl;
	}
}

// Test
int main() {
	
	// Input and output array
	fftw_complex x[N], y[N];
	
	// Fill array some data
	for (int i = 0; i < N; ++i) {
		x[i][REAL] = i + 1; // i.e., { 1, 2, 3, 4, 5, 6, 7, 8 }
		x[i][IMAG] = 0;
	}

	// Compute direct FFT of in and store the results in out
	direct_fftw(x, y);
	// Display the results
	std::cout << "Direct_FFT = " << std::endl;
	displayComplexFFTW(y);

	// Compute inverse FFT of out and store the results in out
	inverse_fftw(y, x);
	std::cout << "\nInverse_FFT = " << std::endl;
	displayReal(x);

	std::cout << std::endl;


	// My FFT test
	complex_data data(N);

	for (int i = 0; i < N; ++i) {
		data[i] = comp(i + 1.0, 0);
	}
	
	std::cout << "input data: ";
	for (int i = 0; i < data.size(); ++i) {
		std::cout << data[i] << ' ';
	}
	std::cout << std::endl;

	radix2_fft(data, false);
	std::cout << "direct fft: ";
	for (int i = 0; i < data.size(); ++i) {
		std::cout << data[i] << ' ';
	}
	std::cout << std::endl;
	
	std::cout << "inverse fft: ";
	radix2_fft(data, true);
	for (int i = 0; i < data.size(); ++i) {
		std::cout << data[i] << ' ';
	}
		
	return 0;
}