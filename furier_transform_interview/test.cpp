#include <iostream>
#include <vector>
#include <complex>

#include "test.h"

void test::direct_fftw(fftw_complex *in, fftw_complex *out) {
	// Createing plan for fft.
	fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	// Execute the plan
	fftw_execute(p);
	// Do some cleaning
	fftw_destroy_plan(p);
	fftw_cleanup();

}

void test::inverse_fftw(fftw_complex *in, fftw_complex *out) {
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

void test::displayComplexFFTW(fftw_complex* y) {
	for (int i = 0; i < N; ++i) {
		if (y[i][IMAG] < 0)
			std::cout << y[i][REAL] << " - " << abs(y[i][IMAG]) << 'i' << std::endl;
		else
			std::cout << y[i][REAL] << " + " << abs(y[i][IMAG]) << 'i' << std::endl;
	}
}

void test::display_real(fftw_complex* y) {
	for (int i = 0; i < N; ++i) {
		std::cout << y[i][REAL] << std::endl;
	}
}

void test::log_fftw_test() {
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
	display_real(x);
}
