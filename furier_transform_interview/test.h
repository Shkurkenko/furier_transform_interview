#pragma once

#include <fftw3.h>

namespace test {

	constexpr auto REAL = 0;
	constexpr auto IMAG = 1;
	constexpr auto N = 3;

	void direct_fftw(fftw_complex* in, fftw_complex* out);

	void inverse_fftw(fftw_complex* in, fftw_complex* out);

	void displayComplexFFTW(fftw_complex* y);

	void display_real(fftw_complex* y);

	void log_fftw_test();
}
