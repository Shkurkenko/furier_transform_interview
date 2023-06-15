#pragma once

#ifndef COMPLEX_TYPE
#define COMPLEX_TYPE double
#endif

typedef std::complex<COMPLEX_TYPE> COMP;
typedef std::vector<COMP> COMPLEX_DATA;

enum class DIRECTION {
	DIRECT = 1,
	INVERSE = -1,
};

enum class COMPLEX_PARTS {
	REAL = 0,
	IMAG = 1,
};

class FFT
{
private:
	bool is_power2(size_t n);

	void fft_radix2_base_transform(std::complex<COMPLEX_TYPE>* y);

	void fft_radix3_base_transform(std::complex<COMPLEX_TYPE>* y);

	void fft_radix5_base_transform(std::complex<COMPLEX_TYPE>* y);

	void fft(
		std::vector<std::complex<COMPLEX_TYPE>>& a, DIRECTION direction, std::function<void(std::complex<COMPLEX_TYPE>*)> base_transform
	);

	void chirp_z_transform(COMPLEX_DATA& a, DIRECTION direction);

	int rev(int num, int lg_n);

public:
	void transform(COMPLEX_DATA& a, DIRECTION direction);

};

