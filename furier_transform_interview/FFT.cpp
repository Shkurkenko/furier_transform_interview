#include <complex>
#include <vector>
#include <iostream>
#include <complex>
#include <functional>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h>

#include "FFT.h"

bool FFT::is_power2(size_t n)
{
	return n && !(n & (n - 1));
}

int FFT::rev(int num, int lg_n)
{
	int res = 0;

	for (int i = 0; i < lg_n; ++i)
	{
		if (num & (1 << i))
		{
			res |= 1 << (lg_n - 1 - i);
		}
	}
	return res;
}

void FFT::fft_radix2_base_transform(std::complex<COMPLEX_TYPE>* y)
{
}

void FFT::fft_radix3_base_transform(std::complex<COMPLEX_TYPE>* y)
{
}

void FFT::fft_radix5_base_transform(std::complex<COMPLEX_TYPE>* y)
{
}

void FFT::fft(
	std::vector<std::complex<COMPLEX_TYPE>>& a, DIRECTION direction, std::function<void(std::complex<COMPLEX_TYPE>*)> base_transform
)
{
	int n = (int) a.size();
	int lg_n = 0;

	while ((1 << lg_n) < n)  ++lg_n;
 
	for (int i = 0; i < n; ++i)
	{
		if (i < rev(i, lg_n)) 
		{
			swap(a[i], a[rev(i, lg_n)]);
		}
	}
 
	for (int len = 2; len <= n; len <<= 1) 
	{
		COMPLEX_TYPE ang = 2 * M_PI / len * static_cast<int>(direction);
		COMP wlen (cos(ang), sin(ang));

		for (int i = 0; i < n; i += len) 
		{
			COMP w (1);

			for (int j = 0; j < len / 2; ++j) 
			{
				COMP u = a[i + j],  v = a[i + j + len / 2] * w;
				a[i + j] = u + v;
				a[i + j + len / 2] = u - v;
				w *= wlen;
			}
		}
	}

	if (direction == DIRECTION::INVERSE) 
	{
		for (int i = 0; i < n; ++i)
		{
			a[i] /= n;
		}
	}
}

void FFT::chirp_z_transform(COMPLEX_DATA& a, DIRECTION direction)
{
	std::cout << "fft_bluestein";
}

void FFT::transform(COMPLEX_DATA& a, DIRECTION direction)
{
	size_t n = a.size();

	if (n < 2)
	{
		return;
	}

	if (is_power2(n))
	{
		return fft(a, direction, fft_radix2_base_transform);
	}
	else if (a.size() % 3 and is_power2(n / 3))
	{
		return fft(a, direction, fft_radix3_base_transform);
	}
	else if (a.size() % 5 and is_power2(n / 5))
	{
		return fft(a, direction, fft_radix5_base_transform);
	}
	else {
		return chirp_z_transform(a, direction);
	}

}

