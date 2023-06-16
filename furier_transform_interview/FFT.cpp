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
	std::complex<COMPLEX_TYPE> t = y[1];
    y[1] = y[0] - t;
    y[0] += t;
}

void FFT::fft_radix3_base_transform(std::complex<COMPLEX_TYPE>* y)
{
	static const COMPLEX_TYPE k = -sin(2.0 * M_PI / 3.0);
    constexpr std::complex<COMPLEX_TYPE> J = std::complex<COMPLEX_TYPE>(0, 1);
    std::complex<COMPLEX_TYPE> a = y[1] + y[2];
    std::complex<COMPLEX_TYPE> t1 = y[0] - a / 2.0;
    std::complex<COMPLEX_TYPE> t2 = (y[1] - y[2]) * k;
    std::complex<COMPLEX_TYPE> b = J * t2;
    y[0] += a;
    y[1] = t1 + b;
    y[2] = t1 - b;
};

void FFT::fft_radix5_base_transform(std::complex<COMPLEX_TYPE>* y)
{
	static const COMPLEX_TYPE k0 = cos(4.0 * M_PI / 5.0);
    static const COMPLEX_TYPE k10 = -sin(4.0 * M_PI / 5.0);
    static const COMPLEX_TYPE t = sin(2.0 * M_PI / 5.0);
    static const COMPLEX_TYPE k11 = t - k10;
    static const COMPLEX_TYPE k12 = -t - k10;
    static const std::complex<COMPLEX_TYPE> J = std::complex<COMPLEX_TYPE>(0, 1);
    std::complex<COMPLEX_TYPE> a1 = y[1] + y[4];
    std::complex<COMPLEX_TYPE> a2 = y[2] + y[3];
    std::complex<COMPLEX_TYPE> a3 = y[2] - y[3];
    std::complex<COMPLEX_TYPE> a4 = y[1] - y[4];
    std::complex<COMPLEX_TYPE> c5 = k0 * (a1 - a2);
    std::complex<COMPLEX_TYPE> c6 = k10 * (a3 + a4);
    std::complex<COMPLEX_TYPE> d1 = y[0] - a1 / 2.0 - c5;
    std::complex<COMPLEX_TYPE> d2 = y[0] - a2 / 2.0 + c5;
    std::complex<COMPLEX_TYPE> t1 = J * (k11 * a3 + c6);
    std::complex<COMPLEX_TYPE> t2 = J * (k12 * a4 + c6);
    y[0] += a1 + a2;
    y[1] = d1 + t2;
    y[2] = d2 + t1;
    y[3] = d2 - t1;
    y[4] = d1 - t2;
};

void FFT::fft(
	std::vector<std::complex<COMPLEX_TYPE>>& a, DIRECTION direction, std::function<void(std::complex<COMPLEX_TYPE>*)> base_transform
)
{
	int n = (int)a.size();
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
		COMP wlen(cos(ang), sin(ang));

		for (int i = 0; i < n; i += len)
		{
			COMP w(1);

			for (int j = 0; j < len / 2; ++j)
			{
				COMP u = a[i + j], v = a[i + j + len / 2] * w;
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
};

void FFT::chirp_z_transform(COMPLEX_DATA& a, DIRECTION direction)
{
	std::cout << "fft_bluestein";
};

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
};

