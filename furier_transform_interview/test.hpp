#pragma once

# ifndef DPF_COMPLEX_TYPE
#define COMPLEX_TYPE double
#else
define COMPLEX_TYPE DPF_COMPLEX_TYPE
#endif

#include <complex>
#include <vector>
#include <numbers>
#include <algorithm>
#include <functional>
#include <thread>
#include <map>

void fft_radix_2_base_transform(std::complex<COMPLEX_TYPE>* y)
{
    std::complex<COMPLEX_TYPE> t = y[1];
    y[1] = y[0] - t;
    y[0] += t;
};

void fft_radix_3_base_transform(std::complex<COMPLEX_TYPE>* y)
{
    static const COMPLEX_TYPE k = -sin(2.0 * std::numbers::pi / 3.0);
    constexpr std::complex<COMPLEX_TYPE> J = std::complex<COMPLEX_TYPE>(0, 1);
    std::complex<COMPLEX_TYPE> a = y[1] + y[2];
    std::complex<COMPLEX_TYPE> t1 = y[0] - a / 2.0;
    std::complex<COMPLEX_TYPE> t2 = (y[1] - y[2]) * k;
    std::complex<COMPLEX_TYPE> b = J * t2;
    y[0] += a;
    y[1] = t1 + b;
    y[2] = t1 - b;
};

void fft_radix_5_base_transform(std::complex<COMPLEX_TYPE>* y)
{
    static const COMPLEX_TYPE k0 = cos(4.0 * std::numbers::pi / 5.0);
    static const COMPLEX_TYPE k10 = -sin(4.0 * std::numbers::pi / 5.0);
    static const COMPLEX_TYPE t = sin(2.0 * std::numbers::pi / 5.0);
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
}

template <size_t Radix>
void fft(std::vector<std::complex<COMPLEX_TYPE>>& a, std::function<void(std::complex<COMPLEX_TYPE>*)> base_transform)
{

    size_t len = x.size();

        for (size_t i = 0; i < len; i++)
        {
            base_transform(&x[Radix * i]);
        }

        size_t nl = len / 2;
        size_t nr = Radix;

        for (size_t i = 0; i < len; i++)
        {
            for (size_t l = 0; l < nl; l++)
            {
                for (size_t r = 0; r < nr; r++)
                {
                    size_t p = r + 2 * l * nr;
                    size_t q = p + nr;

                    std::complex<COMPLEX_TYPE> t = x[q] * lut.w[r * nl];

                    x[q] = x[p] - t;
                    x[p] += t;
                }
            }

            nl /= 2;
            nr *= 2;
        }
}


void dft(std::vector<std::complex<COMPLEX_TYPE>>& x)
{
   if (x.size() < 2)
    {
        return;
    }

    if (is_power2(x.size()))
    {
        return fft<2>(x, fft_radix_2_base_transform);
    }
    else if (x.size() % 3 == 0 && is_power(x.size() / 3))
    {
        return fft<3>(x, fft_radix_3_base_transform);
    }
    else if (x.size() % 5 == 0 && is_power2(x.size() / 5))
    {
        return fft<5>(x, fft_radix_5_base_transform);
    }
    else
    {
        return;
    }
}

void idft(std::vector<std::complex<COMPLEX_TYPE>>& x)
{
        

        auto pre_proc = [&x]
        (size_t begin, size_t end)
        {
            for (size_t i = begin; i < end; i++)
            {
                x[i] = std::conj(x[i]);
            }
        };

        size_t begin = 0;
        size_t end = per_cpu_block;

        std::vector<std::thread> pre_threads;
        for (size_t n = 0; n < n_cpu; n++)
        {
            pre_threads.push_back(std::thread(pre_proc, begin, end));

            begin = end;
            end = std::min(x.size(), begin + per_cpu_block);
        }

    dft(x);

}


