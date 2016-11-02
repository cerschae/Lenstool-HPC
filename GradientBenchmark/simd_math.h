#ifndef SIMD_MATH
#define SIMD_MATH_
//
#include <immintrin.h>
//
inline __m256d RCP(const __m256d d)
{
        const __m128 b    = _mm256_cvtpd_ps(d);
        const __m128 rcp = _mm_rcp_ps (b);
        __m256d x0       = _mm256_cvtps_pd(rcp);
        //
        return x0;
}
//
//
//
inline __m256d RCP_1NR(const __m256d d)
{
        const __m128 b    = _mm256_cvtpd_ps(d);
        const __m128 rcp = _mm_rcp_ps (b);
        __m256d x0       = _mm256_cvtps_pd(rcp);
        //
        x0 = x0 + x0 - d*x0*x0;
        //                        //
        return x0;
}
//
//
//
inline __m256d RCP_2NR(const __m256d d)
{
        const __m128 b    = _mm256_cvtpd_ps(d);
        const __m128 rcp = _mm_rcp_ps (b);
        __m256d x0       = _mm256_cvtps_pd(rcp);
        //
        x0 = x0 + x0 - d*x0*x0;
        x0 = x0 + x0 - d*x0*x0;
        //
        return x0;
}


inline __m256d SQRT(const __m256d d)
{
        const __m128 b    = _mm256_cvtpd_ps(d);
        const __m128 rcp = _mm_sqrt_ps (b);
        __m256d x0       = _mm256_cvtps_pd(rcp);

        return x0;
}
//
//
//
inline __m256d SQRT_1NR(const __m256d d)
{
        const __m128 b   = _mm256_cvtpd_ps(d);
        const __m128 rcp = _mm_sqrt_ps (b);
        __m256d x0       = _mm256_cvtps_pd(rcp);

        __m256d half     = _mm256_set1_pd(0.5);
        __m256d three    = _mm256_set1_pd(3.);

        __m256d a        = RCP_1NR(d);

        x0 = half*x0*(three - x0*x0*a);

        return x0;
}
//
//
//
inline __m256d SQRT_2NR(const __m256d d)
{
        const __m128 b   = _mm256_cvtpd_ps(d);
        const __m128 rcp = _mm_sqrt_ps (b);
        __m256d x0       = _mm256_cvtps_pd(rcp);

        __m256d half     = _mm256_set1_pd(0.5);
        __m256d three    = _mm256_set1_pd(3.);

        __m256d a        = RCP_2NR(d);

        x0 = half*x0*(three - x0*x0*a);
        x0 = half*x0*(three - x0*x0*a);

        return x0;
}




#endif
