#ifndef SIMD_MATH
#define SIMD_MATH_
//

#ifdef __AVX512F__ 
#include <immintrin.h>
//
//
//

#ifdef __INTEL_COMPILER
#warning "simd_math AVX512F detected"
#define ADD(x,y) _mm512_add_pd(x, y)
#define SUB(x,y) _mm512_sub_pd(x, y)
#define MUL(x,y) _mm512_mul_pd(x, y)
#define SET(x)   _mm512_set1_pd(x)
//
inline __m512d operator + (__m512d a, __m512d b) {return _mm512_add_pd(a, b);}
inline __m512d operator - (__m512d a, __m512d b) {return _mm512_sub_pd(a, b);}
inline __m512d operator * (__m512d a, __m512d b) {return _mm512_mul_pd(a, b);}
inline __m512d operator / (__m512d a, __m512d b) {return _mm512_div_pd(a, b);}
#endif




#define __INV RCP_AVX512
//#define __INV RCP_1NR_AVX512
//#define __INV RCP_2NR_AVX512
//#define __INV _mm512_rcp28_pd
//#define __INV _mm512_rcp14_pd
#define __SQRT _mm512_sqrt_pd
//
//
inline __m512d RCP_AVX512(const __m512d d)
{
#warning "RCP_AVX512"
        __m512d x0;
        //x0 = _mm512_div_pd(_mm512_set1_pd(1.), d);
        //return x0;
         x0 = _mm512_rcp28_pd (d);
        //
        return x0;
}

//
inline __m512d RCP_1NR_AVX512(const __m512d d)
{
#warning "RCP_1NR_AVX512"

        __m512d x0;
        //x0 = _mm512_div_pd(_mm512_set1_pd(1.), d);
        //return x0;
         x0 = _mm512_rcp28_pd (d);
        //
        x0 = SUB(ADD(x0, x0), MUL(d, MUL(x0, x0)));
        //
        return x0;
}
//
//
inline __m512d RCP_2NR_AVX512(const __m512d d)
{
#warning "RCP_2NR_AVX512"

	__m512d x0;
	//x0 = _mm512_div_pd(_mm512_set1_pd(1.), d);
	//return x0;
         x0 = _mm512_rcp28_pd (d);
        //
        x0 = SUB(ADD(x0, x0), MUL(d, MUL(x0, x0)));
        x0 = SUB(ADD(x0, x0), MUL(d, MUL(x0, x0)));
        //
        return x0;
}
//
#endif
//
#endif
