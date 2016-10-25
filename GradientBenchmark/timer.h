#pragma once
//
#include <sys/time.h>
extern "C"{
double myseconds()
{
        struct timeval  tp;
        struct timezone tzp;
        //
        int i = gettimeofday(&tp,&tzp);
        //
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
}
