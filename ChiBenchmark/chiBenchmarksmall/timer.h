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


static inline unsigned long long cycles()
{
	unsigned long long u;
	__asm__ volatile ("rdtscp;shlq $32,%%rdx;orq %%rdx,%%rax;movq %%rax,%0":"=q"(u)::"%rax", "%rdx", "rcx");
	return u;
}
