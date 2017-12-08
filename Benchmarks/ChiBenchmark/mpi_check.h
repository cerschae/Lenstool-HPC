#pragma once

#include <mpi.h>

#ifndef MPI_CHECK
#define MPI_CHECK(stmt)							\
    do									\
    {									\
	const int code = stmt;						\
									\
	if (code != MPI_SUCCESS)					\
	{								\
	    char error_string[2048];					\
	    int length_of_error_string = sizeof(error_string);		\
	    MPI_Error_string(code, error_string, &length_of_error_string); \
	    								\
	    fprintf(stderr,						\
		    "ERROR!\n" #stmt " mpiAssert: %s %d %s\n",		\
		    __FILE__, __LINE__, error_string);			\
	    fflush(stderr);						\
	    								\
	    MPI_Abort(MPI_COMM_WORLD, code);				\
	}								\
    }									\
    while(0)
#endif
