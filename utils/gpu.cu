#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>

#if _WIN32
#include <Windows.h>
#endif

//--------------------------------------------------------------------------------------------------
/// @brief Pause the application if run within its own console window.
void pause()
{
    bool pause = false;
    #if _WIN32
    HWND consoleWnd = GetConsoleWindow();
    DWORD dwProcessId;
    GetWindowThreadProcessId(consoleWnd, &dwProcessId);
    if (GetCurrentProcessId() == dwProcessId) pause = true;
    #endif

    if (pause)
    {
        std::cout << std::endl << "Press [ENTER] to exit..." << std::flush;
        std::cin.ignore(1);
    }
} // end pause

//--------------------------------------------------------------------------------------------------
/// @brief Macro to check for and report Cuda errors.
/// @param name - Name (or title) of the call being exeucted
/// @param call - The cuda function the macro should execute and check
#define CUDA_CHECK(name, call) \
{ \
    const cudaError_t r = (call); \
    std::cout << name << ": " << cudaGetErrorString(r) << std::endl; \
    if (r != cudaSuccess) \
    { \
        pause(); \
        return 1; \
    } \
} // end CUDA_CHECK

//--------------------------------------------------------------------------------------------------
int main()
{
    // determine the number of CUDA GPUs
    int count = 0;
    CUDA_CHECK("Device count", cudaGetDeviceCount(&count));
    std::cout << count << " CUDA devices available" << std::endl;

    // display stats on each GPU
    cudaDeviceProp prop;
    size_t memory;
    for (int i = 0; i < count; ++i)
    {
        CUDA_CHECK("Device properties", cudaGetDeviceProperties(&prop, i));
        CUDA_CHECK("Set GPU", cudaSetDevice(i));
        CUDA_CHECK("Memory info", cudaMemGetInfo(NULL, &memory));

        int cores = 0;
        switch (prop.major)
        {
        case 1: cores = prop.multiProcessorCount * 8; break;                            // Tesla (not supported starting CUDA 7.0)
        case 2: cores = prop.multiProcessorCount * (prop.minor == 0 ? 32 : 48); break;  // Fermi (not supported starting CUDA 9.2)
        case 3: cores = prop.multiProcessorCount * 192; break;                          // Kepler
        case 5: cores = prop.multiProcessorCount * 128; break;                          // Maxwell
        case 6: cores = prop.multiProcessorCount * (prop.minor == 0 ? 64 : 128); break; // Pascal
        case 7: cores = prop.multiProcessorCount * 64; break;                           // Volta
        }

        std::cout << "GPU #" << i << ": " << prop.name << std::endl
                  << "        Compute capability: " << prop.major << "." << prop.minor << std::endl
                  << "        Multiprocessors:    " << prop.multiProcessorCount << std::endl
                  << "        Cores:              " << cores << std::endl
                  << "        Clock rate:         " << prop.clockRate * 1.0e-6 << " GHz" << std::endl // KHz -> GHz
                  << "        Memory:             " << memory * 1e-9 << " GB" << std::endl
                  << "        Ratio 32 vs. 64:    " << prop.singleToDoublePrecisionPerfRatio << ":1" << std::endl
                  << "        maxThreads:	  " << prop.maxThreadsDim[0] << " " << prop.maxThreadsDim[1] << " " << prop.maxThreadsDim[2] << std::endl
		  << "        PCI Bus ID: 	  " << prop.pciBusID << std::endl;
    }

    // display driver version and application runtime version
    int driver, runtime;
    CUDA_CHECK("Driver version", cudaDriverGetVersion(&driver));
    CUDA_CHECK("Runtime version", cudaRuntimeGetVersion(&runtime));
    std::cout << "Driver version:  " << driver <<std::endl
              << "Runtime version: " << runtime << std::endl;

    pause();
    return 0;
}
