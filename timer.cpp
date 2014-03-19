#if defined(_WIN32)
#include <windows.h>
	long double Timer()
	{
		LARGE_INTEGER ticksPerSecond;
		LARGE_INTEGER tick;
		QueryPerformanceFrequency(&ticksPerSecond);
		QueryPerformanceCounter(&tick);
		return (long double)tick.QuadPart/(long double)ticksPerSecond.QuadPart;
	}
#else
#include <sys/time.h>
#include <stddef.h>
	long double Timer()
	{
		long double t1;
		struct timeval time;
		gettimeofday(&(time), NULL);
		t1 =  (long double)time.tv_sec + (long double)time.tv_usec/(1000.0*1000.0);
		return t1;
	}
#endif
