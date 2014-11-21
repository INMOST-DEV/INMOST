
#if defined(_WIN32)
#include <windows.h>
	double Timer()
	{
		LARGE_INTEGER ticksPerSecond;
		LARGE_INTEGER tick;
		QueryPerformanceFrequency(&ticksPerSecond);
		QueryPerformanceCounter(&tick);
		return (double)tick.QuadPart/(double)ticksPerSecond.QuadPart;
	}
#else
#include <sys/time.h>
#include <stddef.h>
	long double Timer()
	{
		double t1;
		struct timeval time;
		gettimeofday(&(time), NULL);
		t1 =  (double)time.tv_sec + (double)time.tv_usec/(1000.0*1000.0);
		return t1;
	}
#endif
