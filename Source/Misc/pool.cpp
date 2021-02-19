#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "inmost.h"

namespace INMOST
{
	static thread_private<memory_pool> _pool;


	memory_pool & get_pool()
	{
		//std::cout << "get pool " << (void *)&(*_pool) << std::endl;
		return *_pool;
	}
}
