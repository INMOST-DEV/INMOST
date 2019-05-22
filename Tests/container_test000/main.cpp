#include "inmost.h"
#include <cstdio>
#include <cmath>
using namespace INMOST;

int main(int argc,char ** argv)
{
	int test = 0;
	if (argc > 1)  test = atoi(argv[1]);

	
	if( test == 0 ) //check order of deallocation does not cause problem
	{
		std::cout << "allocs: " << get_pool().allocations() << " byte " << get_pool().last_byte() << std::endl;
		pool_array<double> * a = new pool_array<double>(2);
		std::cout << "allocs: " << get_pool().allocations() << " byte " << get_pool().last_byte() << std::endl;
		pool_array<int> * b = new pool_array<int>(3);
		std::cout << "allocs: " << get_pool().allocations() << " byte " << get_pool().last_byte() << std::endl;
		pool_array<char> * c = new pool_array<char>(6);
		std::cout << "allocs: " << get_pool().allocations() << " byte " << get_pool().last_byte() << std::endl;
		delete a;
		std::cout << "allocs: " << get_pool().allocations() << " byte " << get_pool().last_byte() << std::endl;
		delete b;
		std::cout << "allocs: " << get_pool().allocations() << " byte " << get_pool().last_byte() << std::endl;
		delete c;
		std::cout << "allocs: " << get_pool().allocations() << " byte " << get_pool().last_byte() << std::endl;
	}
	return 0;
}
