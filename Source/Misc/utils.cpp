#include "utils.h"
#include <cstdlib>
#if defined(USE_ZLIB)
#include "zlib.h"
#endif
// temporary fix for GeRa

bool IsMaster() {
    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank == 0;
}

namespace INMOST {

    template<>
    int from_string<int>(std::string str) {
        return atoi(str.c_str());
    }

    template<>
    double from_string<double>(std::string str) {
        return atof(str.c_str());
    }

    template<>
    std::string from_string<std::string>(std::string str) {
        return str;
    }

    std::string string_to_lower(const std::string &str) {
        std::string lower = std::string(str);
        for (size_t    i     = 0; i < lower.length(); i++)
            lower[i] = (char) tolower(lower[i]);
        return lower;
    }

    void MPIBarrier() {
#if defined(USE_MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    int MPIGetRank() {
        int rank = 0;
#if defined(USE_MPI)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        return rank;
    }

	void split_string(const std::string & str, std::vector<std::string> & str_list, const std::string& delims)
	{
		std::string::size_type end = str.find_first_not_of(delims, 0);
		std::string::size_type pos = str.find_first_of(delims, end);
		while (std::string::npos != pos || std::string::npos != end) 
		{
			str_list.push_back(str.substr(end, pos - end));
			end = str.find_first_not_of(delims, pos);
			pos = str.find_first_of(delims, end);
		}
	}


	INMOST_DATA_REAL_TYPE Integrate(INMOST_DATA_REAL_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_REAL_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_REAL_TYPE, MPI_SUM, comm);
#else//USE_MPI
		(void)input;
#endif//USE_MPI
		return output;
	}
	INMOST_DATA_ENUM_TYPE Integrate(INMOST_DATA_ENUM_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_ENUM_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_SUM, comm);
#else//USE_MPI
		(void)input;
#endif//USE_MPI
		return output;
	}
	INMOST_DATA_INTEGER_TYPE Integrate(INMOST_DATA_INTEGER_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_INTEGER_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_INTEGER_TYPE, MPI_SUM, comm);
#else//USE_MPI
		(void)input;
#endif//USE_MPI
		return output;
	}
	void Integrate(INMOST_DATA_REAL_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_REAL_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_REAL_TYPE, MPI_SUM, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI

	}
	void Integrate(INMOST_DATA_ENUM_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_ENUM_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_ENUM_TYPE, MPI_SUM, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI
	}
	void Integrate(INMOST_DATA_INTEGER_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_INTEGER_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_INTEGER_TYPE, MPI_SUM, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI
	}
	INMOST_DATA_INTEGER_TYPE ExclusiveSum(INMOST_DATA_INTEGER_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_INTEGER_TYPE output = 0;
#if defined(USE_MPI)
		MPI_Scan(&input, &output, 1, INMOST_MPI_DATA_INTEGER_TYPE, MPI_SUM, comm);
		output -= input;
#else//USE_MPI
		(void)input;
#endif//USE_MPI
		return output;
	}
	INMOST_DATA_ENUM_TYPE ExclusiveSum(INMOST_DATA_ENUM_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_ENUM_TYPE output = 0;
#if defined(USE_MPI)
		MPI_Scan(&input, &output, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_SUM, comm);
		output -= input;
#else//USE_MPI
		(void)input;
#endif//USE_MPI
		return output;
	}
	INMOST_DATA_REAL_TYPE AggregateMax(INMOST_DATA_REAL_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_REAL_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_REAL_TYPE, MPI_MAX, comm);
#else //USE_MPI
		(void)input;
#endif //USE_MPI
		return output;
	}
	INMOST_DATA_ENUM_TYPE AggregateMax(INMOST_DATA_ENUM_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_ENUM_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_MAX, comm);
#else //USE_MPI
		(void)input;
#endif //USE_MPI
		return output;
	}
	INMOST_DATA_INTEGER_TYPE AggregateMax(INMOST_DATA_INTEGER_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_INTEGER_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_INTEGER_TYPE, MPI_MAX, comm);
#else //USE_MPI
		(void)input;
#endif //USE_MPI
		return output;
	}
	void AggregateMax(INMOST_DATA_REAL_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_REAL_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_REAL_TYPE, MPI_MAX, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI
	}
	void AggregateMax(INMOST_DATA_ENUM_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_ENUM_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_ENUM_TYPE, MPI_MAX, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI
	}
	void AggregateMax(INMOST_DATA_INTEGER_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_INTEGER_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_INTEGER_TYPE, MPI_MAX, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI
	}
	INMOST_DATA_REAL_TYPE AggregateMin(INMOST_DATA_REAL_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_REAL_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_REAL_TYPE, MPI_MIN, comm);
#else //USE_MPI
		(void)input;
#endif //USE_MPI
		return output;
	}
	INMOST_DATA_ENUM_TYPE AggregateMin(INMOST_DATA_ENUM_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_ENUM_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_MIN, comm);
#else //USE_MPI
		(void)input;
#endif //USE_MPI
		return output;
	}
	INMOST_DATA_INTEGER_TYPE AggregateMin(INMOST_DATA_INTEGER_TYPE input, INMOST_MPI_Comm comm)
	{
		INMOST_DATA_INTEGER_TYPE output = input;
#if defined(USE_MPI)
		MPI_Allreduce(&input, &output, 1, INMOST_MPI_DATA_INTEGER_TYPE, MPI_MIN, comm);
#else //USE_MPI
		(void)input;
#endif //USE_MPI
		return output;
	}
	void AggregateMin(INMOST_DATA_REAL_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_REAL_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_REAL_TYPE, MPI_MIN, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI
	}
	void AggregateMin(INMOST_DATA_ENUM_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_ENUM_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_ENUM_TYPE, MPI_MIN, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI

	}
	void AggregateMin(INMOST_DATA_INTEGER_TYPE* input, INMOST_DATA_ENUM_TYPE size, INMOST_MPI_Comm comm)
	{
#if defined(USE_MPI)
		std::vector<INMOST_DATA_INTEGER_TYPE> temp(size);
		if (!temp.empty())
		{
			std::copy(input, input + size, temp.begin());
			MPI_Allreduce(&temp[0], input, size, INMOST_MPI_DATA_INTEGER_TYPE, MPI_MIN, comm);
		}
#else//USE_MPI
		(void)input;
		(void)size;
#endif//USE_MPI

	}
	bool zcompress(const void* buffer, size_t dsize, void*& buffer_out, size_t& dsize_out)
	{
#if defined(USE_ZLIB)
		uLongf bufsize = compressBound((uLongf)dsize);
		void* zbuffer = malloc(bufsize);
		if (zbuffer)
		{
			int res = compress2(static_cast<Bytef*>(zbuffer), &bufsize, static_cast<const Bytef*>(buffer), (uLongf)dsize, 9);
			if (res != Z_OK)
				std::cout << __FILE__ << ":" << __LINE__ << " fail " << res << std::endl;
			else
			{
				buffer_out = zbuffer;
				dsize_out = bufsize;
				return true;
			}
		}
		else std::cout << __FILE__ << ":" << __LINE__ << " allocation of " << bufsize << " failed" << std::endl;
		return false;
#else
		return false;
#endif
	}
	bool zuncompress(const void* zbuffer, size_t zdsize, void* buffer_out, size_t& dsize_out)
	{
#if defined(USE_ZLIB)
		uLongf outsize = (uLongf)dsize_out, insize = (uLongf)zdsize;
		int res = uncompress(static_cast<Bytef*>(buffer_out), &outsize, static_cast<const Bytef*>(zbuffer), insize);
		if (res != Z_OK)
		{
			std::cout << __FILE__ << ":" << __LINE__ << " uncompress fail " << res << std::endl;
			return false;
		}
		if (dsize_out != outsize)
			std::cout << "Unexpected uncompressed data size: " << outsize << " expected: " << dsize_out << std::endl;
		return true;
#else
		return false;
#endif
	}
}
