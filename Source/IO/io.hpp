
#ifndef _IO_H
#define _IO_H	

/*
	This is the part of INMOST or INMOST library

	conversation was tested on intel core i7 (LittleEndian, IEEE) on windows (visual studio 2012), ubuntu (g++-4.6)?
	short, int, long; float, double, long double (8-byte)
	bluegene/p (BigEndian, IEEE) on aix (xlc++, v9.0)
	short, int, long; float, double, long double (16-byte)
	and between intel core i7 and bluegene/p

	due to absence of hardware:
		conversations for MiddleEndian were not tested
		conversations for IBM4, IBM8 and VAX4 were not tested

	VAX8 not implemented (probably wrong format test, no algorithm in convert_bytes_to_float)
	IEEE16, VAX16,CRAY4,CRAY8,CRAY16 not implemented (no format test, no algorithm in convert_bytes_to_float)

	unknown formats with   4 bytes will be converted to IEEE4
	unknown formats with >=8 bytes will be converted to IEEE8 (explicitly set in write_fValue, read_fValue)

	it's supposed that integers and floats have the same endianess
*/

#include <math.h>
#include <ostream>
#include <istream>
#include <vector>
#include <stdint.h>
namespace INMOST
{
	
#define INT_CONST                   "\x01\x02\x03\x04"
#define LITTLE_ENDIAN_BYTE_ORDER    0x04030201
#define MIDDLE_ENDIAN_BYTE_ORDER    0x02010403
#define BIG_ENDIAN_BYTE_ORDER       0x01020304

#define FLOAT_CONST -10.0 //use const that have precise representation in memory
#define IEEE_SINGLE "\x00\x00\x20\xC1"
#define IBM_SINGLE  "\x00\x00\xA0\xC1"
#define VAX_SINGLE  "\x00\x00\x20\xC2" //not tested
#define CRAY_SINGLE "\x00\x00\x00\x00" //not implemented


#define IEEE_DOUBLE "\x00\x00\x00\x00\x00\x00\x24\xC0"
#define IBM_DOUBLE  "\x00\x00\x00\x00\x00\x00\xA0\xC1"
#define VAX_DOUBLE  "\x00\x00\x00\x00\x00\x00\x20\xC2" //not tested
#define CRAY_DOUBLE "\x00\x00\x00\x00\x00\x00\x00\x00" //not implemented
//long double?

#define IEEE_QUAD "\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x24\xC0"
#define IBM_QUAD  "\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xA0\xC1"
#define VAX_QUAD  "\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x20\xC2" //not tested
#define CRAY_QUAD "\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00" //not implemented

	template<typename iType, typename fType>
	class io_converter
	{
	public:
		enum iByteOrder {LittleEndian = 0x01,BigEndian = 0x02,MiddleEndian = 0x04, iUnknown = 0x00}; 
		enum fByteOrder {IEEE = 0x01, IBM = 0x02, VAX = 0x03, CRAY = 0x04, fUnknown = 0x00};
	private:
		iByteOrder local_iorder, source_iorder;
		fByteOrder local_forder, source_forder;
		unsigned char local_ibytes, local_fbytes;
		unsigned char source_ibytes, source_fbytes;
		unsigned char temp[255];
		void flip_bytes(unsigned char size)
		{
			for (unsigned char i = 0; i < (size >> 1); i++)
			{
				unsigned char t = temp[i];
				temp[i] = temp[size-1-i];
				temp[size-1-i] = t;
			}
		}
		template<typename float_type>
		fType extract_double()
		{
			const uint64_t a0 = 0x7FFFFFFFFFFFFFFFull;
			const uint64_t a1 = 0x000FFFFFFFFFFFFFull;
			const uint64_t a2 = 0x0010000000000000ull;
			const uint64_t a3 = 0x7FF0000000000000ull;
			const uint64_t a4 = 0x8000000000000000ull;
			uint64_t * tmp = reinterpret_cast<uint64_t *>(&temp[0]);
			float_type f = 0.0;
			if ( ((*tmp) & a0) != 0 ) 
				f = static_cast<float_type>(ldexp( static_cast<long double>(((*tmp) & a1) | a2), static_cast<int>( ((*tmp) & a3) >> 52 )- 1022 - 53 ));
			if ( ((*tmp) & a4) != 0 )
				f = -f;
			return static_cast<fType>(f);
		}
		template<typename float_type>
		fType extract_float()
		{
			uint32_t * tmp = reinterpret_cast<uint32_t *>(&temp[0]);
			float_type f = 0.0;
			if ( ((*tmp) & 0x7FFFFFFF) != 0 ) 
				f = static_cast<float_type>(ldexp( static_cast<long double>(((*tmp) & 0x007FFFFF) | 0x00800000), static_cast<int>( ((*tmp) & 0x7F800000) >> 23 )- 126 - 24 ));
			if ( ((*tmp) & 0x80000000) != 0 )
				f = -f;
			return static_cast<fType>(f);
		}
		void convert_endianess(unsigned char size,  iByteOrder input,  iByteOrder output)
		{
			if( input != output )
			{
				if( (input & (BigEndian | LittleEndian)) && (output & (BigEndian | LittleEndian) ) )
					flip_bytes(size);
				else if( input & MiddleEndian )
				{
					for(unsigned char k = 0; k < size; k+=2) // to Big-Endian
					{
						unsigned char t = temp[k];
						temp[k] = temp[k+1];
						temp[k+1] = t;
					}
					if( output == LittleEndian ) flip_bytes(size);
				}
				else if( output & MiddleEndian )
				{
					if( input == LittleEndian ) flip_bytes(size); // to Big-Endian
					for(unsigned char k = 0; k < size; k+=2) // to Middle-Endian
					{
						unsigned char t = temp[k];
						temp[k] = temp[k+1];
						temp[k+1] = t;
					}
				}
			}
		}
		template<typename float_type>
		fType convert_bytes_to_float(unsigned char * fbytes, unsigned char size,   fByteOrder forder)
		{
			float_type ret = 0;
			if( size == 4 )
			{
				unsigned char S; int E; unsigned long F;
				unsigned char b1 = fbytes[3],b2 = fbytes[2],b3  = fbytes[1], b4 = fbytes[0];
				float_type M, F1, A, B, C, D2, e23, e24;
				switch(forder)
				{
				case IEEE:
					S = ( b1 & 0x80 ) >> 7;
					E = ( ( b1 & 0x7f ) << 1 ) + ( ( b2 & 0x80 ) >> 7 );
					F = ( ( b2 & 0x7f ) << 16 ) + ( b3 << 8 ) + b4;
					A = 2.0;
					B = 127.0;
					C = 1.0;	
					D2 = -126.0;
					e23 = 8388608.0;		// 2^23
					M = (float_type) F / e23;
					if ( S == 0) F1 = 1.0; else F1 = -1.0;		
					if ( 0 < E && E < 255 ) ret = F1 * ( C + M ) * ::pow ( A, E - B ) ;	
					else if ( E == 0 && F != 0 ) ret = F1 * M * ::pow ( A, D2 );
					else if ( E == 0 && F == 0 && S == 1 ) ret = -0;
    				else if ( E == 0 && F == 0 && S == 0 ) ret = 0;	
					else if ( E == 255 && F != 0) ret = 0; // Not a number
					else if ( E == 255 && F == 0 && S == 1 ) ret = 0; // -Infinity
					else if ( E == 255 && F == 0 && S == 0 ) ret = 0; // Infinity	
					break;
				case IBM:
					S = ( b1 & 0x80 ) >> 7;
					E = ( b1 & 0x7f );
					F = ( b2 << 16 ) + ( b3 << 8 ) + b4;
					A = 16.0;
					B = 64.0;
					//D = static_cast<float_type>(0.69314718055994529);	// log2
					e24 = 16777216.0;		// 2^24
					M = (float_type) F / e24;
					if ( S == 0) F1 = 1.0;
					else F1 = -1.0;
					if ( S == 0 && E == 0 && F == 0 ) ret = 0;
					else ret = F1 * M * ::pow ( A, E - B ) ;
					break;
				case VAX:
					S = ( b1 & 0x80 ) >> 7;
					E = ( ( b1 & 0x7f ) << 1 ) + ( ( b2 & 0x80 ) >> 7 );
					F = ( ( b2 & 0x7f ) << 16 ) + ( b3 << 8 ) + b4;
					A = 2.0;
					B = 128.0;
					C = 0.5;	
					e24 = 16777216.0;		// 2^24
					M = (float_type) F / e24;	
					if ( S == 0 ) F1 = 1.0;
					else F1 = -1.0;
					if ( 0 < E ) ret = F1 * ( C + M ) * ::pow ( A, E - B ) ;	
                    else if ( E == 0 && S == 0 ) ret = 0;
					else if ( E == 0 && S == 1 ) ret = 0; // reserved
					break;
				default: throw "not implemented";
				}
			}
			else if( size == 8 )
			{
				unsigned char S; int E; unsigned long L1,L2;
				unsigned char b1 = fbytes[7],b2 = fbytes[6],b3 = fbytes[5], b4 = fbytes[4];
				unsigned char b5 = fbytes[3],b6 = fbytes[2],b7 = fbytes[1],b8 = fbytes[0];
				float_type F1, M, M1, M2, A, B, C, D1, D2, e20, e24, e52, e56;
				switch(forder)
				{
				case IEEE:
					S = ( b1 & 0x80 ) >> 7;
					E = ( ( b1 & 0x7f ) << 4 ) + ( ( b2 & 0xf0 ) >> 4 );
					L1 = ( ( b2 & 0x0f ) << 16 ) + ( b3 << 8 ) + b4;
					L2 = ( b5 << 24 ) + ( b6 << 16 ) + ( b7 << 8 ) + b8;
					e20 = 1048576.0;			// 2^20
					e52 = 4503599627370496.0;	// 2^52
					A = 2.0;
					B = 1023.0;
					C = 1.0;
					D1 = 2047;
					D2 = 1022.0;
					M1 = (float_type) L1 / e20 ;
					M2 = (float_type) L2 / e52 ;
					M = M1 + M2 ;
					if ( S == 0) F1 = 1.0;
					else F1 = -1.0;
					if ( 0 < E && E < D1 ) ret = F1 * ( C + M ) * ::pow ( A, E - B ) ;		
					else if ( E == 0 && M != 0 ) ret = F1 * M * ::pow ( A, D2 );
					else if ( E == 0 && M == 1 ) ret = F1 * 0;
					else if ( E == D1 && M == 0 ) ret = 0;
					else if ( E == D1 && M != 0 ) ret = 0; // Not a number
					break;
				case IBM:
					S = ( b1 & 0x80 ) >> 7;
					E = ( b1 & 0x7f );
					L1 = ( b2 << 16 ) + ( b3 << 8 ) + b4;
					L2 = ( b5 << 24 ) + ( b6 << 16 ) + ( b7 << 8 ) + b8;
					A = 16.0;
					B = 64.0;
					e24 = 16777216.0;			    // 2^24
					e56 = 72057594037927936.0;		// 2^56
					M1 = (float_type) L1 / e24 ;
					M2 = (float_type) L2 / e56 ;
					M = M1 + M2 ;
					if ( S == 0) F1 = 1.0;
					else F1 = -1.0;
					if ( S == 0 && E == 0 && M == 0 ) ret = 0;
					else ret = F1 * M * ::pow ( A, E - B );
					break;
				default: throw NotImplemented;
				}
			}
            else
            {
                std::cout << "Converter for float of size " << size << " do not exist" << std::endl;
                throw NotImplemented;
            }
			return static_cast<fType>(ret);
		}
		fType load_fp_test_bytes(const char * bytes, unsigned char size)
		{
			memcpy(&temp[0],bytes,size);
			convert_endianess(size,LittleEndian,local_iorder);
			//flip_bytes(size);
			return (*(fType *)(&temp[0]));
		}
	public:
		io_converter()
		{
			local_iorder = get_iByteOrder();
			local_forder = get_fByteOrder();
			local_ibytes = get_iByteSize();
			local_fbytes = get_fByteSize();
			source_iorder = iUnknown;
			source_forder = fUnknown;
			source_ibytes = 255;
			source_fbytes = 255;
		}
		io_converter(const io_converter & other)
		{
			local_iorder = get_iByteOrder();
			local_forder = get_fByteOrder();
			local_ibytes = get_iByteSize();
			local_fbytes = get_fByteSize();
			source_iorder = other.source_iorder;
			source_forder = other.source_forder;
			source_ibytes = other.source_ibytes;
			source_fbytes = other.source_fbytes;
		}
		io_converter & operator = (io_converter const & other)
		{
			local_iorder = get_iByteOrder();
			local_forder = get_fByteOrder();
			local_ibytes = get_iByteSize();
			local_fbytes = get_fByteSize();
			source_iorder = other.source_iorder;
			source_forder = other.source_forder;
			source_ibytes = other.source_ibytes;
			source_fbytes = other.source_fbytes;
			return *this;
		}


		unsigned char get_iByteSize() {return sizeof(iType);}
		unsigned char get_fByteSize() {return sizeof(fType);}
		unsigned char get_source_iByteSize() {return source_ibytes;}
		unsigned char get_source_fByteSize() {return source_fbytes;}
		iByteOrder get_source_iByteOrder() {return source_iorder;}
		fByteOrder get_source_fByteOrder() {return source_forder;}
		iByteOrder get_iByteOrder()
		{
			if( (*(uint32_t *)INT_CONST) == LITTLE_ENDIAN_BYTE_ORDER)
				return LittleEndian;
			else if( (*(uint32_t *)INT_CONST) == BIG_ENDIAN_BYTE_ORDER)
				return BigEndian;
			else if( (*(uint32_t *)INT_CONST) == MIDDLE_ENDIAN_BYTE_ORDER)
				return MiddleEndian;
			else
				return iUnknown;
		}
		fByteOrder get_fByteOrder()
		{
			fType test;
			if( get_fByteSize() == 4 )
			{
				test = load_fp_test_bytes(IEEE_SINGLE,4);
				if( test == FLOAT_CONST )
					return IEEE;
				test = load_fp_test_bytes(IBM_SINGLE,4);
				if( test == FLOAT_CONST )
					return IBM;
				test = load_fp_test_bytes(VAX_SINGLE,4);
				if( test == FLOAT_CONST )
					return VAX;
				test = load_fp_test_bytes(CRAY_SINGLE,4);
				if( test == FLOAT_CONST )
					return CRAY;
			}
			else if( get_fByteSize() == 8 )
			{
				test = load_fp_test_bytes(IEEE_DOUBLE,8);
				if( test == FLOAT_CONST )
					return IEEE;
				test = load_fp_test_bytes(IBM_DOUBLE,8);
				if( test == FLOAT_CONST )
					return IBM;
				test = load_fp_test_bytes(VAX_DOUBLE,8);
				if( test == FLOAT_CONST )
					return VAX;
				test = load_fp_test_bytes(CRAY_DOUBLE,8);
				if( test == FLOAT_CONST )
					return CRAY;
			}
			return fUnknown;
		}
		static const char * str_iByteOrder(iByteOrder order)
		{
			switch(order)
			{
			case LittleEndian: return "Little-endian";
			case BigEndian: return "Big-endian";
			case MiddleEndian: return "Middle-endian";
			default: return "Unknown";
			}
		}
		static const char * str_fByteOrder(fByteOrder order)
		{
			switch(order)
			{
			case IEEE: return "IEEE";
			case IBM: return "IBM";
			case VAX: return "VAX";
			case CRAY: return "CRAY";
			default: return "Unknown";
			}
		}
		std::ostream& write_iByteOrder(std::ostream& dest)
		{
			unsigned char c;
			switch(local_iorder)
			{
				case LittleEndian: c = 0x01; break;
				case BigEndian:    c = 0x02; break;
				case MiddleEndian: c = 0x03; break;
				default:           c = 0x00; break;
			}
			dest.put(c);
			return dest;
		}
		std::ostream& write_iByteSize(std::ostream& dest)
		{
			dest.put(local_ibytes);
			return dest;
		}
		std::ostream& write_fByteOrder(std::ostream& dest)
		{
			unsigned char c;
			switch(local_forder)
			{
				case IEEE: c = 0x01; break;
				case IBM:  c = 0x02; break;
				case VAX:  c = 0x03; break;
				case CRAY:  c = 0x04; break;
				default: c = 0x00; break;
			}
			dest.put(c);
			return dest;
		}
		std::ostream& write_fByteSize(std::ostream& dest)
		{
			dest.put(local_fbytes);
			return dest;
		}
		std::istream& read_iByteOrder( std::istream& source)
		{
			char c;
			source.get(c);
			switch(c)
			{
				case 0x01: source_iorder = LittleEndian; break;
				case 0x02: source_iorder = BigEndian;    break;
				case 0x03: source_iorder = MiddleEndian;    break;
				default: source_iorder = iUnknown;      break;
			}
			return source;
		}
		std::istream& read_iByteSize(std::istream& source)
		{
			source.get(reinterpret_cast<char &>(source_ibytes));
			return source;
		}
		std::istream& read_fByteOrder( std::istream& source)
		{
			char c;
			source.get(c);
			switch(c)
			{
				case 0x01: source_forder = IEEE; break;
				case 0x02: source_forder = IBM; break;
				case 0x03: source_forder = VAX; break;
				case 0x04: source_forder = CRAY; break;
				default: source_forder = fUnknown;      break;
			}
			return source;
		}
		std::istream& read_fByteSize(std::istream& source)
		{
			source.get(reinterpret_cast<char &>(source_fbytes));
			return source;
		}
		void set_iBytesSize(unsigned char size) {source_ibytes = size;}
		void set_fBytesSize(unsigned char size) {source_fbytes = size;}
		void set_iBytesOrder(iByteOrder order) {source_iorder = order;}
		void set_fBytesOrder(fByteOrder order) {source_forder = order;}

		
		std::ostream & write_iValue(std::ostream & dest, iType value)
		{

			dest.write(reinterpret_cast<const char *>(&value),local_ibytes);
			return dest;
		}

		std::istream & read_iValue(std::istream & source, iType & value)
		{
			unsigned char * bytes = reinterpret_cast<unsigned char *>(&value);
			memset(bytes,0,local_ibytes);
			source.read(reinterpret_cast<char *>(&temp[0]),source_ibytes); // read all bytes to temporary place
			convert_endianess(source_ibytes,source_iorder,local_iorder);
			unsigned char min_ibytes = std::min(local_ibytes,source_ibytes);
			if( local_iorder & LittleEndian )
				for(unsigned char i = 0; i < min_ibytes; i++) bytes[i] = temp[i]; //copy bytes to output
			else //that should be fine for middle-endian
				for(unsigned char i = source_ibytes-1; i >= source_ibytes-min_ibytes; i--) bytes[i+local_ibytes-source_ibytes] = temp[i]; //copy bytes to output
			return source;
		}
		
		std::ostream & write_fValue(std::ostream& dest, fType value)
		{
			switch(local_forder)
			{
			case IEEE:
			case IBM:
			case VAX:
			//case CRAY:
				dest.write(reinterpret_cast<const char *>(&value),local_fbytes);
				break;
			default: //store in ieee format with little-endian order
				{
					//if( local_ibytes == 16 ) //cray
					if( local_fbytes >= 8 )
					{
						bool isNeg = value < 0;
						if ( isNeg ) value = - value;
						int exp;
						if ( value == 0.0 ) exp = 0;
						else 
						{
							value = ldexp( frexp( value, &exp ), 53 );
							exp += 1022;
						}
						uint64_t mant = static_cast<uint64_t>( value );
						temp[7] = static_cast<unsigned char>( (isNeg ? 0x80 : 0x00) | exp >> 4 );
						temp[6] = static_cast<unsigned char>(((exp << 4) & 0xF0) | ((mant >> 48) & 0x0F));
						temp[5] = static_cast<unsigned char>(mant >> 40);
						temp[4] = static_cast<unsigned char>(mant >> 32);
						temp[3] = static_cast<unsigned char>(mant >> 24);
						temp[2] = static_cast<unsigned char>(mant >> 16);
						temp[1] = static_cast<unsigned char>(mant >> 8 );
						temp[0] = static_cast<unsigned char>(mant      );
						dest.write(reinterpret_cast<const char *>(&temp[0]),8);

					}
					else if( local_fbytes == 4 )
					{
						bool isNeg = value < 0;
						if ( isNeg ) value = - value;
						int exp;
						if ( value == 0.0 ) exp = 0;
						else 
						{
							value = ldexp( frexp( value, &exp ), 24 );
							exp += 126;
						}
						uint32_t mant = static_cast< uint32_t >( value );
						temp[3] =  (isNeg ? 0x80 : 0x00) | exp >> 1;
						temp[2] = ((exp << 7) & 0x80) | ((mant >> 16) & 0x7F) ;
						temp[1] = ( mant >> 8 );
						temp[0] = ( mant );
						dest.write(reinterpret_cast<const char *>(&temp[0]),4);
					}
				}
			}
			return dest;
		}
	
	
		


		std::istream& read_fValue(std::istream& source, fType & value)
		{
			switch(source_forder)
			{
			case IEEE:
			case IBM:
			case VAX:
			//case CRAY:
				{
					source.read(reinterpret_cast<char *>(&temp[0]),source_fbytes);
					if( source_forder == local_forder && source_fbytes == local_fbytes )
					{
						convert_endianess(source_fbytes,source_iorder, local_iorder);
						value = *reinterpret_cast<fType *>(&temp[0]);
					}
					else
					{
						convert_endianess(source_fbytes,source_iorder, LittleEndian); //next algorithm expect little-endian byte order
						if( sizeof(float) == source_fbytes ) value = convert_bytes_to_float<float>(&temp[0],source_fbytes,source_forder);
						else if( sizeof(double) == source_fbytes ) value = convert_bytes_to_float<double>(&temp[0],source_fbytes,source_forder);
						else if( sizeof(long double) == source_fbytes ) value = convert_bytes_to_float<long double>(&temp[0],source_fbytes,source_forder);
					}
				}
				break;
			default:
				{
					unsigned char read_fbytes = std::min(static_cast<unsigned char>(8),source_fbytes);
					source.read(reinterpret_cast<char *>(&temp[0]),read_fbytes);
					convert_endianess(read_fbytes,LittleEndian, local_iorder); //we always save in little-endian order
					if( read_fbytes == 8 )
					{
						if( sizeof(float) == 8 ) value = extract_double<float>();
						else if( sizeof(double) == 8 ) value = extract_double<double>();
						else if( sizeof(long double) == 8 ) value = extract_double<long double>();
						else throw "No matching datatype";
					}
					else if( source_fbytes == 4 )
					{
						if( sizeof(float) == 4 ) value = extract_float<float>();
						else if( sizeof(double) == 4 ) value = extract_float<double>();
						else if( sizeof(long double) == 4 ) value = extract_float<long double>();
						else throw "No matching datatype";
					}
				}
			}
			return source;
		}
	};
	
}
#endif //_IO_H
