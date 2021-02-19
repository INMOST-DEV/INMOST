#ifndef _CONTAINER_HPP
#define _CONTAINER_HPP

#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <cmath>
#include <cassert>
#include <memory>
#include <new>
#include <iterator>
#include <limits>
#include "inmost_common.h"

//#define DEBUGMEM

#define NEW_CHUNKS
#define __VOLATILE volatile
//#define OUT_OF_RANGE

#if defined(_OPENMP)
#include <omp.h>
#endif

//TODO
// 1. change to uniform size_type instead of size_t, make it INMOST_DATA_ENUM_TYPE
/*
template<class element, class T1> struct isInputRandomIterators
{
	static void constraints(T1 a, T1 b) { ++a; (void)a++; a==a; a!=a; a-b; }
	isInputRandomIterators() { void(*p)(T1,T1) = constraints; (void)p; }
};

template<class element, class T1> struct isInputForwardIterators
{
	static void constraints(T1 a) {++a; (void)a++; }
	isInputForwardIterators() { void(*p)(T1) = constraints; (void)p; }
};
*/

namespace INMOST
{
	//add more templates as needed
	template<class T> struct make_unsigned;
	template<> struct make_unsigned<char> {typedef unsigned char type;};
	template<> struct make_unsigned<short> {typedef unsigned short type;};
	template<> struct make_unsigned<int> {typedef unsigned int type;};
	template<> struct make_unsigned<long> {typedef unsigned long type;};
	template<> struct make_unsigned<long long> {typedef unsigned long long type;};
	template<> struct make_unsigned<unsigned char> {typedef unsigned char type;};
	template<> struct make_unsigned<unsigned short> {typedef unsigned short type;};
	template<> struct make_unsigned<unsigned int> {typedef unsigned int type;};
	template<> struct make_unsigned<unsigned long> {typedef unsigned long type;};
	template<> struct make_unsigned<unsigned long long> {typedef unsigned long long type;};


#if defined(PACK_ARRAY)
#pragma pack(push,r1,4)
#endif
	// notice: array class would not properly support classes that contain self-references
	//         like std::map
	// notice: next class shell have to implement same algorithms as array
	template<typename element>//, typename enumerator = unsigned int>
	class array
	{
	public:
		typedef unsigned size_type;
		typedef make_unsigned<size_type>::type uenum;
		template<typename etype>
		class _iterator
		{
		private:
			etype * e;
		public:
			typedef etype * pointer;
			typedef etype & reference;
			typedef etype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_iterator():e(NULL){}
			_iterator(etype * i):e(i){}
			_iterator(const _iterator & other){e = other.e;}
			~_iterator() {};
			_iterator operator -(size_t n) const { return _iterator(e-n); }
			_iterator & operator -=(size_t n) { e-=n; return *this; }
			_iterator operator +(size_t n) const { return _iterator(e+n); }
			_iterator & operator +=(size_t n) { e+=n; return *this; }
			_iterator & operator ++(){ ++e; return *this;}
			_iterator operator ++(int){ return _iterator(e++); }
			_iterator & operator --(){ --e; return *this; }
			_iterator operator --(int){ return _iterator(e--); }
			ptrdiff_t operator -(const _iterator & other) const {return e-other.e;}
			etype & operator *() { return *e; }
			const etype & operator *() const { return *e; }
			etype * operator ->() { return e; }
			_iterator & operator =(_iterator const & other) { e = other.e; return *this; }
			bool operator ==(const _iterator & other) const { return e == other.e;}
			bool operator !=(const _iterator & other) const { return e != other.e;}
			bool operator <(const _iterator & other) const { return e < other.e;}
			bool operator >(const _iterator & other) const { return e > other.e;}
			bool operator <=(const _iterator & other) const { return e <= other.e;}
			bool operator >=(const _iterator & other) const { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
			operator const void *() const {return const_cast<const void *> (e);}
		};
		typedef _iterator<element> iterator;
		typedef _iterator<const element> const_iterator;
		template<typename etype>
		class _reverse_iterator
		{
		private:
			etype * e;
		public:
			typedef etype * pointer;
			typedef etype & reference;
			typedef etype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_reverse_iterator():e(NULL){}
			_reverse_iterator(etype * i):e(i){}
			_reverse_iterator(const _reverse_iterator & other){e = other.e;}
			~_reverse_iterator() {};
			_reverse_iterator operator -(size_t n) const { return _reverse_iterator(e+n); }
			_reverse_iterator & operator -=(size_t n) { e+=n; return *this; }
			_reverse_iterator operator +(size_t n) const {return _reverse_iterator(e-n); }
			_reverse_iterator & operator +=(size_t n) { e-=n; return *this; }
			_reverse_iterator & operator ++(){ --e; return *this;}
			_reverse_iterator operator ++(int){ return _reverse_iterator(e--); }
			_reverse_iterator & operator --(){ ++e; return *this; }
			_reverse_iterator operator --(int){ return _reverse_iterator(e++); }
			ptrdiff_t operator -(const _reverse_iterator & other) const {return other.e-e;}
			etype & operator *() { return *e; }
			const etype & operator *() const { return *e; }
			etype * operator ->() { return e; }
			_reverse_iterator & operator =(_reverse_iterator const & other) { e = other.e; return *this;}
			bool operator ==(const _reverse_iterator & other) const { return e == other.e;}
			bool operator !=(const _reverse_iterator & other) const { return e != other.e;}
			bool operator <(const _reverse_iterator & other) const { return e < other.e;}
			bool operator >(const _reverse_iterator & other) const { return e > other.e;}
			bool operator <=(const _reverse_iterator & other) const { return e <= other.e;}
			bool operator >=(const _reverse_iterator & other) const { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
			operator const void *() const {return static_cast<const void *> (e);}
		};
		typedef _reverse_iterator<element> reverse_iterator;
		typedef _reverse_iterator<const element> const_reverse_iterator;
	private:

		element * m_arr;
		size_type m_size;

		//notice: push_back, pop_back, insert, erase are optimized for current growth_formula
		__INLINE static size_type growth_formula(size_type future_size)
		{
			uenum v = static_cast<uenum>(future_size);
			v--;
			v|= (v>>1);
			v|= (v>>2);
			v|= (v>>4);
			v|= (v>>8);
			v|= (v>>16);
			//TODO produces compiler warning
			//if( sizeof(size_type) > 4 ) v|= (v>>32); //must be compile time brench
			v++;
			return static_cast<size_type>(v);
		}
	protected:
		static element * realloc(element * ptr, size_type old_size, size_type new_size)
		{
			size_type gf = growth_formula(new_size);
			if( gf != growth_formula(old_size) )
			{
				element * tmp = new element[gf];
#if defined(DEBUGMEM)
				if( tmp == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				std::copy(ptr,ptr+std::min(old_size,new_size),tmp);
				delete [] ptr;
				ptr = tmp;
				assert(ptr != NULL);
			}
			return ptr;
		}
		static element * alloc(size_type init_size)
		{
			element * tmp = new element[growth_formula(init_size)];
#if defined(DEBUGMEM)
			if( tmp == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
			assert(tmp != NULL);
			return tmp;
		}
		static element * copy(const element * ibeg, const element * iend, element * obeg) 
		{
			if (std::less<const element *>()(ibeg, obeg)) 
			{
				obeg += (iend - ibeg);
				std::copy_backward(ibeg, iend, obeg);
				return obeg;
			} else return std::copy(ibeg, iend, obeg);
		}
	public:
		__INLINE element * data() {return m_arr;}
		__INLINE const element * data() const {return m_arr;}
		array() {m_arr = NULL; m_size = 0;}
		array(size_type n,element c = element())
		{
			m_size = n;
			m_arr = alloc(n);
			std::fill(m_arr,m_arr+m_size,c);
		}
		template<class InputIterator>
		array(InputIterator first, InputIterator last)
		{
			m_size = static_cast<size_type>(std::distance(first,last));
			m_arr = alloc(m_size);
			std::copy(first,last,m_arr);
		}
		array(const array & other)
		{
			m_size = other.m_size;
			if( m_size )
			{
				m_arr = alloc(m_size);
				std::copy(other.m_arr,other.m_arr+other.m_size,m_arr);
			}
			else m_arr = NULL;
		}
		~array() {clear();}
		__INLINE const element & operator [] (size_type n) const
		{
			assert(n < m_size);
			return m_arr[n];
		}
		__INLINE element & operator [] (size_type n)
		{
			assert(n < m_size);
			return m_arr[n];
		}
		__INLINE const element & at (size_type n) const
		{
			assert(n < m_size);
			return m_arr[n];
		}
		__INLINE element & at (size_type n)
		{
			assert(n < m_size);
			return m_arr[n];
		}
		__INLINE element & at_safe (size_type n)
		{
			if( n >= m_size ) resize(n+1);
			return m_arr[n];
		}
		array & operator =(array const & other)
		{
			if( this != &other )
			{
				if(m_arr != NULL )
					clear();
				if( other.m_arr != NULL )
				{
					m_size = other.m_size;
					m_arr = alloc(m_size);
					std::copy(other.m_arr,other.m_arr+other.m_size,m_arr);
				}
			}
			return *this;
		}
		void push_back(const element & e)
		{
			m_arr = realloc(m_arr,m_size,m_size+1);
			m_arr[m_size++] = e;
		}
		void pop_back()
		{
			assert(m_arr != NULL);
			if( m_size > 0)
			{
				m_arr = realloc(m_arr,m_size,m_size-1);
				m_size--;
			}
			else clear();
		}
		__INLINE element & back()
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[m_size-1];
		}
		__INLINE const element & back() const
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[m_size-1];
		}
		__INLINE element & front()
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[0];
		}
		__INLINE const element & front() const
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[0];
		}
		__INLINE size_type capacity() { return growth_formula(m_size); }
		__INLINE bool empty() const { if( m_size ) return false; return true; }
		void resize(size_type n, element c = element() )
		{
			size_type oldsize = m_size;
			m_size = n;
			if( m_size > 0 )
			{
				m_arr = realloc(m_arr,oldsize,m_size);
				if( oldsize < m_size ) std::fill(m_arr+oldsize,m_arr+m_size,c);
			}
			else clear();
		}
		__INLINE size_type size() const {return m_size;}
		__INLINE size_type capacity() const {return growth_formula(m_size);}
		void clear()
		{
			if( m_arr ) delete [] m_arr;
			m_arr = NULL;
			m_size = 0;
		}
		void swap(array<element> & other)
		{
			element * t_m_arr = m_arr;
			size_type t_m_size = m_size;
			m_arr = other.m_arr;
			m_size = other.m_size;
			other.m_arr = t_m_arr;
			other.m_size = t_m_size;
		}
		__INLINE iterator begin() { return m_arr; }
		__INLINE iterator end() { return m_arr+m_size; }
		__INLINE const_iterator begin() const { return m_arr; }
		__INLINE const_iterator end() const { return m_arr+m_size; }
		__INLINE reverse_iterator rbegin() { return reverse_iterator(m_arr+m_size-1); }
		__INLINE reverse_iterator rend() { return reverse_iterator(m_arr-1); }
		__INLINE const_reverse_iterator rbegin() const { return const_reverse_iterator(m_arr+m_size-1); }
		__INLINE const_reverse_iterator rend() const { return const_reverse_iterator(m_arr-1); }
		iterator erase(iterator pos)
		{
			ptrdiff_t d = pos-begin();
			ptrdiff_t s = iterator(m_arr+(m_size-1))-pos;
			m_size--;
			if( m_size > 0 )
			{
				if( s ) copy(m_arr+d+1, m_arr+d+1+s, m_arr+d);
				m_arr = realloc(m_arr, m_size+1, m_size);
			}
			else clear();
			return m_arr+d;
		}
		iterator erase(iterator b, iterator e)
		{
			ptrdiff_t d = b-begin();
			ptrdiff_t s = end()-e;
			ptrdiff_t n = e-b;
			m_size -= n;
			if( m_size > 0 )
			{
				if( s ) copy(m_arr+d+1,m_arr+d+1+s,m_arr+d);
				m_arr = realloc(m_arr,m_size+n,m_size);
				
			}
			else clear();
			return m_arr+d;
		}
		iterator insert(iterator pos, const element & x)
		{
			if( static_cast<void *>(pos) == NULL)
			{
				assert(m_arr == NULL);
				m_size = 1;
				m_arr = alloc(m_size);
				m_arr[0] = x;
				return m_arr;
			}
			else
			{
				ptrdiff_t d = pos-begin();
				ptrdiff_t s = end()-pos;
				m_arr = realloc(m_arr,m_size,m_size+1);
				if( s ) copy(m_arr+d,m_arr+d+s,m_arr+d+1);
				m_arr[d] = x;
				return m_arr+d;
			}
		}
		void insert(iterator pos, size_type n, const element & x)
		{
			if( n > 0 )
			{
				if( static_cast<void *>(pos) == NULL)
				{
					assert(m_arr == NULL);
					m_arr = alloc(n);
					m_size = n;
					std::fill(m_arr,m_arr+m_size,x);
				}
				else
				{
					ptrdiff_t d = pos-begin();
					ptrdiff_t s = end()-pos;
					m_arr = realloc(m_arr,m_size,m_size+n);
					m_size += n;
					if( s ) copy(m_arr+d,m_arr+d+s,m_arr+d+n);
					std::fill(m_arr+d,m_arr+d+n,x);
				}
			}
		}
		template <class InputIterator>
		void insert(iterator pos, InputIterator first, InputIterator last)
		{
			size_type n = static_cast<size_type>(std::distance(first,last));
			if( n > 0 )
			{
				if( static_cast<void *>(pos) == NULL)
				{
					assert(m_arr == NULL);
					m_arr = alloc(n);
					m_size = n;
					std::copy(first,last,m_arr);
				}
				else
				{
					ptrdiff_t d = pos-begin();
					ptrdiff_t s = end()-pos;
					m_arr = realloc(m_arr,m_size,m_size+n);
					m_size += n;
					if( s ) copy(m_arr+d,m_arr+d+s,m_arr+d+n);
					std::copy(first,last,m_arr+d);
				}
			}
		}
		template <class InputIterator>
		void replace(iterator m_first, iterator m_last, InputIterator first, InputIterator last)
		{
			assert( m_size >= 0 );
			size_type n = static_cast<size_type>(std::distance(first,last));
			if( static_cast<void *>(m_first) == NULL && m_arr == NULL)
			{
				assert(m_arr == NULL);
				m_arr = alloc(n);
				m_size = n;
				std::copy(first,last,m_arr);
			}
			else
			{
				ptrdiff_t q = m_last-m_first;
				ptrdiff_t d = m_first-iterator(m_arr);
				ptrdiff_t s = iterator(m_arr+m_size)-m_last;
				if( n-q != 0 )
				{
					m_arr = realloc(m_arr,m_size,m_size+static_cast<size_type>(n-q));	
					m_size+=static_cast<size_type>(n-q);
					if( s ) copy(m_arr+d+q,m_arr+d+q+s,m_arr+d+n);
				}
				std::copy(first,last,m_arr+d);
			}
		}
		template <class InputIterator>
		void assign(InputIterator first, InputIterator last)
		{
			replace(begin(),end(),first,last);
		}
		template<class> friend class shell;
	};
#if defined(PACK_ARRAY)
#pragma pack(pop,r1)
#endif
	template<typename element>
	class shell
	{
	public:
		typedef typename array<element>::size_type size_type;
		template<typename dtype>
		class _iterator
		{
		private:
			dtype * e;
		public:
			typedef dtype * pointer;
			typedef dtype & reference;
			typedef dtype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_iterator():e(NULL){}
			_iterator(dtype * i):e(i){}
			_iterator(const _iterator & other){e = other.e;}
			~_iterator() {};
			_iterator operator -(size_t n) const { return _iterator(e-n); }
			_iterator & operator -=(size_t n) { e-=n; return *this; }
			_iterator operator +(size_t n) const { return _iterator(e+n); }
			_iterator & operator +=(size_t n) { e+=n; return *this; }
			_iterator & operator ++(){ ++e; return *this;}
			_iterator operator ++(int){ return _iterator(e++); }
			_iterator & operator --(){ --e; return *this; }
			_iterator operator --(int){ return _iterator(e--); }
			ptrdiff_t operator -(const _iterator & other) const {return e-other.e;}
			dtype & operator *() { return *e; }
			const dtype & operator *() const { return *e; }
			dtype * operator ->() { return e; }
			_iterator & operator =(_iterator const & other) { e = other.e; return *this; }
			bool operator ==(const _iterator & other) { return e == other.e;}
			bool operator !=(const _iterator & other) { return e != other.e;}
			bool operator <(const _iterator & other) { return e < other.e;}
			bool operator >(const _iterator & other) { return e > other.e;}
			bool operator <=(const _iterator & other) { return e <= other.e;}
			bool operator >=(const _iterator & other) { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
		};
		typedef _iterator<element> iterator;
		typedef _iterator<const element> const_iterator;
		template<typename dtype>
		class _reverse_iterator
		{
		private:
			dtype * e;
		public:
			typedef dtype * pointer;
			typedef dtype & reference;
			typedef dtype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_reverse_iterator():e(NULL){}
			_reverse_iterator(dtype * i):e(i){}
			_reverse_iterator(const _reverse_iterator & other){e = other.e;}
			~_reverse_iterator() {};
			_reverse_iterator operator -(size_t n) const { return _reverse_iterator(e+n); }
			_reverse_iterator & operator -=(size_t n) { e+=n; return *this; }
			_reverse_iterator operator +(size_t n) const {return _reverse_iterator(e-n); }
			_reverse_iterator & operator +=(size_t n) { e-=n; return *this; }
			_reverse_iterator & operator ++(){ --e; return *this;}
			_reverse_iterator operator ++(int){ return _reverse_iterator(e--); }
			_reverse_iterator & operator --(){ ++e; return *this; }
			_reverse_iterator operator --(int){ return _reverse_iterator(e++); }
			ptrdiff_t operator -(const _reverse_iterator & other) const {return other.e-e;}
			dtype & operator *() { return *e; }
			const dtype & operator *() const { return *e; }
			dtype * operator ->() { return e; }
			_reverse_iterator & operator =(_reverse_iterator const & other) { e = other.e; return *this;}
			bool operator ==(const _reverse_iterator & other) { return e == other.e;}
			bool operator !=(const _reverse_iterator & other) { return e != other.e;}
			bool operator <(const _reverse_iterator & other) { return e < other.e;}
			bool operator >(const _reverse_iterator & other) { return e > other.e;}
			bool operator <=(const _reverse_iterator & other) { return e <= other.e;}
			bool operator >=(const _reverse_iterator & other) { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
		};
		typedef _reverse_iterator<element> reverse_iterator;
		typedef _reverse_iterator<const element> const_reverse_iterator;
	private:
		element ** m_arr;
		size_type * m_size;
		element * local_link;
		size_type local_size;
		bool fixed;
	public:
		__INLINE element * data() {return *m_arr;}
		__INLINE const element * data() const {return *m_arr;}
		shell() :m_arr(NULL), m_size(NULL), local_link(NULL), local_size(0), fixed(false) { }
		shell(array<element> & arr) //dynamic
			:m_arr(&arr.m_arr), m_size(&arr.m_size), local_link(NULL), local_size(0), fixed(false) {}
		shell(element * link, size_type size) //fixed
			:m_arr(&local_link), m_size(&local_size), local_link(NULL),local_size(0), fixed(true)
		{
			*m_arr = link;
			*m_size = size;
		}
		shell(const shell & other)
			:m_arr(&local_link), m_size(&local_size), local_link(NULL), local_size(0),fixed(other.fixed)
		{
			if( fixed )
			{
				*m_size = *other.m_size;
				*m_arr = *other.m_arr;
			}
			else
			{
				m_size = other.m_size;
				m_arr = other.m_arr;
			}
		}
		~shell()
		{
		}
		__INLINE const element & operator [] (size_type n) const
		{
			assert(n < *m_size );
			return *((*m_arr)+n);
		}
		__INLINE element & operator [] (size_type n)
		{
			assert(n < *m_size );
			return *((*m_arr)+n);
		}
		__INLINE const element & at (size_type n) const
		{
			assert(n < *m_size );
			return *((*m_arr)+n);
		}
		__INLINE element & at (size_type n)
		{
			assert(n < *m_size );
			return *((*m_arr)+n);
		}
		shell & operator =(shell const & other)
		{
			fixed = other.fixed;
			if( fixed )
			{
				m_size = &local_size;
				m_arr = &local_link;
				*m_size = *other.m_size;
				*m_arr = *other.m_arr;
			}
			else
			{
				m_size = other.m_size;
				m_arr = other.m_arr;
			}
			return *this;
		}
		void push_back(const element & e)
		{
			assert( !fixed ); // array size is fixed
			*m_arr = array<element>::realloc(*m_arr,*m_size,(*m_size)+1);
			(*m_arr)[(*m_size)++] = e;
		}
		void pop_back()
		{
			assert( !fixed ); // array size is fixed
			assert((*m_arr) != NULL);
			if( *m_size > 0)
			{
				*m_arr = array<element>::realloc(*m_arr,*m_size,(*m_size)-1);
				(*m_size)--;
			}
			else clear();
		}
		__INLINE element & back()
		{
			assert(*m_arr != NULL);
			assert(*m_size > 0 );
			return (*m_arr)[(*m_size)-1];
		}
		__INLINE const element & back() const
		{
			assert(*m_arr != NULL);
			assert(*m_size > 0 );
			return (*m_arr)[(*m_size)-1];
		}
		__INLINE element & front()
		{
			assert(*m_arr != NULL);
			assert(*m_size > 0 );
			return (*m_arr)[0];
		}
		__INLINE const element & front() const
		{
			assert(*m_arr != NULL);
			assert(*m_size > 0 );
			return (*m_arr)[0];
		}
		__INLINE size_type capacity() { return array<element>::growth_formula(*m_size); }
		__INLINE bool empty() const { if( *m_size ) return false; return true; }
		void resize(size_type n, element c = element() )
		{
			assert( !fixed ); // array size is fixed
			size_type oldsize = *m_size;
			*m_size = n;
			if( *m_size > 0 )
			{
				*m_arr = array<element>::realloc(*m_arr,oldsize,*m_size);
				if( oldsize < *m_size ) std::fill((*m_arr)+oldsize,(*m_arr)+(*m_size),c);
			}
			else clear();
		}
		__INLINE size_type size() const {return *m_size;}
		void clear()
		{
			assert( !fixed ); // array size is fixed
			if( *m_arr ) delete [] *m_arr;
			*m_arr = NULL;
			*m_size = 0;
		}
		void swap(shell<element> & other)
		{
			element * t_m_arr = *m_arr;
			size_type t_m_size = *m_size;
			*m_arr = *other.m_arr;
			*m_size = *other.m_size;
			*other.m_arr = t_m_arr;
			*other.m_size = t_m_size;
		}
		__INLINE iterator begin() { return *m_arr; }
		__INLINE iterator end() { return (*m_arr)+(*m_size); }
		__INLINE const_iterator begin() const { return (*m_arr); }
		__INLINE const_iterator end() const { return (*m_arr)+(*m_size); }
		__INLINE reverse_iterator rbegin() { return reverse_iterator((*m_arr)+(*m_size)-1); }
		__INLINE reverse_iterator rend() { return reverse_iterator((*m_arr)-1); }
		__INLINE const_reverse_iterator rbegin() const { return const_reverse_iterator((*m_arr)+(*m_size)-1); }
		__INLINE const_reverse_iterator rend() const { return const_reverse_iterator((*m_arr)-1); }
		iterator erase(iterator pos)
		{
			assert( !fixed ); // array size is fixed
			ptrdiff_t d = pos-begin();
			ptrdiff_t s = iterator((*m_arr)+(*m_size)-1)-pos;
			(*m_size)--;
			if( *m_size > 0 )
			{
				if( s ) array<element>::copy((*m_arr)+d+1, (*m_arr)+d+1+s, (*m_arr)+d);
				*m_arr = array<element>::realloc(*m_arr, (*m_size)+1, *m_size);
			}
			else clear();
			return (*m_arr)+d;
		}
		iterator erase(iterator b, iterator e)
		{
			assert( !fixed ); // array size is fixed
			ptrdiff_t d = b-begin();
			ptrdiff_t s = end()-e;
			ptrdiff_t n = e-b;
			(*m_size) -= n;
			if( *m_size > 0 )
			{
				if( s ) array<element>::copy((*m_arr)+d+1,(*m_arr)+d+1+s,(*m_arr)+d);
				*m_arr = array<element>::realloc(*m_arr,(*m_size)+n,*m_size);
				
			}
			else clear();
			return (*m_arr)+d;
		}
		iterator insert(iterator pos, const element & x)
		{
			assert( !fixed ); // array size is fixed
			if( static_cast<void *>(pos) == NULL )
			{
				assert((*m_arr) == NULL);
				*m_size = 1;
				*m_arr = array<element>::alloc(*m_size);
				(*m_arr)[0] = x;
				return *m_arr;
			}
			else
			{
				ptrdiff_t d = pos-begin();
				ptrdiff_t s = end()-pos;
				*m_arr = array<element>::realloc(*m_arr,*m_size,(*m_size)+1);
				(*m_size)++;
				if( s ) array<element>::copy((*m_arr)+d,(*m_arr)+d+s,(*m_arr)+d+1);
				(*m_arr)[d] = x;
				return (*m_arr)+d;
			}
		}
		void insert(iterator pos, size_type n, const element & x)
		{
			assert( !fixed ); // array size is fixed
			if( n )
			{
				if( static_cast<void *>(pos) == NULL)
				{
					assert((*m_arr) == NULL);
					*m_arr = array<element>::alloc(n);
					*m_size = n;
					std::fill(*m_arr,(*m_arr)+(*m_size),x);
				}
				else
				{
					ptrdiff_t d = pos-iterator(*m_arr);
					ptrdiff_t s = iterator((*m_arr)+(*m_size))-pos;
					*m_arr = array<element>::realloc(*m_arr,*m_size,(*m_size)+n);
					*m_size += n;
					if( s ) array<element>::copy((*m_arr)+d,(*m_arr)+d+s,(*m_arr)+d+n);
					std::fill((*m_arr)+d,(*m_arr)+d+n,x);
				}
			}
		}
		template <class InputIterator>
		void insert(iterator pos, InputIterator first, InputIterator last)
		{
			assert( !fixed ); // array size is fixed
			size_type n = static_cast<size_type>(std::distance(first,last));
			if( n )
			{
				if( static_cast<void *>(pos) == NULL)
				{
					assert((*m_arr) == NULL);
					*m_arr = array<element>::alloc(n);
					*m_size = n;
					std::copy(first,last,*m_arr);
				}
				else
				{
					ptrdiff_t d = pos-iterator(*m_arr);
					ptrdiff_t s = iterator((*m_arr)+(*m_size))-pos;
					*m_arr = array<element>::realloc(*m_arr,*m_size,(*m_size)+n);
					*m_size += n;
					if( s ) array<element>::copy((*m_arr)+d,(*m_arr)+d+s,(*m_arr)+d+n);
					std::copy(first,last,(*m_arr)+d);
				}
			}
		}
		template <class InputIterator>
		void replace(iterator m_first, iterator m_last, InputIterator first, InputIterator last)
		{
			size_type n = static_cast<size_type>(std::distance(first,last));
			if( static_cast<void *>(m_first) == NULL)
			{
				assert((*m_arr)==NULL);
				*m_arr = array<element>::alloc(n);
				*m_size = n;
				std::copy(first,last,*m_arr);
			}
			else
			{
				ptrdiff_t q = m_last-m_first;
				ptrdiff_t d = m_first-iterator(*m_arr);
				ptrdiff_t s = iterator((*m_arr)+(*m_size))-m_last;
				if( n-q != 0 )
				{
					*m_arr = array<element>::realloc(*m_arr,*m_size,(*m_size)+static_cast<size_type>(n-q));	
					(*m_size)+=static_cast<size_type>(n-q);
					if( s ) array<element>::copy((*m_arr)+d+q,(*m_arr)+d+q+s,(*m_arr)+d+n);
				}
				std::copy(first,last,(*m_arr)+d);
			}
		}
		template <class InputIterator>
		void assign(InputIterator first, InputIterator last)
		{
			replace(begin(),end(),first,last);
		}
	};

	template<typename IndType,typename ValType>
	class interval
	{
	public:
		typedef ValType * iterator;
		typedef ValType const * const_iterator;
	private:
		ValType * parray;
		IndType beg_index, end_index;
	public:
		void clear()
		{
			if (!empty()) delete [] (parray + beg_index);
			parray = NULL;
			end_index = beg_index;
		}
		void swap(interval<IndType, ValType> & other)
		{
			{
				IndType tmp = beg_index;
				beg_index = other.beg_index;
				other.beg_index = tmp;
				tmp = end_index;
				end_index = other.end_index;
				other.end_index = tmp;
			}
			{
				ValType * tmp = parray;
				parray = other.parray;
				other.parray = tmp;
			}
		}
		interval()
		{
			beg_index = 0;
			end_index = 0;
			parray = NULL;
		}
		interval(IndType beg)
		{
			beg_index = beg;
			end_index = beg_index;
			parray = NULL;
		}
		interval(IndType beg, IndType end, ValType c = ValType())
		{
			beg_index = beg;
			end_index = end;
			if (beg != end)
			{
				ValType * tmp = new ValType[end_index-beg_index];
#if defined(DEBUGMEM)
				if( tmp == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				parray = tmp;
				assert(parray != NULL);
				parray = parray - beg_index;
				std::fill(parray+beg_index,parray+end_index,c);
			}
			else parray = NULL;
		}
		interval(const interval & other)
		{
			beg_index = other.beg_index;
			end_index = other.end_index;
			if( beg_index != end_index )
			{
				ValType * tmp = new ValType[end_index-beg_index];
#if defined(DEBUGMEM)
				if( tmp == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				parray = static_cast<ValType *>(tmp);
				assert(parray != NULL);
				parray = parray - beg_index;
				std::copy(other.parray+beg_index,other.parray+end_index,parray+beg_index);
			}
			else parray = NULL;
		}
		~interval()
		{
			clear();
		}
		interval & operator =(interval const & other)
		{
			if( &other != this )
			{
				clear();
				beg_index = other.beg_index;
				end_index = other.end_index;
				if( beg_index != end_index )
				{
					ValType * tmp = new ValType[end_index-beg_index];
#if defined(DEBUGMEM)
					if( tmp == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
					parray = static_cast<ValType *>(tmp);
					assert(parray != NULL);
					parray = parray - beg_index;
					std::copy(other.parray+beg_index,other.parray+end_index,parray+beg_index);
				}
			}
			return *this;
		}
		ValType & at(IndType row)
		{
			assert(row >= beg_index);
			assert(row < end_index);
			return parray[row];
		}
		const ValType & at(IndType row) const
		{
			assert(row >= beg_index);
			assert(row < end_index);
			return parray[row];
		}
		ValType & operator [](IndType row)
		{
			assert(row >= beg_index );
			assert(row < end_index );
			return parray[row];
		}
		const ValType & operator [](IndType row) const
		{
			assert(row >= beg_index );
			assert(row < end_index );
			return parray[row];
		}
		void set_interval_beg(IndType beg)
		{
			IndType shift = beg-beg_index;
			shift_interval(shift);
		}
		void set_interval_end(IndType end, const ValType & c = ValType())
		{
			if( end == end_index ) return;
			if( beg_index != end )
			{
				ValType * parray_new = new ValType[end-beg_index];
#if defined(DEBUGMEM)
				if( parray_new == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				assert(parray_new != NULL);
				parray_new = parray_new - beg_index;
				//IndType mbeg = beg_index;
				IndType mend = std::min(end,end_index);
				if( !empty() ) 
				{
					//~ std::cout << beg_index << ":" << end_index << " new end " << end << " old ptr " << (void*)parray << " new ptr " << (void *)parray_new << std::endl;
					std::copy(parray+beg_index,parray+mend,parray_new+beg_index);
					clear();
				}
				if( mend < end ) std::fill(parray_new+mend,parray_new+end,c);
				parray = parray_new;
			}
			else clear();
			end_index = end;
		}

		void shift_interval(IndType shift)
		{
			parray = parray + beg_index;
			beg_index += shift;
			end_index += shift;
			parray = parray - beg_index;
		}
		iterator begin() {return parray+beg_index;}
		const_iterator begin() const {return parray+beg_index;}
		iterator end() {return parray + end_index;}
		const_iterator end() const {return parray + end_index;}
		IndType get_interval_beg() const { return beg_index; }
		IndType get_interval_end() const { return end_index; }
		int size() const {return end_index - beg_index;}
		bool empty() const {return beg_index == end_index;}
	};

	template<typename element, unsigned int stacked>
	class dynarray
	{
	public:
		typedef size_t size_type;
		template<typename dtype>
		class _iterator
		{
		private:
			dtype * e;
		public:
			typedef dtype * pointer;
			typedef dtype & reference;
			typedef dtype value_type;
			typedef size_type difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_iterator():e(NULL){}
			_iterator(dtype * i):e(i){}
			_iterator(const _iterator & other){e = other.e;}
			~_iterator() {};
			_iterator operator -(size_type n) const { return _iterator(e-n); }
			_iterator & operator -=(size_type n) { e-=n; return *this; }
			_iterator operator +(size_type n) const { return _iterator(e+n); }
			_iterator & operator +=(size_type n) { e+=n; return *this; }
			_iterator & operator ++(){ ++e; return *this;}
			_iterator operator ++(int){ return _iterator(e++); }
			_iterator & operator --(){ --e; return *this; }
			_iterator operator --(int){ return _iterator(e--); }
			size_type operator -(const _iterator & other) const {return static_cast<size_type>(e-other.e);}
			dtype & operator *() { return *e; }
			const dtype & operator *() const { return *e; }
			dtype * operator ->() { return e; }
			_iterator & operator =(_iterator const & other) { e = other.e; return *this; }
			bool operator ==(const _iterator & other) const { return e == other.e;}
			bool operator !=(const _iterator & other) const { return e != other.e;}
			bool operator <(const _iterator & other) const { return e < other.e;}
			bool operator >(const _iterator & other) const { return e > other.e;}
			bool operator <=(const _iterator & other) const { return e <= other.e;}
			bool operator >=(const _iterator & other) const { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
			operator const void *() const {return static_cast<const void *> (e);}
		};
		typedef _iterator<element> iterator;
		typedef _iterator<const element> const_iterator;
		template<typename dtype>
		class _reverse_iterator
		{
		private:
			dtype * e;
		public:
			typedef dtype * pointer;
			typedef dtype & reference;
			typedef dtype value_type;
			typedef size_type difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_reverse_iterator():e(NULL){}
			_reverse_iterator(dtype * i):e(i){}
			_reverse_iterator(const _reverse_iterator & other){e = other.e;}
			~_reverse_iterator() {};
			_reverse_iterator operator -(size_type n) const { return _reverse_iterator(e+n); }
			_reverse_iterator & operator -=(size_type n) { e+=n; return *this; }
			_reverse_iterator operator +(size_type n) const {return _reverse_iterator(e-n); }
			_reverse_iterator & operator +=(size_type n) { e-=n; return *this; }
			_reverse_iterator & operator ++(){ --e; return *this;}
			_reverse_iterator operator ++(int){ return _reverse_iterator(e--); }
			_reverse_iterator & operator --(){ ++e; return *this; }
			_reverse_iterator operator --(int){ return _reverse_iterator(e++); }
			size_type operator -(const _reverse_iterator & other) const {return static_cast<size_type>(other.e-e);}
			dtype & operator *() { return *e; }
			const dtype & operator *() const { return *e; }
			dtype * operator ->() { return e; }
			_reverse_iterator & operator =(_reverse_iterator const & other) { e = other.e; return *this;}
			bool operator ==(const _reverse_iterator & other) const { return e == other.e;}
			bool operator !=(const _reverse_iterator & other) const { return e != other.e;}
			bool operator <(const _reverse_iterator & other) const { return e < other.e;}
			bool operator >(const _reverse_iterator & other) const { return e > other.e;}
			bool operator <=(const _reverse_iterator & other) const { return e <= other.e;}
			bool operator >=(const _reverse_iterator & other) const { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
			operator const void *() const {return static_cast<const void *> (e);}
		};
		typedef _reverse_iterator<element> reverse_iterator;
		typedef _reverse_iterator<const element> const_reverse_iterator;
	private:
		element stack[stacked];
		element * pbegin;
		element * pend;
		element * preserved;
		void preallocate(size_type n)
		{
			if( n <= static_cast<size_type>(stacked) )
			{
				pbegin = stack;
				pend = pbegin + n;
				preserved = stack+static_cast<size_type>(stacked);
			}
			else
			{
				pbegin = new element[n];
#if defined(DEBUGMEM)
				if( pbegin == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				assert(pbegin != NULL);
				pend = pbegin+n;
				preserved = pbegin+n;
			}
		}
		static element * copy(const element * ibeg, const element * iend, element * obeg) 
		{
			if (std::less<const element *>()(ibeg, obeg)) 
			{
				obeg += (iend - ibeg);
				std::copy_backward(ibeg, iend, obeg);
				return obeg;
			} else return std::copy(ibeg, iend, obeg);
		}
	public:
		__INLINE element * data() {return pbegin;}
		__INLINE const element * data() const {return pbegin;}
		void report_addr()
		{
			std::cout << "stack:     " << &stack << std::endl;
			std::cout << "pbegin:    " << pbegin << std::endl;
			std::cout << "pend:      " << pend << std::endl;
			std::cout << "preserved: " << preserved << std::endl;
			std::cout << "size:      " << pend-pbegin << std::endl;
			std::cout << "reserved:  " << preserved-pbegin << std::endl;
		}
		void reserve(size_type n)
		{
			//std::cout << n << std::endl;
			size_type k = size();
			if( n > static_cast<size_type>(stacked) )
			{
				if( pbegin == stack )
				{
					pbegin = new element[n];
#if defined(DEBUGMEM)
					if( pbegin == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
					assert(pbegin != NULL);
					std::copy(stack,stack+k,pbegin);
				}
				else //if( n > capacity() )
				{
					element * pbegin_new = new element[n];
#if defined(DEBUGMEM)
					if( pbegin_new == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
					assert(pbegin_new != NULL);
					std::copy(pbegin,pbegin+k,pbegin_new);
					delete [] pbegin;
					pbegin = pbegin_new;
				}
				pend = pbegin+ (k < n ? k : n);
				preserved = pbegin + n;
			}
		}
		dynarray() : stack()
		{
			pbegin = pend = stack;
			preserved = stack+static_cast<size_type>(stacked);
		}
		dynarray(size_type n,element c = element()) : stack()
		{
			preallocate(n);
			std::fill(pbegin,pend,c);
		}
		template<class InputIterator>
		dynarray(InputIterator first, InputIterator last) : stack()
		{
			preallocate(static_cast<size_type>(std::distance(first,last)));
			std::copy(first,last,pbegin);
		}
		dynarray(const dynarray & other)
		{
			preallocate(other.size());
			std::copy(other.pbegin,other.pend,pbegin);
		}
		~dynarray()
		{
			if( pbegin != stack ) delete [] pbegin;
		}
		dynarray & operator =(dynarray const & other)
		{
			if( this != &other )
			{
				if(pbegin != stack) delete [] pbegin;
				preallocate(other.size());
				std::copy(other.pbegin,other.pend,pbegin);
			}
			return *this;
		}
		__INLINE const element & operator [] (size_type n) const
		{
			assert(pbegin+n < pend);
			return pbegin[n];
		}
		__INLINE element & operator [] (size_type n)
		{
			assert(pbegin+n < pend);
			return pbegin[n];
		}
		__INLINE const element & at (size_type n) const
		{
			assert(pbegin+n < pend);
			return pbegin[n];
		}
		__INLINE element & at (size_type n)
		{
			assert(pbegin+n < pend);
			return pbegin[n];
		}

		void push_back(const element & e)
		{
			if( pend == preserved ) reserve(capacity()*2);
			(*pend) = e;
			pend++;
		}
		void pop_back()
		{
			pend--;
		}
		__INLINE element & back() {return *(pend-1);}
		__INLINE const element & back() const {return *(pend-1);}
		__INLINE element & front() {return *pbegin;}
		__INLINE const element & front() const {return *pbegin;}
		__INLINE size_type capacity() { return static_cast<size_type>(preserved-pbegin); }
		__INLINE bool empty() const { return pend==pbegin; }
		void resize(size_type n, element c = element() )
		{
			size_type oldsize = size();
			size_type capsize0 = capacity();
			size_type capsize = std::max<size_type>(capsize0,1);
			while( capsize < n ) capsize *= 2;
			reserve(capsize);
			if( n > oldsize ) std::fill(pbegin+oldsize,pbegin+n,c);
			pend = pbegin + n;
		}
		__INLINE size_type size() const {return static_cast<size_type>(pend-pbegin);}
		void clear()
		{
			pend = pbegin;
			//pbegin = pend = stack;
			//preserved = stack+static_cast<size_type>(stacked);
		}
		void move(dynarray<element,stacked> & other)
		{
			size_type k = size(), n = other.size();
			if( n > static_cast<size_type>(stacked) ) delete [] other.pbegin;//free(other.pbegin);
			if( k > static_cast<size_type>(stacked) )
			{
				other.pbegin = pbegin;
				other.pend = pend;
				other.preserved = preserved;
			}
			else
			{
				std::copy(stack,stack+k,other.stack);
				other.pbegin = other.stack;
				other.pend = other.stack+k;
				other.preserved = other.stack+static_cast<size_type>(stacked);
			}
			pbegin = pend = stack;
			preserved = stack+static_cast<size_type>(stacked);
		}
		void swap(dynarray<element,stacked> & other)
		{
			size_type k = size(), n = other.size();
			bool konstack = (pbegin == stack);
			bool nonstack = (other.pbegin == other.stack);
			if( konstack && nonstack )
			{
				element temp[stacked];
				std::copy(stack,stack+k,temp);
				std::copy(other.stack,other.stack+n,stack);
				std::copy(temp,temp+k,other.stack);
				other.pend = other.pbegin+k;
				pend = pbegin+n;
			}
			else if( konstack && !nonstack )
			{
				std::copy(stack,stack+k,other.stack);
				pbegin = other.pbegin;
				pend = other.pend;
				preserved = other.preserved;
				other.pbegin = other.stack;
				other.pend = other.stack+k;
				other.preserved = other.stack+static_cast<size_type>(stacked);
			}
			else if( !konstack && nonstack )
			{
				std::copy(other.stack,other.stack+n,stack);
				other.pbegin = pbegin;
				other.pend = pend;
				other.preserved = preserved;
				pbegin = stack;
				pend = stack+n;
				preserved = stack+static_cast<size_type>(stacked);
			}
			else
			{
				element * temp;
				temp = pbegin;
				pbegin = other.pbegin;
				other.pbegin = temp;
				temp = pend;
				pend = other.pend;
				other.pend = temp;
				temp = preserved;
				preserved = other.preserved;
				other.preserved = temp;
			}
		}
		__INLINE iterator begin() { return pbegin; }
		__INLINE iterator end() { return pend; }
		__INLINE const_iterator begin() const { return pbegin; }
		__INLINE const_iterator end() const { return pend; }
		__INLINE reverse_iterator rbegin() { return reverse_iterator(pend-1); }
		__INLINE reverse_iterator rend() { return reverse_iterator(pbegin-1); }
		__INLINE const_reverse_iterator rbegin() const { return const_reverse_iterator(pend-1); }
		__INLINE const_reverse_iterator rend() const { return const_reverse_iterator(pbegin-1); }
		iterator erase(iterator pos)
		{
			ptrdiff_t d = pos-begin();
			ptrdiff_t s = (end()-pos)-1;
			if( s ) copy(pbegin+d+1,pbegin+d+1+s,pbegin+d);
			pend--;
			return pbegin+d;
		}
		iterator erase(iterator b, iterator e)
		{
			ptrdiff_t d = b-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-e;
			ptrdiff_t n = e-b;
			if( s ) copy(pbegin+d+n,pbegin+d+n+s,pbegin+d);
			pend -= n;
			return pbegin+d;
		}
		iterator insert(iterator pos, const element & x)
		{
			ptrdiff_t d = pos-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-pos;
			if( pend == preserved ) reserve(capacity()*2);
			pend++;
			//if( s > 0 ) memmove((void *)(pbegin+d+1),(void *)(pbegin+d),sizeof(element)*s);
			//new (pbegin+d) element(x);
			if( s ) copy(pbegin+d,pbegin+d+s,pbegin+d+1);
			pbegin[d] = s;
			return pbegin+d;
		}
		void insert(iterator pos, size_type n, const element & x)
		{
			ptrdiff_t d = pos-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-pos;
			size_type capsize = std::max<size_type>(capacity(),1);
			while( capsize < n + size() ) capsize *= 2;
			reserve(capsize);
			//if( s > 0 ) memmove((void *)(pbegin+d+n),(void *)(pbegin+d),sizeof(element)*s);
			if( s ) copy(pbegin+d,pbegin+d+s,pbegin+d+n);
			pend+=n;
			//for(size_type i = 0; i < n; i++) new (pbegin+d+i) element(x);
			std::fill(pbegin+d,pbegin+d+n,x);
		}
		template <class InputIterator>
		void insert(iterator pos, InputIterator first, InputIterator last)
		{
			ptrdiff_t n = static_cast<ptrdiff_t>(std::distance(first,last));
			ptrdiff_t d = pos-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-pos;
			size_type capsize = std::max<size_type>(capacity(),1);
			while( capsize < size()+n ) capsize *= 2;
			reserve(capsize);
			if( s ) copy(pbegin+d,pbegin+d+s,pbegin+d+n);
			std::copy(first,last,pbegin+d);
			pend+=n;
		}
	};
	
	//this is supposed to provide thread-safe storage type for chunks in chunk_array
	//the idea is that pointer to block is never moved
	//so that access to any entry before size() is correct on resize
	//critical section is required on resize()
	template<typename element, size_t base = 512>
	class linked_array
	{
		element e[base]; //storage of current node
		size_t ne; //number of elements in array
		linked_array * next; //link to next node
		linked_array(const linked_array & b) {} //don't copy me
		linked_array & operator = (linked_array const & b) { return *this; } //don't copy me
	public:
		linked_array() : ne(0), next(NULL) {}
		void clear()
		{
			if (next)
				next->clear();
			for (size_t k = 0; k < ne; ++k)
				e[k].~element();
			ne = 0;
		}
		element & operator [](size_t n)
		{
			if (n < base)
			{
				assert(n < ne);
				return e[n];
			}
			else
			{
				assert(next);
				return next->operator [](n-base);
			}
		}
		const element & operator [](size_t n) const
		{
			if (n < base)
			{
				assert(n < ne);
				return e[n];
			}
			else
			{
				assert(next);
				return next->operator [](n - base);
			}
		}
		size_t size() const
		{
			if( ne == base && next )
				return base+next()->size();
			else
				return ne;
		}
		void resize(size_t new_n, const element & c =  element())
		{
			if (new_n <= base)
			{
				for (size_t k = new_n; k < ne; ++k)
					e[k].~element();
				for (size_t k = ne; k < new_n; ++k)
					new (&e[k]) element(c);
				ne = new_n;
				if( next ) next->clear();
			}
			else
			{
				if( !next )
					next = new linked_array;
				if (ne < base)
				{
					for (size_t k = ne; k < base; ++k)
						new (&e[k]) element(c);
					ne = base;
				}
				next->resize(new_n-base);
			}
		}
		~linked_array()
		{
			clear();
			if( next ) delete next;
		}
	};
	
	template<typename element, int block_bits>
	class chunk_array
	{
	public:
		typedef size_t size_type; //need signed type for reverse iterator? (reverse iterators commented)
		typedef make_unsigned<size_type>::type uenum; //this is required for right bit shifts not to create leading 1s
	private:
		static size_type const block_bits_mask = (1 << (block_bits)) - 1;
		static size_type const block_val = 1 << block_bits;
		static size_type const block_size = sizeof(element)*block_val;
		static size_type const fwd_alloc_chunk_bits = 6;
		static size_type const fwd_alloc_chunk_bits_mask = (1 << (fwd_alloc_chunk_bits))-1;
		static size_type const fwd_alloc_chunk_val = 1 << fwd_alloc_chunk_bits;
		static size_type const fwd_alloc_chunk_size = sizeof(element *)*fwd_alloc_chunk_val;

		linked_array<element *> chunks;


		size_type m_size;
		//This neads static_cast to unsigned to
		__INLINE size_type GetChunkNumber(size_type k) const {return static_cast<uenum>(k) >> block_bits;}
		__INLINE size_type GetElementNumber(size_type k) const {return (k & block_bits_mask);}
		__INLINE element * access_block(size_type k) {return chunks[GetChunkNumber(k)];}
		__INLINE const element * access_block(size_type k) const {return chunks[GetChunkNumber(k)];}
		__INLINE element & access_element(size_type k) {return access_block(k)[GetElementNumber(k)];}
		__INLINE const element & access_element(size_type k) const {return access_block(k)[GetElementNumber(k)];}
	public:
		__INLINE element & operator [] (size_type i)
		{
			assert(i < size());
			return access_element(i);
		}
		__INLINE const element & operator [] (size_type i) const
		{
			assert(i < size());
			return access_element(i);
		}
		__INLINE element & at(size_type i)
		{
			assert(i < size());
			return access_element(i);
		}
		__INLINE const element & at(size_type i) const
		{
			assert(i < size());
			return access_element(i);
		}




		class iterator
		{
		private:
			chunk_array<element, block_bits> * link;
			size_type pos;
		public:
			typedef element * pointer;
			typedef element & reference;
			typedef element value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			iterator(){pos = 0; link = NULL;}
			iterator(chunk_array<element,block_bits> * c, size_type pos):link(c),pos(pos){}
			iterator(const iterator & other){link = other.link; pos = other.pos;}
			~iterator() {};
			iterator operator -(size_type n) const { return iterator(link,pos-n); }
			iterator & operator -=(size_type n) { pos-=n; return *this; }
			iterator operator +(size_type n) const { return iterator(link,pos+n); }
			iterator & operator +=(size_type n) { pos+=n; return *this; }
			iterator & operator ++(){ ++pos; return *this;}
			iterator operator ++(int){ return iterator(link,pos++); }
			iterator & operator --(){ --pos; return *this; }
			iterator operator --(int){ return iterator(link,pos--); }
			ptrdiff_t operator -(const iterator & other) const {return pos-other.pos;}
			element & operator *() { return link->at(pos); }
			const element & operator *() const { return link->at(pos); }
			element * operator ->() { return &link->at(pos); }
			iterator & operator =(iterator const & other) { link = other.link; pos = other.pos; return *this; }
			bool operator ==(const iterator & other) const { assert(link == other.link); return pos == other.pos;}
			bool operator !=(const iterator & other) const { assert(link == other.link); return pos != other.pos;}
			bool operator <(const iterator & other) const { assert(link == other.link); return pos < other.pos;}
			bool operator >(const iterator & other) const { assert(link == other.link); return pos > other.pos;}
			bool operator <=(const iterator & other) const { assert(link == other.link); return pos <= other.pos;}
			bool operator >=(const iterator & other) const { assert(link == other.link); return pos >= other.pos;}
		};

		class const_iterator
		{
		private:
			const chunk_array<element,block_bits> * const link;
			size_type pos;
		public:
			typedef element * pointer;
			typedef element & reference;
			typedef element value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			const_iterator(){pos = 0; link = NULL;}
			const_iterator(const chunk_array<element,block_bits> * c, size_type pos):link(c),pos(pos){}
			const_iterator(const const_iterator & other) :link(other.link) {pos = other.pos;}
			~const_iterator() {};
			const_iterator operator -(size_type n) const { return const_iterator(link,pos-n); }
			const_iterator & operator -=(size_type n) { pos-=n; return *this; }
			const_iterator operator +(size_type n) const { return const_iterator(link,pos+n); }
			const_iterator & operator +=(size_type n) { pos+=n; return *this; }
			const_iterator & operator ++(){ ++pos; return *this;}
			const_iterator operator ++(int){ return const_iterator(link,pos++); }
			const_iterator & operator --(){ --pos; return *this; }
			const_iterator operator --(int){ return const_iterator(link,pos--); }
			ptrdiff_t operator -(const const_iterator & other) const {return pos-other.pos;}
			const element & operator *() { return link->at(pos); }
			const element * operator ->() { return &link->at(pos); }
			const_iterator & operator =(const_iterator const & other) { link = other.link; pos = other.pos; return *this; }
			bool operator ==(const const_iterator & other) const { assert(link == other.link); return pos == other.pos;}
			bool operator !=(const const_iterator & other) const { assert(link == other.link); return pos != other.pos;}
			bool operator <(const const_iterator & other) const { assert(link == other.link); return pos < other.pos;}
			bool operator >(const const_iterator & other) const { assert(link == other.link); return pos > other.pos;}
			bool operator <=(const const_iterator & other) const { assert(link == other.link); return pos <= other.pos;}
			bool operator >=(const const_iterator & other) const { assert(link == other.link); return pos >= other.pos;}
		};

		void inner_resize(size_type new_size)
		{
			size_type oldnchunks2 = (static_cast<uenum>(m_size) >> block_bits) + ( (m_size & block_bits_mask) ? 1 : 0);
			size_type oldn = (static_cast<uenum>(oldnchunks2) >> fwd_alloc_chunk_bits) + ( (oldnchunks2 & fwd_alloc_chunk_bits_mask) ? 1 : 0);
			size_type newnchunks2 = (static_cast<uenum>(new_size) >> block_bits) + ( (new_size & block_bits_mask)? 1 : 0);
			size_type newn = (static_cast<uenum>(newnchunks2) >> fwd_alloc_chunk_bits) + ( (newnchunks2 & fwd_alloc_chunk_bits_mask) ? 1 : 0);
			for(size_type q = newnchunks2; q < oldnchunks2; q++)
			{
				assert(chunks[q] != NULL);
				delete [] chunks[q];
				chunks[q] = NULL;
			}
			if( newn != oldn )
			{
				if( newn > 0 )
					chunks.resize(fwd_alloc_chunk_size*newn);
				else
					chunks.clear();
			}
			for(size_type q = oldnchunks2; q < newnchunks2; q++)
			{
				assert(chunks[q] == NULL);
				chunks[q] = new element[block_size];//(element *)malloc(block_size);
#if defined(DEBUGMEM)
				if( chunks[q] == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				assert(chunks[q] != NULL);
			}
			//for(size_type q = m_size; q < new_size; q++) new (&access_element(q)) element();
		}
	public:

		size_type capacity() const
		{
			size_type chunks = ((static_cast<uenum>(m_size)>>block_bits) + ((m_size & block_bits_mask)? 1 : 0));
			return chunks*block_val;
		}
		size_type size() const {return m_size;}
		bool empty() const {return size() == 0;}
		void clear()
		{
			size_type cend = (static_cast<uenum>(m_size) >> block_bits) + ((m_size & block_bits_mask)? 1 : 0);
			for(size_type q = 0; q < cend; q++)
			{
				//free(chunks[q]);
				delete [] chunks[q];
				chunks[q] = NULL;
			}
			chunks.clear();
			m_size = 0;
		}
		chunk_array()
		{
			m_size = 0;
		}
		chunk_array(const chunk_array & other)
		{
			m_size = 0;
			inner_resize(other.size());
			for(size_type k = 0; k < other.size(); k++)
				new(&access_element(k)) element(other.access_element(k));
			m_size = other.size();
		}
		chunk_array & operator =(chunk_array const & other)
		{
			if( this != &other )
			{
				clear();
				inner_resize(other.size());
				for(size_type k = 0; k < other.size(); k++)
					new(&access_element(k)) element(other.access_element(k));
				m_size = other.size();
			}
			return *this;
		}
		~chunk_array()
		{
			clear();
		}

		element & front() { assert(!empty()); return access_element(0);}
		element & back() { assert(!empty()); return access_element(size()-1);}
		const element & front() const { assert(!empty()); return access_element(0);}
		const element & back() const { assert(!empty()); return access_element(size()-1);}
		void pop_back()
		{
			assert(!empty());
			inner_resize(m_size-1);
			m_size--;
		}
		void push_back(const element & e)
		{
			inner_resize(m_size+1);
			new (&access_element(m_size)) element(e);
			m_size++;
		}
		void resize(size_type n, const element & e = element())
		{
			inner_resize(n);
			for(size_type k = m_size; k < n; k++)
				new (&access_element(k)) element(e);
			m_size = n;
		}

		iterator erase(iterator pos)
		{
			iterator it = pos, jt = it++;
			while(it != end()) (*jt++) = (*it++);
			inner_resize(m_size-1);
			m_size--;
			return pos;
		}

		iterator begin() {return iterator(this,0);}
		iterator end() {return iterator(this,m_size);}
		const_iterator begin() const {return const_iterator(this,0);}
		const_iterator end() const {return const_iterator(this,m_size);}
	};


	template<int block_bits>
	class chunk_bulk_array
	{
	public:
		typedef size_t size_type;
		typedef make_unsigned<size_type>::type uenum; //this is required for right shift not to create leading 1s
	private:
		static size_type const block_bits_mask = (1 << (block_bits)) - 1;
		static size_type const block_val = 1 << block_bits;
		static size_type const block_size = sizeof(char)*block_val;
		static size_type const fwd_alloc_chunk_bits = 6;
		static size_type const fwd_alloc_chunk_bits_mask = (1 << (fwd_alloc_chunk_bits))-1;
		static size_type const fwd_alloc_chunk_val = 1 << fwd_alloc_chunk_bits;
		static size_type const fwd_alloc_chunk_size = sizeof(char *)*fwd_alloc_chunk_val;
		linked_array<char *> chunks;

		size_type record_size;
		size_type m_size;

		__INLINE size_type GetChunkNumber(size_type k) const {return static_cast<uenum>(k) >> block_bits;}
		__INLINE size_type GetElementNumber(size_type k) const {return (k & block_bits_mask) * record_size;}
		__INLINE char * access_block(size_type k) {return chunks[GetChunkNumber(k)];}
		__INLINE const char * access_block(size_type k) const {return chunks[GetChunkNumber(k)];}
		__INLINE char & access_element(size_type k) {return access_block(k)[GetElementNumber(k)];}
		__INLINE const char & access_element(size_type k) const {return access_block(k)[GetElementNumber(k)];}
	public:
		__INLINE char & operator [] (size_type i)
		{
			assert(i < size());
			return access_element(i);
		}
		__INLINE const char & operator [] (size_type i) const
		{
			assert(i < size());
			return access_element(i);
		}
		__INLINE char & at(size_type i)
		{
			assert(i < size());
			return access_element(i);
		}
		__INLINE const char & at(size_type i) const
		{
			assert(i < size());
			return access_element(i);
		}
		void inner_resize(size_type new_size)
		{
			size_type oldnchunks2 = (static_cast<uenum>(m_size) >> block_bits) + ( (m_size & block_bits_mask) ? 1 : 0);
			size_type oldn = (static_cast<uenum>(oldnchunks2) >> fwd_alloc_chunk_bits) + ( (oldnchunks2 & fwd_alloc_chunk_bits_mask) ? 1 : 0);
			size_type newnchunks2 = (static_cast<uenum>(new_size) >> block_bits) + ( (new_size & block_bits_mask)? 1 : 0);
			size_type newn = (static_cast<uenum>(newnchunks2) >> fwd_alloc_chunk_bits) + ( (newnchunks2 & fwd_alloc_chunk_bits_mask) ? 1 : 0);
			for(size_type q = newnchunks2; q < oldnchunks2; q++)
			{
				assert(chunks[q] != NULL);
				//free(chunks[q]);
				delete [] chunks[q];
				chunks[q] = NULL;
			}
			if( newn != oldn )
			{
				if( newn > 0 )
					chunks.resize(fwd_alloc_chunk_size*newn);
				else
					chunks.clear();
			}
			for(size_type q = oldnchunks2; q < newnchunks2; q++)
			{
				assert(chunks[q] == NULL);
				chunks[q] = new char[block_size*record_size];//static_cast<char *>(malloc(block_size*record_size));
#if defined(DEBUGMEM)
				if( chunks[q] == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				assert(chunks[q] != NULL);
				std::fill(chunks[q],chunks[q]+block_size*record_size,char());
				//memset(chunks[q],0,block_size*record_size);
			}
		}
	public:

		size_type capacity() const
		{
			size_type nchunks = ((static_cast<uenum>(m_size)>>block_bits) + ((m_size & block_bits_mask)? 1 : 0));
			return nchunks*block_val*record_size;
		}
		size_type size() const {return m_size;}
		bool empty() const {return size() == 0;}
		void clear()
		{
			size_type nchunks = ((static_cast<uenum>(m_size)>>block_bits) + ((m_size & block_bits_mask)? 1 : 0));
			for(size_type q = 0; q < nchunks; q++)
			{
				//free(chunks[q]);
				delete [] chunks[q];
				chunks[q] = NULL;
			}
			chunks.clear();
			m_size = 0;
		}
		chunk_bulk_array(size_type set_record_size = 1)
		{
			record_size = set_record_size;
			m_size = 0;
		}
		chunk_bulk_array(const chunk_bulk_array & other)
		{
			record_size = other.record_size;
			m_size = 0;
			inner_resize(other.size());
			size_type nchunks = ((static_cast<uenum>(m_size)>>block_bits) + ((m_size & block_bits_mask)? 1 : 0));
			for(size_type k = 0; k < nchunks; k++) memcpy(chunks[k],other.chunks[k],block_size*record_size);
			m_size = other.size();
		}
		chunk_bulk_array & operator =(chunk_bulk_array const & other)
		{
			if( this != &other )
			{
				clear();
				record_size = other.record_size;
				inner_resize(other.size());
				size_type nchunks = ((static_cast<uenum>(m_size)>>block_bits) + ((m_size & block_bits_mask)? 1 : 0));
				for(size_type k = 0; k < nchunks; k++) memcpy(chunks[k],other.chunks[k],block_size*record_size);
				m_size = other.size();
			}
			return *this;
		}
		~chunk_bulk_array() {clear();}
		void resize(size_type n)
		{
			inner_resize(n);
			m_size = n;
		}
	};

#if defined(_OPENMP)
#define PADDING_SIZE 4096
	template<typename T>
	struct thread_private_item
	{
		T item;
		char padding[PADDING_SIZE-sizeof(T)];
		thread_private_item() :item() {}
		thread_private_item(const thread_private_item & b) :item(b.item){}
		thread_private_item & operator =(thread_private_item const & b)
		{
			item = b.item;
			return *this;
		}
	};

	/// This class is used to replace #pragma omp threadprivate
	/// Functionality that is not supported on many older systems.
	template<typename T>
	class thread_private
	{
		//std::vector< thread_private_item<T> > items;
		//T * items;
		thread_private_item<T> * items;
	public:
		thread_private()
		{
			//std::cout << "constructor " << this << std::endl;
			items = new thread_private_item<T>[omp_get_max_threads()];
			//for(int k = 0; k < omp_get_max_threads(); ++k)
			//{
			//	std::cout << (void *)&items[k] << std::endl;
			//}
		}
		thread_private(const T & b)
		{
			//std::cout << "T copy constructor " << this << std::endl;
			items = new thread_private_item<T>[omp_get_max_threads()];
			for(int k = 0; k < omp_get_max_threads(); ++k)
				items[k].item = b;
		}
		thread_private(const thread_private & b)
		{
			//std::cout << "copy constructor " << this << std::endl;
			items = new thread_private_item<T>[omp_get_max_threads()];
			for(int k = 0; k < omp_get_max_threads(); ++k)
				items[k].item = b.get(k);
		}
		~thread_private()
		{
			//std::cout << "destructor " << this << std::endl;
		}
		thread_private & operator =(thread_private const & b)
		{
			if( omp_in_parallel() )
			{
				items[omp_get_thread_num()].item = b.get();
			}
			else
			{
#pragma omp parallel
				items[omp_get_thread_num()].item = b.get();
			}
			return *this;
		}
		T & operator *() {return items[omp_get_thread_num()].item;}
		const T & operator *() const {return items[omp_get_thread_num()].item;}
		//operator T & () {return items[omp_get_thread_num()].item;}
		//operator const T & () const {return items[omp_get_thread_num()].item;}
		//operator T () {return items[omp_get_thread_num()].item;}
		//operator T () const {return items[omp_get_thread_num()].item;}
		//template <typename B>
		//T & operator = (B const & b) {items[omp_get_thread_num()].item = b; return items[omp_get_thread_num()].item;}
		//template <typename B>
		//T & operator += (B const & b) {items[omp_get_thread_num()].item += b; return items[omp_get_thread_num()].item;}
		//template <typename B>
		//T & operator -= (B const & b) {items[omp_get_thread_num()].item -= b; return items[omp_get_thread_num()].item;}
		//template <typename B>
		//T & operator *= (B const & b) {items[omp_get_thread_num()].item *= b; return items[omp_get_thread_num()].item;}
		//template <typename B>
		//T & operator /= (B const & b) {items[omp_get_thread_num()].item /= b; return items[omp_get_thread_num()].item;}
		T & get() {return items[omp_get_thread_num()].item;}
		const T & get() const {return items[omp_get_thread_num()].item;}
		T & get(int k) {return items[k].item;}
		const T & get(int k) const {return items[k].item;}
		T * operator ->() {return &items[omp_get_thread_num()].item;}
		const T * operator ->() const {return &items[omp_get_thread_num()].item;}
	};
#else //_OPENMP
	template<typename T>
	class thread_private
	{
		T item;
	public:
		thread_private() {}
		~thread_private() {}
		thread_private(const T & b) {item = b;}
		thread_private(const thread_private & b) {item = b();}
		thread_private & operator = (thread_private const & b) {item = b(); return *this;}
		T & operator *() {return item;}
		const T & operator *() const {return item;}
		//operator T & () {return item;}
		//operator const T & () const {return item;}
		//operator T () {return items[omp_get_thread_num()].item;}
		//operator T () const {return items[omp_get_thread_num()].item;}
		//template <typename B>
		//T & operator = (B const & b) {item = b; return item;}
		//template <typename B>
		//T & operator += (B const & b) {item += b; return item;}
		//template <typename B>
		//T & operator -= (B const & b) {item -= b; return item;}
		//template <typename B>
		//T & operator *= (B const & b) {item *= b; return item;}
		//template <typename B>
		//T & operator /= (B const & b) {item /= b; return item;}
		T & get() {return item;}
		const T & get() const {return item;}
		T & get(int k) {return item;}
		const T & get(int k) const {return item;}
		T * operator ->() {return &item;}
		const T * operator ->() const {return &item;}
	};
#endif //_OPENMP
	
	class memory_pool
	{
		static const unsigned pool_size_bits = 16;
		//typedef char pool_type[pool_size];
		typedef std::map<void *,unsigned> page_fault_type;
		std::vector< char * > pool;           ///< Data storage
		std::vector<bool>     inuse; ///< marks block as in use, for unordered deallocation
		std::vector<unsigned> last_alloc; ///< Last allocated block position, used to track allocation
		page_fault_type page_fault;
	public:
		unsigned last_byte() const {return last_alloc.back();}
		unsigned allocations() const {return (unsigned)(inuse.size()-1); }
		memory_pool()
		{
			//void * tmp = malloc(sizeof(char)*(1 << pool_size_bits));
			char * tmp = new char[1 << pool_size_bits];
#if defined(DEBUGMEM)
			if( tmp == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
			pool.push_back((char*)tmp); 
			last_alloc.push_back(0); 
			inuse.push_back(true); //never delete
			//std::cout << "mempool " << (void *)this << " constructor, addr " << (void *)pool.back() << std::endl;
		}
		//memory_pool(const memory_pool & b) : pool(b.pool), last_alloc(b.last_alloc) {}
		//memory_pool & operator = (memory_pool const & b) {pool = b.pool; last_alloc = b.last_alloc; return *this;}
		template<typename T>
		void * allocate(unsigned n, const T & c = T())
		{
			if( sizeof(T)*n < (1 << pool_size_bits) )
			{
				unsigned oldpos = last_alloc.size() > 0 ? last_alloc.back() : 0;
				unsigned newpos = oldpos + sizeof(T)*n;
				unsigned pageold = oldpos >> pool_size_bits;
				unsigned pagepos = newpos >> pool_size_bits;
				unsigned datapos;
				if( pagepos == pool.size() )
				{
					//std::cout << "position from " << oldpos << " to " << newpos << " need new page " << pagepos << std::endl;
					//void * tmp = malloc(sizeof(char)*(1 << pool_size_bits));
					char * tmp = new char[1 << pool_size_bits];
#if defined(DEBUGMEM)
					if( tmp == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
					pool.push_back((char*)tmp);
				}
				
				if( pagepos != pageold || last_alloc.empty() )
				{
					//std::cout << "add page " << pagepos << " start marker" << std::endl;
					last_alloc.push_back(pagepos << pool_size_bits);
					inuse.push_back(false);
					oldpos = last_alloc.back();
					newpos = oldpos + sizeof(T)*n;
				}
				datapos = oldpos%(1<<pool_size_bits);
				void * data = (void *)&pool[pagepos][datapos];
				for(unsigned i = 0; i < n; ++i) new (&static_cast<T *>(data)[i]) T(c);
				last_alloc.push_back(newpos);
				inuse.push_back(true);
				//std::cout << "allocated " << sizeof(T)*n << " bytes from " << oldpos << " to " << newpos << " at page " << pagepos << " at " << data << std::endl;
				//std::cout << this << " last_alloc[" << last_alloc.size() << "]:";
				//for(unsigned i = 0; i < last_alloc.size(); ++i) std::cout << " " << last_alloc[i] << "(" << (inuse[i] ? '+' : '-') << ")";
				//std::cout << std::endl;
				return data;
			}
			else
			{
				T * data = new T[n];
				//for(unsigned i = 0; i < n; ++i) data[i] = c;
				//void * data = malloc(sizeof(T)*n);
#if defined(DEBUGMEM)
				if( data == NULL ) {std::cout << __FILE__ << ":" << __LINE__ << "allocation returns NULL\n";}
#endif
				//for(unsigned i = 0; i < n; ++i) new (&static_cast<T *>(data)[i]) T(c);
				std::fill(data,data+n,c);
				page_fault[(void *)data] = n;
				//std::cout << "page fault for " << sizeof(T)*n << " bytes allocated at " << data << std::endl;
				//std::cout << this << " last_alloc[" << last_alloc.size() << "]:";
				//for(unsigned i = 0; i < last_alloc.size(); ++i) std::cout << " " << last_alloc[i] << "(" << (inuse[i] ? '+' : '-') << ")";
				//std::cout << std::endl;
				return (void *)data;
			}
		}
		template<typename T>
		void deallocate(T * mem)
		{
			bool find = false;
			unsigned checkpos = (unsigned)last_alloc.size(), oldpos, newpos, pagepos, datapos;
			if( checkpos > 1 )
			{
				while( !find && checkpos > 1 )
				{
					oldpos = last_alloc[checkpos-2];
					newpos = last_alloc[checkpos-1];
					pagepos = newpos >> pool_size_bits;
					datapos = oldpos % (1 << pool_size_bits);
					if( (((T*)&pool[pagepos][datapos]) == mem) )
					{
						unsigned n = (newpos - oldpos)/sizeof(T);
						for(unsigned i = 0; i < n; ++i) mem[i].~T();
						inuse[checkpos-1] = false;
						find = true;
					}
					checkpos--;
				}
			}
			//if( find )
			//	std::cout << "deallocate " << mem << " from " << oldpos << " to " << newpos << " page " << pagepos <<  std::endl;
			//else
			//	std::cout << "deallocate " << mem << " position not found " << std::endl;
			//std::cout << this << " last_alloc[" << last_alloc.size() << "]:";
			//for(unsigned i = 0; i < last_alloc.size(); ++i) std::cout << " " << last_alloc[i] << "(" << (inuse[i] ? '+' : '-') << ")";
			//std::cout << std::endl;
			if( !find )
			{
				page_fault_type::iterator it = page_fault.find((void *)mem);
				assert(it != page_fault.end() && "deallocated block does not match neither last allocated block nor page fault");
				if(it == page_fault.end() )
				{
				//	std::cout << "deallocated block does not match neither last allocated block nor page fault";
				//	throw -1;
				}
				else
				{
					//unsigned n = it->second;
					//std::cout << "deallocate page fault of " << sizeof(T)*n << " bytes at " << mem << std::endl;
					delete [] mem;
					//for(unsigned i = 0; i < n; ++i)	mem[i].~T();
					//free(mem);
					page_fault.erase(it);
				}
			}
			
			while( !inuse.empty() && inuse.back() == false )
			{
				//std::cout << "pop " << last_alloc.back() << "(" << (inuse.back() ? '+' : '-') << ") sizes " << last_alloc.size() << " " << inuse.size() << std::endl;
				inuse.pop_back();
				last_alloc.pop_back();
				//std::cout << "deallocate from " << oldpos << " to " << newpos << std::endl;
				//if( last_alloc.back() != 0 && last_alloc.back()%(1<<pool_size_bits) == 0 )
				//{
				//	std::cout << "next should remove page " << (last_alloc.back() >> pool_size_bits) << " start marker " << std::endl;
					//last_alloc.pop_back();
				//}
			}
		}
		~memory_pool()
		{
			//std::cout << "mempool destructor " << (void *)this << " in pool: ";
			//if( last_alloc.back() != 0 ) std::cout << "warning: memory pool not empty on deallocation!!!" << std::endl;
			for(unsigned k = 0; k < pool.size(); ++k)
			{
				//std::cout << (void *)pool[k] << " ";
				//free(pool[k]);
				delete [] pool[k];
			}
			//std::cout << std::endl;
			pool.clear();
			if( !page_fault.empty() )
			{
				std::cout << "warning: memory pool's page fault not empty on deallocation!!!" << std::endl;
				//for(page_fault_type::iterator it = page_fault.begin(); it != page_fault.end(); ++it)
					//free(it->first);
					//delete [] static_cast<T *>(it->first);
			}
		}
	};
	
	//static thread_private<memory_pool> _pool;
	memory_pool & get_pool();
	
	
	
	template<typename element>//, typename enumerator = unsigned int>
	class pool_array
	{
	public:
		typedef unsigned size_type;
		typedef make_unsigned<size_type>::type uenum;
		template<typename etype>
		class _iterator
		{
		private:
			etype * e;
		public:
			typedef etype * pointer;
			typedef etype & reference;
			typedef etype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_iterator():e(NULL){}
			_iterator(etype * i):e(i){}
			_iterator(const _iterator & other){e = other.e;}
			~_iterator() {};
			_iterator operator -(size_t n) const { return _iterator(e-n); }
			_iterator & operator -=(size_t n) { e-=n; return *this; }
			_iterator operator +(size_t n) const { return _iterator(e+n); }
			_iterator & operator +=(size_t n) { e+=n; return *this; }
			_iterator & operator ++(){ ++e; return *this;}
			_iterator operator ++(int){ return _iterator(e++); }
			_iterator & operator --(){ --e; return *this; }
			_iterator operator --(int){ return _iterator(e--); }
			ptrdiff_t operator -(const _iterator & other) const {return e-other.e;}
			etype & operator *() { return *e; }
			const etype & operator *() const { return *e; }
			etype * operator ->() { return e; }
			_iterator & operator =(_iterator const & other) { e = other.e; return *this; }
			bool operator ==(const _iterator & other) const { return e == other.e;}
			bool operator !=(const _iterator & other) const { return e != other.e;}
			bool operator <(const _iterator & other) const { return e < other.e;}
			bool operator >(const _iterator & other) const { return e > other.e;}
			bool operator <=(const _iterator & other) const { return e <= other.e;}
			bool operator >=(const _iterator & other) const { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
			operator const void *() const {return const_cast<const void *> (e);}
		};
		typedef _iterator<element> iterator;
		typedef _iterator<const element> const_iterator;
		template<typename etype>
		class _reverse_iterator
		{
		private:
			etype * e;
		public:
			typedef etype * pointer;
			typedef etype & reference;
			typedef etype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::random_access_iterator_tag iterator_category;
			_reverse_iterator():e(NULL){}
			_reverse_iterator(etype * i):e(i){}
			_reverse_iterator(const _reverse_iterator & other){e = other.e;}
			~_reverse_iterator() {};
			_reverse_iterator operator -(size_t n) const { return _reverse_iterator(e+n); }
			_reverse_iterator & operator -=(size_t n) { e+=n; return *this; }
			_reverse_iterator operator +(size_t n) const {return _reverse_iterator(e-n); }
			_reverse_iterator & operator +=(size_t n) { e-=n; return *this; }
			_reverse_iterator & operator ++(){ --e; return *this;}
			_reverse_iterator operator ++(int){ return _reverse_iterator(e--); }
			_reverse_iterator & operator --(){ ++e; return *this; }
			_reverse_iterator operator --(int){ return _reverse_iterator(e++); }
			ptrdiff_t operator -(const _reverse_iterator & other) const {return other.e-e;}
			etype & operator *() { return *e; }
			const etype & operator *() const { return *e; }
			etype * operator ->() { return e; }
			_reverse_iterator & operator =(_reverse_iterator const & other) { e = other.e; return *this;}
			bool operator ==(const _reverse_iterator & other) const { return e == other.e;}
			bool operator !=(const _reverse_iterator & other) const { return e != other.e;}
			bool operator <(const _reverse_iterator & other) const { return e < other.e;}
			bool operator >(const _reverse_iterator & other) const { return e > other.e;}
			bool operator <=(const _reverse_iterator & other) const { return e <= other.e;}
			bool operator >=(const _reverse_iterator & other) const { return e >= other.e;}
			operator void *() {return static_cast<void *> (e);}
			operator const void *() const {return static_cast<const void *> (e);}
		};
		typedef _reverse_iterator<element> reverse_iterator;
		typedef _reverse_iterator<const element> const_reverse_iterator;
	private:
		
		element * m_arr;
		size_type m_size;
	public:
		__INLINE element * data() {return m_arr;}
		__INLINE const element * data() const {return m_arr;}
		pool_array()
		{
			m_size = 0;
			m_arr = NULL;
		}
		pool_array(size_type n,element c = element())
		{
			m_size = n;
			m_arr = (element *)get_pool().allocate(n,c);
			assert(m_arr != NULL);
		}
		template<class InputIterator>
		pool_array(InputIterator first, InputIterator last)
		{
			//isInputForwardIterators<element,InputIterator>();
			m_size = static_cast<size_type>(std::distance(first,last));
			m_arr = (element *)get_pool().allocate(m_size,element());
			assert(m_arr != NULL);
			{
				size_type i = 0;
				InputIterator it = first;
				while(it != last) new (m_arr+i++) element(*it++);
			}
		}
		pool_array(const pool_array & other)
		{
			m_size = other.m_size;
			if( m_size )
			{
				m_arr = (element *)get_pool().allocate(m_size,element());
				assert(m_arr != NULL);
			}
			else m_arr = NULL;
			for(size_type i = 0; i < m_size; i++) new (m_arr+i) element(other.m_arr[i]);
		}
		~pool_array()
		{
			get_pool().deallocate(m_arr);
			m_arr = NULL;
			m_size = 0;
		}
		__INLINE const element & operator [] (size_type n) const
		{
			assert(n < m_size);
			return m_arr[n];
		}
		__INLINE element & operator [] (size_type n)
		{
			assert(n < m_size);
			return m_arr[n];
		}
		__INLINE const element & at (size_type n) const
		{
			assert(n < m_size);
			return m_arr[n];
		}
		__INLINE element & at (size_type n)
		{
			assert(n < m_size);
			return m_arr[n];
		}
		void resize(size_type n, element c = element())
		{
			if( m_size == 0 && m_arr == NULL )
			{
				m_arr = (element *)get_pool().allocate(n,c);
				assert(m_arr != NULL);
				m_size = n;
			}
			else if( n != m_size )
			{
				assert(false && "resize of non-empty pool_array");
			}
		}
		__INLINE element & back()
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[m_size-1];
		}
		__INLINE const element & back() const
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[m_size-1];
		}
		__INLINE element & front()
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[0];
		}
		__INLINE const element & front() const
		{
			assert(m_arr != NULL);
			assert(m_size > 0);
			return m_arr[0];
		}
		__INLINE size_type capacity() { return m_size; }
		__INLINE bool empty() const { if( m_size ) return false; return true; }
		__INLINE size_type size() const {return m_size;}
		__INLINE size_type capacity() const {return m_size;}
		void clear()
		{
			get_pool().deallocate(m_arr);
			m_arr = NULL;
			m_size = 0;
		}
		void swap(pool_array<element> & other)
		{
			element * t_m_arr = m_arr;
			size_type t_m_size = m_size;
			m_arr = other.m_arr;
			m_size = other.m_size;
			other.m_arr = t_m_arr;
			other.m_size = t_m_size;
		}
		__INLINE iterator begin() { return m_arr; }
		__INLINE iterator end() { return m_arr+m_size; }
		__INLINE const_iterator begin() const { return m_arr; }
		__INLINE const_iterator end() const { return m_arr+m_size; }
		__INLINE reverse_iterator rbegin() { return reverse_iterator(m_arr+m_size-1); }
		__INLINE reverse_iterator rend() { return reverse_iterator(m_arr-1); }
		__INLINE const_reverse_iterator rbegin() const { return const_reverse_iterator(m_arr+m_size-1); }
		__INLINE const_reverse_iterator rend() const { return const_reverse_iterator(m_arr-1); }
	};
}

#endif
