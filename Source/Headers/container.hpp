
#ifndef _CONTAINER_HPP
#define _CONTAINER_HPP


#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <memory>
#include <new>
#include <iterator>
#include <cmath>
#include <assert.h>
#include <limits>
#include "inmost_common.h"

//#define OUT_OF_RANGE

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
	
#define USE_OPTIMIZED_ARRAY_ALLOCATION
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
	public:
		__INLINE element * data() {return m_arr;}
		__INLINE const element * data() const {return m_arr;}
		array() {m_arr = NULL; m_size = 0;}
		array(size_type n,element c = element())
		{
			m_size = n;
			m_arr = static_cast<element *>(malloc(sizeof(element)*growth_formula(m_size)));
			assert(m_arr != NULL);
			for(size_type i = 0; i < m_size; i++) new (m_arr+i) element(c);
		}
		template<class InputIterator>
		array(InputIterator first, InputIterator last)
		{
			//isInputForwardIterators<element,InputIterator>();
			m_size = static_cast<size_type>(std::distance(first,last));
			m_arr = static_cast<element *>(malloc(sizeof(element)*growth_formula(m_size)));
			assert(m_arr != NULL);
			{
				size_type i = 0;
				InputIterator it = first;
				while(it != last) new (m_arr+i++) element(*it++);
			}
		}
		array(const array & other)
		{
			m_size = other.m_size;
			if( m_size ) 
			{
				m_arr = static_cast<element *>(malloc(sizeof(element)*growth_formula(m_size)));
				assert(m_arr != NULL);
			}
			else m_arr = NULL;
			for(size_type i = 0; i < m_size; i++) new (m_arr+i) element(other.m_arr[i]);
		}
		~array()
		{
			for(size_type i = 0; i < m_size; i++) m_arr[i].~element();
			if( m_arr != NULL ) free(m_arr);
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
		__INLINE element & at_safe (size_type n) 
		{
			if( n >= m_size ) resize(n+1);
			return m_arr[n];
		}
		array & operator =(array const & other)
		{
			if( this != &other )
			{
				for(size_type i = 0; i < m_size; i++) m_arr[i].~element();
				if(m_arr != NULL ) 
				{
					free(m_arr);
					m_arr = NULL;
					m_size = 0;
				}
				if( other.m_arr != NULL )
				{
					m_size = other.m_size;
					m_arr = static_cast<element *>(malloc(sizeof(element)*growth_formula(m_size)));
					assert(m_arr != NULL);
					memcpy(m_arr,other.m_arr,sizeof(element)*m_size);
				}
			}
			return *this;
		}
		void push_back(const element & e)
		{
#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
			//unoptimized variant
			if( m_size+1 > growth_formula(m_size) )
				m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*growth_formula(++m_size)));
			else ++m_size;
#else
			//optimized for current growth_formula
			if( m_size < 2 )
				m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*(++m_size)));
			else if( ((m_size+1) & (m_size-1)) == 1 )
				m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*(m_size++ << 1)));
			else m_size++;
#endif
			assert(m_arr != NULL);
			new (m_arr+m_size-1) element(e);
		}
		void pop_back()
		{
			assert(m_arr != NULL);
			m_arr[m_size--].~element();
			if( m_size > 0)
			{
#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
				//unoptimized variant
				size_type gf = growth_formula(m_size);
				if( m_size+1 > gf )
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*gf));
#else
				if( ((m_size+1) & (m_size-1)) == 1 || m_size == 1)
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*m_size));
#endif
				assert(m_arr != NULL);
			} 
			else 
			{
				free(m_arr);
				m_arr = NULL;
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
		void resize(size_type n, element c = element() )
		{
			size_type oldsize = m_size;
			m_size = n;
			for(size_type i = m_size; i < oldsize; i++) m_arr[i].~element(); //delete elements, located over the size
			if( m_size > 0 )
			{
				if( growth_formula(oldsize) != growth_formula(m_size) )
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*growth_formula(m_size)));
				assert(m_arr != NULL);
				for(size_type i = oldsize; i < m_size; i++) new (m_arr+i) element(c); //initialize extra entities
			}
			else
			{
				free(m_arr);
				m_arr = NULL;
			}
		}
		__INLINE size_type size() const {return m_size;}
		__INLINE size_type capacity() const {return growth_formula(m_size);}
		void clear() 
		{ 
			for(size_type i = 0; i < m_size; i++) m_arr[i].~element();
			m_size = 0; 
			if( m_arr ) free(m_arr); 
			m_arr = NULL; 
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
			ptrdiff_t s = iterator(m_arr+m_size-1)-pos;
			(*pos).~element();
			m_size--;
			if( m_size > 0 )
			{
				if( s > 0 ) memmove(m_arr+d,m_arr+d+1,sizeof(element)*s);
				//unoptimized variant
#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
				size_type gf = growth_formula(m_size);
				if( m_size+1 > gf )
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*gf));
#else
				if( ((m_size+1) & (m_size-1)) == 1 || m_size == 1)
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*m_size));
#endif
				assert(m_arr != NULL);
			}
			else
			{
				free(m_arr);
				m_arr = NULL;
			}
			return m_arr+d;
		}
		iterator erase(iterator b, iterator e)
		{
			ptrdiff_t d = b-begin();
			ptrdiff_t s = end()-e;
			ptrdiff_t n = e-b;
			for(iterator i = b; i != e; i++) (*i).~element();
			m_size -= n;
			if( m_size > 0 )
			{
				if( s > 0 ) memmove(m_arr+d,m_arr+d+1,sizeof(element)*s);
				size_type gf = growth_formula(m_size);
				if( m_size+n > gf )
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*gf));
				assert(m_arr != NULL);
			}
			else
			{
				free(m_arr);
				m_arr = NULL;
			}
			return m_arr+d;
		}
		iterator insert(iterator pos, const element & x)
		{
			if( static_cast<void *>(pos) == NULL)
			{
				assert(m_arr == NULL);
				pos = iterator(m_arr = static_cast<element *>(malloc(sizeof(element))));
				assert(m_arr != NULL);
			}
			ptrdiff_t d = pos-begin();
			ptrdiff_t s = end()-pos;

			//unoptimized variant
#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
			if( m_size+1 > growth_formula(m_size) )
				m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*growth_formula(++m_size)));
			else ++m_size;
#else
			//optimized for current growth_formula
			if( m_size < 2 )
				m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*(++m_size)));
			else if( ((m_size+1) & (m_size-1)) == 1 )
				m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*(m_size++ << 1)));
			else ++m_size;
#endif
			assert(m_arr != NULL);
			if( s > 0 ) memmove(m_arr+d+1,m_arr+d,sizeof(element)*s);
			new (m_arr+d) element(x);
			return m_arr+d;
		}
		void insert(iterator pos, size_type n, const element & x)
		{
			if( n > 0 )
			{
				if( static_cast<void *>(pos) == NULL)
				{
					assert(m_arr == NULL);
					pos = iterator(m_arr = static_cast<element *>(malloc(sizeof(element))));
					assert(m_arr != NULL);
				}
				ptrdiff_t d = pos-begin();
				ptrdiff_t s = end()-pos;

				if( m_size+n > growth_formula(m_size) )
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*growth_formula(m_size+n)));
				m_size+=n;

				assert(m_arr != NULL);
				if( s > 0 ) memmove(m_arr+d+n,m_arr+d,sizeof(element)*s);
				for(size_type i = 0; i < n; i++) new (m_arr+d+i) element(x);
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
					pos = iterator(m_arr = static_cast<element *>(malloc(sizeof(element))));
					assert(m_arr != NULL);
				}
				ptrdiff_t d = pos-begin();
				ptrdiff_t s = end()-pos;

				if( m_size+n > growth_formula(m_size) )
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*growth_formula(m_size+n)));
				m_size+=n;
				
				assert(m_arr != NULL);
				if( s > 0 ) memmove(m_arr+d+n,m_arr+d,sizeof(element)*s);
				{
					InputIterator it = first;
					size_type i = 0;
					while(it != last) new (m_arr+d+i++) element(*it++);
				}
			}
		}
		template <class InputIterator>
		void replace(iterator m_first, iterator m_last, InputIterator first, InputIterator last)
		{
			assert( m_size >= 0 );
			ptrdiff_t n = static_cast<ptrdiff_t>(std::distance(first,last));
			if( static_cast<void *>(m_first) == NULL && m_arr == NULL)
			{
				assert(m_arr == NULL);
				m_first = m_last = iterator(m_arr = static_cast<element *>(malloc(sizeof(element))));
				assert(m_arr != NULL);
			}
			ptrdiff_t q = m_last-m_first; 
			ptrdiff_t d = m_first-iterator(m_arr);
			ptrdiff_t s = iterator(m_arr+m_size)-m_last;
			for(iterator it = m_first; it != m_last; it++) (*it).~element();
			if( n-q != 0 )
			{
				size_type gf = growth_formula(m_size+static_cast<size_type>(n-q));
				if( gf != growth_formula(m_size) )
					m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*gf));
				m_size+=static_cast<size_type>(n-q);
				if( s > 0 ) memmove(m_arr+d+n,m_arr+d+q,sizeof(element)*s);
			}
			{
				InputIterator it = first;
				size_type i = 0;
				while(it != last) new (m_arr+d+i++) element(*it++);
			}
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
		shell() {m_arr = NULL; m_size = NULL; fixed = false;}
		shell(array<element> & arr) //dynamic
		{
			m_arr = &arr.m_arr;
			m_size = &arr.m_size;
			fixed = false;
		}
		shell(element * link, size_type size) //fixed
		{
			m_arr = &local_link;
			*m_arr = link;
			m_size = &local_size;
			*m_size = size;
			fixed = true;
		}
		shell(const shell & other)
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
#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
			//unoptimized variant
			if( (*m_size)+1 > array<element>::growth_formula(*m_size) )
				*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*array<element>::growth_formula(++(*m_size))));
			else ++(*m_size);
#else
			//optimized for current growth_formula
			if( *m_size < 2 )
				*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*(++(*m_size))));
			else if( (((*m_size)+1) & ((*m_size)-1)) == 1 )
				*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*((*m_size)++ << 1)));
			else ++(*m_size);
#endif
			assert((*m_arr) != NULL);
			new ((*m_arr)+(*m_size)-1) element(e);
		}
		void pop_back()
		{
			assert( !fixed ); // array size is fixed
			assert((*m_arr) != NULL);
			(*m_arr)[(*m_size)--].~element();
			if( (*m_size) > 0 )
			{
#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
				//unoptimized variant
				size_type gf = array<element>::growth_formula(*m_size);
				if( (*m_size)+1 > gf )
					*m_arr = static_cast<element *>(realloc(m_arr,sizeof(element)*gf));
#else
				if( (((*m_size)+1) & ((*m_size)-1)) == 1 || (*m_size) == 1)
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*(*m_size)));
#endif
				assert( (*m_arr) != NULL );
			}
			else
			{
				free(*m_arr);
				(*m_arr) = NULL;
			}
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
			for(size_type i = *m_size; i < oldsize; i++) (*m_arr)[i].~element(); //delete elements, located over the size
			if( *m_size > 0 )
			{
				if( array<element>::growth_formula(oldsize) != array<element>::growth_formula(*m_size) )
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*array<element>::growth_formula(*m_size)));
				assert( (*m_arr) != NULL );
				for(size_type i = oldsize; i < *m_size; i++) new ((*m_arr)+i) element(c); //initialize extra entities
			}
			else
			{
				free(*m_arr);
				*m_arr = NULL;
			}
		}
		__INLINE size_type size() const {return *m_size;}
		void clear() 
		{ 
			assert( !fixed ); // array size is fixed
			for(size_type i = 0; i < *m_size; i++) (*m_arr)[i].~element();
			*m_size = 0; 
			if( *m_arr ) free(*m_arr); 
			*m_arr = NULL; 
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
		__INLINE iterator end() { return *m_arr+(*m_size); }
		__INLINE const_iterator begin() const { return *m_arr; }
		__INLINE const_iterator end() const { return *m_arr+(*m_size); }
		__INLINE reverse_iterator rbegin() { return reverse_iterator(*m_arr+(*m_size)-1); }
		__INLINE reverse_iterator rend() { return reverse_iterator(*m_arr-1); }
		__INLINE const_reverse_iterator rbegin() const { return const_reverse_iterator(*m_arr+(*m_size)-1); }
		__INLINE const_reverse_iterator rend() const { return const_reverse_iterator(*m_arr-1); }
		iterator erase(iterator pos) 
		{ 
			assert( !fixed ); // array size is fixed
			ptrdiff_t d = pos-begin();
			ptrdiff_t s = iterator(*m_arr+(*m_size)-1)-pos;
			(*pos).~element();
			(*m_size)--;
			if( (*m_size) > 0 )
			{
				if( s > 0 )memmove(*m_arr+d,*m_arr+d+1,sizeof(element)*s);
#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
				size_type gf = array<element>::growth_formula(*m_size);
				if( (*m_size)+1 > gf )
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*gf));
#else
				if( (((*m_size)+1) & ((*m_size)-1)) == 1 || (*m_size) == 1)
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*(*m_size)));
#endif
				assert((*m_arr) != NULL);
			}
			else
			{
				free(*m_arr);
				*m_arr = NULL;
			}
			return (*m_arr)+d;
		}
		iterator erase(iterator b, iterator e)
		{
			assert( !fixed ); // array size is fixed
			ptrdiff_t d = b-begin();
			ptrdiff_t s = end()-e;
			ptrdiff_t n = e-b;
			for(iterator i = b; i != e; i++) (*i).~element();
			(*m_size) -= n;
			if( (*m_size) > 0 )
			{
				if( s > 0 ) memmove(*m_arr+d,*m_arr+d+n,sizeof(element)*s);	
				size_type gf = array<element>::growth_formula(*m_size);
				if( (*m_size)+n > gf )
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*gf));
				assert((*m_arr) != NULL);
			}
			else
			{
				free(*m_arr);
				*m_arr = NULL;
			}
			return (*m_arr)+d;
		}
		iterator insert(iterator pos, const element & x)
		{
			assert( !fixed ); // array size is fixed
			if( static_cast<void *>(pos) == NULL )
			{
				assert((*m_arr) == NULL);
				pos = iterator((*m_arr) = static_cast<element *>(malloc(sizeof(element))));
				assert((*m_arr) != NULL);
			}
			ptrdiff_t d = pos-begin();
			ptrdiff_t s = end()-pos;

#if !defined(USE_OPTIMIZED_ARRAY_ALLOCATION)
			//unoptimized variant
			if( (*m_size)+1 > array<element>::growth_formula(*m_size) )
				*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*array<element>::growth_formula(++(*m_size))));
			else ++(*m_size);
#else
			//optimized for current growth_formula
			if( *m_size < 2 )
				*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*(++(*m_size))));
			else if( (((*m_size)+1) & ((*m_size)-1)) == 1 )
				*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*((*m_size)++ << 1)));
			else ++(*m_size);
#endif
			assert((*m_arr) != NULL);
			if( s > 0 ) memmove((*m_arr)+d+1,(*m_arr)+d,sizeof(element)*s);
			new ((*m_arr)+d) element(x);
			return (*m_arr)+d;
		}
		void insert(iterator pos, size_type n, const element & x)
		{
			assert( !fixed ); // array size is fixed
			if( n > 0 )
			{
				if( static_cast<void *>(pos) == NULL)
				{
					assert((*m_arr) == NULL);
					pos = iterator((*m_arr) = static_cast<element *>(malloc(sizeof(element))));
					assert((*m_arr) != NULL);
				}
				ptrdiff_t d = pos-iterator(*m_arr);
				ptrdiff_t s = iterator((*m_arr)+(*m_size))-pos;

				if( (*m_size)+n > array<element>::growth_formula(*m_size) )
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*array<element>::growth_formula((*m_size)+n)));
				(*m_size)+=n;
				
				assert((*m_arr) != NULL);
				if( s > 0 ) memmove((*m_arr)+d+n,(*m_arr)+d,sizeof(element)*s);
				for(size_type i = 0; i < n; i++) new ((*m_arr)+d+i) element(x);
			}
		}
		template <class InputIterator>
		void insert(iterator pos, InputIterator first, InputIterator last)
		{
			assert( !fixed ); // array size is fixed
			ptrdiff_t n = static_cast<ptrdiff_t>(std::distance(first,last));
			if( n > 0 )
			{
				if( static_cast<void *>(pos) == NULL)
				{
					assert((*m_arr) == NULL);
					pos = iterator((*m_arr) = static_cast<element *>(malloc(sizeof(element))));
					assert((*m_arr) != NULL);
				}
				ptrdiff_t d = pos-iterator(*m_arr);
				ptrdiff_t s = iterator((*m_arr)+(*m_size))-pos;


				if( (*m_size)+n > array<element>::growth_formula(*m_size) )
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*array<element>::growth_formula((*m_size)+static_cast<size_type>(n))));
				(*m_size)+=static_cast<size_type>(n);

				assert((*m_arr) != NULL);
				if( s > 0 ) memmove((*m_arr)+d+n,(*m_arr)+d,sizeof(element)*s);
				{
					InputIterator it = first;
					size_type i = 0;
					while(it != last) new ((*m_arr)+d+i++) element(*it++);
				}
			}
		}
		template <class InputIterator>
		void replace(iterator m_first, iterator m_last, InputIterator first, InputIterator last)
		{
			ptrdiff_t n = static_cast<ptrdiff_t>(std::distance(first,last));
			if( static_cast<void *>(m_first) == NULL)
			{
				assert((*m_arr)==NULL);
				m_first = m_last = iterator((*m_arr) = static_cast<element *>(malloc(sizeof(element))));
				assert((*m_arr)!=NULL);
			}
			ptrdiff_t q = m_last-m_first; 
			ptrdiff_t d = m_first-iterator(*m_arr);
			ptrdiff_t s = iterator((*m_arr)+(*m_size))-m_last;
			for(iterator it = m_first; it != m_last; it++) (*it).~element();
			if( n-q != 0 )
			{
				assert( !fixed ); // array size is fixed
				size_type gf = array<element>::growth_formula((*m_size)+static_cast<size_type>(n-q));
				if( gf != array<element>::growth_formula(*m_size) )
					*m_arr = static_cast<element *>(realloc(*m_arr,sizeof(element)*gf));
				(*m_size)+=static_cast<size_type>(n-q);
			}
			if( s > 0 ) 
				memmove((*m_arr)+d+n,(*m_arr)+d+q,sizeof(element)*s);
			{
				InputIterator it = first;
				size_type i = 0;
				while(it != last) new ((*m_arr)+d+i++) element(*it++);
			}
		}
	};


	
	template <typename IndType,typename ValType>	
	class sparse_data
	{
	public:
		typedef unsigned enumerate;
		static const enumerate prealloc = 4;
		typedef struct pair_t
		{
			IndType first;
			ValType second; 
			pair_t() :first(0),second(0.0) {}
			pair_t(IndType first, ValType second) :first(0), second(0.0) {}
		} pair;
		typedef pair * iterator;
		typedef const iterator const_iterator;
	private:
		static int comparator(const void * pa, const void * pb)
		{
			pair * a = (pair *)pa;
			pair * b = (pair *)pb;
			return a->first - b->first;
		}
		pair * array;
		enumerate arr_size;
		enumerate arr_alloc;
		void test_allocate()
		{
			if( arr_size > arr_alloc )
			{
				enumerate old_arr_alloc = arr_alloc;
				while(arr_size > arr_alloc) arr_alloc = arr_alloc << 1;
				array = static_cast<pair *>(realloc(array,arr_alloc*sizeof(pair)));
				assert(array != NULL);
				for (enumerate i = old_arr_alloc; i < arr_alloc; i++) array[i].first = std::numeric_limits<IndType>::max();
				//memset(array+old_arr_alloc,0xff,sizeof(pair)*(arr_alloc-old_arr_alloc));
			}
		}
	public:
		void swap(sparse_data<IndType, ValType> & other)
		{
			pair * tmp = array;
			array = other.array;
			other.array = tmp;
			enumerate itmp = arr_size;
			arr_size = other.arr_size;
			other.arr_size = itmp;
			itmp = arr_alloc;
			arr_alloc = other.arr_alloc;
			other.arr_alloc = itmp;
		}
		void reserve(enumerate size)
		{
			enumerate new_alloc = 1;
			while( new_alloc < size ) new_alloc = new_alloc << 1;
			array = static_cast<pair *>(realloc(array,new_alloc*sizeof(pair)));
			assert(array != NULL);
			for (enumerate i = arr_alloc; i < new_alloc; i++) array[i].first = std::numeric_limits<IndType>::max();
			//memset(array+arr_alloc,0xff,sizeof(pair)*(new_alloc-arr_alloc));
			arr_alloc = new_alloc;
		}
		iterator lower_bound(IndType ind)
		{
			/*
			if( arr_size < 16 )
			{	
				for(enumerate i = 0; i < arr_size; i++)
					if( array[i].first >= ind ) return array+i;
				return array+arr_size;
			}
			*/			
			unsigned k = 0;
			for(unsigned b = arr_alloc >> 1; b ; b = b >> 1)
			{
				unsigned j = k | b;
				if( array[j].first <= ind ) k = j;
			}
			if( array[k].first < ind ) k++;	
			
			return array+k;
		}
		iterator find(IndType ind)
		{
			iterator k = lower_bound(ind);
			if( k == end() ) return end();
			if( k->first == ind ) return k;
			return end();
		}
		iterator insert(iterator pos, const IndType & x)
		{
			assert(pos == end() || x < pos->first);//check here that we don't break the order
			ptrdiff_t d = pos-array;
			ptrdiff_t s = arr_size-d;
			arr_size++;
			test_allocate();
			if( s ) memmove(array+d+1,array+d,sizeof(pair)*s);
			(array+d)->first = x;
			new (&(array[d].second)) ValType();
			return array+d;
		}
		void push_back(const pair & in)
		{
			assert(arr_size == 0 || in.first < (array+arr_size-1)->first);//check here that we don't break the order
			arr_size++;
			test_allocate();
			(array+arr_size-1)->first = in.first;
			(array+arr_size-1)->second = in.second;
		}
		iterator erase(iterator pos)
		{ 
			ptrdiff_t d = pos-array;
			ptrdiff_t s = (array+arr_size-1)-pos;
			(pos->second).~ValType();
			memmove(array+d,array+d+1,sizeof(pair)*s);
			arr_size--;
			return array+d;
		}
		sparse_data(pair * first, pair * last)
		{
			arr_size = last-first;
			arr_alloc = static_cast<enumerate>(prealloc);
			if( arr_size <= arr_alloc )
				array = static_cast<pair *>(malloc(arr_alloc*sizeof(pair)));
			else
			{
				array = NULL;
				test_allocate();
			}
			assert(array != NULL);
			memcpy(array,first,sizeof(pair)*arr_size);
			for (enumerate i = arr_alloc; i < arr_size; i++) array[i].first = std::numeric_limits<IndType>::max();
			//memset(array+arr_size,0xff,sizeof(pair)*(arr_alloc-arr_size));
			bool need_sort = false;
			for(enumerate k = 1; k < arr_size; k++)
				if( array[k].first < array[k-1].first )
				{
					need_sort = true;
					break;
				}
			if( need_sort ) qsort(array,sizeof(pair),arr_size,comparator);
		}
		sparse_data()
		{
			arr_size = 0;
			arr_alloc = static_cast<enumerate>(prealloc);
			array = static_cast<pair *>(malloc(sizeof(pair)*arr_alloc));
			assert(array != NULL);
			for (enumerate i = 0; i < arr_alloc; i++) array[i].first = std::numeric_limits<IndType>::max();
			//memset(array,0xff,sizeof(pair)*arr_alloc);
		}
		sparse_data(const sparse_data & other)
		{
			arr_size = other.arr_size;
			arr_alloc = other.arr_alloc;
			array = static_cast<pair *>(malloc(arr_alloc*sizeof(pair)));
			assert(array != NULL);
			memcpy(array,other.array,other.arr_alloc*sizeof(pair));
		}
		~sparse_data()
		{
			for(iterator i = begin(); i != end(); i++) (i->second).~ValType();
			free(array);
			arr_size = arr_alloc = 0;
		}
		sparse_data & operator =(sparse_data const & other)
		{
			if( &other != this )
			{
				for(iterator i = begin(); i != end(); i++) (i->second).~ValType();
				arr_size = other.arr_size;
				arr_alloc = other.arr_alloc;
				array = static_cast<pair *>(realloc(array,arr_alloc*sizeof(pair)));
				assert(array != NULL);
				memcpy(array,other.array,arr_alloc*sizeof(pair));
			}
			return *this;
		}
		ValType & operator [](IndType row)
		{
			iterator q = lower_bound(row);
			if( q != end() && q->first == row ) return q->second;
			return insert(q,row)->second;
		}
		ValType operator [](IndType row) const
		{
			iterator q = lower_bound(row);
			assert(q != end() && q->first == row);
			return q->second;
		}
		enumerate size() const { return arr_size; }
		bool empty() const {return size() == 0; }
		iterator begin() { return array; }
		const_iterator begin() const { return array; }
		iterator end() { return array+arr_size; }
		const_iterator end() const { return array+arr_size; }
		void clear()
		{
			for(iterator i = begin(); i != end(); i++) (i->second).~ValType();
			for (enumerate i = 0; i < arr_size; i++) array[i].first = std::numeric_limits<IndType>::max();
			//memset(array,0xff,sizeof(pair)*arr_size);
			arr_size = 0;	
		}
		enumerate capacity() {return arr_alloc;}
		
	};
	
	
	
	
	template<typename IndType,typename ValType>
	class interval
	{
	public:
		typedef ValType * iterator;
		typedef ValType const * const_iterator;
		//typedef const iterator const_iterator;
	private:
		ValType * array;
		IndType beg_index, end_index;
	public:
		void clear()
		{
			for (IndType i = beg_index; i < end_index; i++) (array[i]).~ValType();
			if (beg_index != end_index) free(array + beg_index);
			array = NULL;
			beg_index = end_index = 0;
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
				ValType * tmp = array;
				array = other.array;
				other.array = tmp;
			}
		}
		interval()
		{
			beg_index = 0;
			end_index = 0;
			array = NULL;//static_cast<ValType *>(malloc(sizeof(ValType)*(end_index-beg_index)));
			//assert(array != NULL);
			//array = array - beg_index;
			//for(IndType i = beg_index; i != end_index; ++i) new (array+i) ValType();
		}
		interval(IndType beg)
		{
			beg_index = beg;
			end_index = beg_index;//+2;
			array = NULL; //static_cast<ValType *>(malloc(sizeof(ValType)*(end_index-beg_index)));
			//assert(array != NULL);
			//array = array - beg_index;
			//for(IndType i = beg_index; i != end_index; ++i) new (array+i) ValType();
		}
		interval(IndType beg, IndType end, ValType c = ValType())
		{
			beg_index = beg;
			end_index = end;
			if (beg != end)
			{
				array = static_cast<ValType *>(malloc(sizeof(ValType)*(end_index - beg_index)));
				assert(array != NULL);
				array = array - beg_index;
				for (IndType i = beg_index; i < end_index; ++i) new (array + i) ValType(c);

				//std::cout << __FUNCTION__ << " address " << array << std::endl;
			}
			else array = NULL;
		}
		interval(const interval & other)
		{
			//std::cout << __FUNCTION__ << std::endl;
			beg_index = other.beg_index;
			end_index = other.end_index;
			if( beg_index != end_index )
			{
				array = static_cast<ValType *>(malloc(sizeof(ValType)*(end_index-beg_index)));
				assert(array != NULL);
				array = array - beg_index;
				for(IndType i = beg_index; i < end_index; ++i) 
				{
					new (array+i) ValType(other.array[i]);
				}
				//std::cout << __FUNCTION__ << " address " << array << std::endl;
			}
			else array = NULL;
		}
		~interval()
		{
			//std::cout << __FUNCTION__ << " delete address " << array << std::endl;
			for(IndType i = beg_index; i < end_index; i++) (array[i]).~ValType();
			if( beg_index != end_index ) free(array+beg_index);
			array = NULL;
		}
		interval & operator =(interval const & other)
		{
			if( &other != this )
			{
				for(iterator i = begin(); i != end(); ++i) (*i).~ValType();
				beg_index = other.beg_index;
				end_index = other.end_index;
				if( beg_index != end_index )
				{
					array = static_cast<ValType *>(realloc(array+beg_index,sizeof(ValType)*(end_index-beg_index)));
					assert(array != NULL);
					array = array - beg_index;
					for(IndType i = beg_index; i < end_index; ++i) new (array+i) ValType(other.array[i]);
					//std::cout << __FUNCTION__ << " address " << array << std::endl;
				}
				else 
				{
					free(array+beg_index);
					array = NULL;
				}
			}
			return *this;
		}
		ValType & at(IndType row)
		{
			assert(row >= beg_index);
			assert(row < end_index);
			return array[row];
		}
		const ValType & at(IndType row) const
		{
			assert(row >= beg_index);
			assert(row < end_index);
			return array[row];
		}
		ValType & operator [](IndType row)
		{
			//std::cout << "pos: " << row << std::endl;
			/*
			if( row >= end_index )
			{
				IndType new_end_index = 1;
				IndType temp = row-beg_index;
				while( new_end_index <= temp ) new_end_index = new_end_index << 1;
				new_end_index += beg_index;
				//std::cout << "end: " << end_index << " new end: " << new_end_index << std::endl;
				array = static_cast<ValType *>(realloc(array,sizeof(ValType)*(new_end_index-beg_index)));
				IndType end = new_end_index-beg_index;
				for(IndType i = end_index-beg_index; i != end; ++i) new (array+i) ValType();
				end_index = new_end_index;
			}
			if( row >= last_index ) last_index = row+1;
			*/
			assert(row >= beg_index );
			assert(row < end_index );
			return array[row];
		}
		const ValType & operator [](IndType row) const
		{
			assert(row >= beg_index );
			assert(row < end_index );
			return array[row];
		}
		void set_interval_beg(IndType beg)
		{
			IndType shift = beg-beg_index;
			shift_interval(shift);
		}
		void set_interval_end(IndType end)
		{
			if( end == end_index ) return;
			if( beg_index != end )
			{
				ValType * array_new = static_cast<ValType *>(malloc(sizeof(ValType)*(end-beg_index)));
				assert(array_new != NULL);
				array_new = array_new - beg_index;
				for(IndType i = beg_index; i < std::min(end,end_index); ++i) new (array_new+i) ValType(array[i]);
				for(IndType i = end_index; i < end; ++i) new (array_new+i) ValType();
				for(IndType i = beg_index; i < end_index; ++i) array[i].~ValType();

				if( array != NULL ) free(array+beg_index);
				array = array_new;
			}
			else
			{
				free(array+beg_index);
				array = NULL;
			}
			end_index = end;
		}
		
		void shift_interval(IndType shift)
		{
			array = array + beg_index;
			beg_index += shift;
			end_index += shift;
			array = array - beg_index;
		}
		iterator begin() {return array+beg_index;}
		const_iterator begin() const {return array+beg_index;}
		//const_iterator begin() {return array;}
		iterator end() {return array + end_index;}
		const_iterator end() const {return array + end_index;}
		//const_iterator end() {return array + (end_index-end_index);}
		IndType get_interval_beg() const { return beg_index; }
		IndType get_interval_end() const { return end_index; }
		int size() const {return end_index - beg_index;}
		bool empty() const {return beg_index == end_index;}
	};
	
	//this version is safe for std::map
	/*
	template<typename IndType,typename ValType>
	class interval
	{
	public:
		typedef ValType * iterator;
		typedef ValType const * const_iterator;
	private:
		ValType * array;
		IndType beg_index, end_index, last_index;
	public:
		interval()
		{
			beg_index = 0;
			last_index = beg_index;
			end_index = 0;
			array = NULL;
		}
		interval(IndType beg)
		{
			beg_index = beg;
			last_index = beg_index;
			end_index = beg_index;
			array = NULL;
		}
		interval(IndType beg, IndType end)
		{
			beg_index = beg;
			last_index = end;
			end_index = end;
			if( end_index-beg_index > 0 )
			{
				array = static_cast<ValType *>(malloc(sizeof(ValType)*(end_index-beg_index)));
				assert(array != NULL);
				IndType cycle_end = end_index-beg_index;
				for(IndType i = 0; i != cycle_end; ++i) new (array+i) ValType();
			}
			else array = NULL;
		}
		interval(const interval & other)
		{
			beg_index = other.beg_index;
			last_index = other.last_index;
			end_index = other.end_index;
			array = static_cast<ValType *>(malloc(sizeof(ValType)*(end_index-beg_index)));
			assert(array != NULL);
			IndType end = end_index-beg_index;
			for(IndType i = 0; i != end; ++i) 
			{
				//std::cout << this << " " << __FILE__  << ":" << __LINE__ << " call constructor " << i << " obj " << array+i << std::endl;
				new (array+i) ValType(other.array[i]);
			}
		}
		~interval()
		{
			//for(iterator i = begin(); i != end(); ++i) (*i).~ValType();
			for(IndType i = 0; i < end_index-beg_index; i++) (array[i]).~ValType();
			free(array);
		}
		interval & operator =(interval const & other)
		{
			if( &other != this )
			{
				for(iterator i = begin(); i != end(); ++i) (*i).~ValType();
				beg_index = other.beg_index;
				last_index = other.last_index;
				end_index = other.end_index;
				array = static_cast<ValType *>(realloc(array,sizeof(ValType)*(end_index-beg_index)));
				assert(array != NULL);
				IndType end = end_index-beg_index;
				for(IndType i = 0; i != end; ++i) 
				{
					//std::cout << this << " " << __FILE__  << ":" << __LINE__ << " call constructor " << i << " obj " << array+i << std::endl;
					new (array+i) ValType(other.array[i]);
				}
			}
			return *this;
		}
		ValType & operator [](IndType row)
		{
			//std::cout << "pos: " << row << std::endl;
			if( row >= end_index )
			{
				IndType new_end_index = 1;
				IndType temp = row-beg_index;
				while( new_end_index <= temp ) new_end_index = new_end_index << 1;
				new_end_index += beg_index;
				ValType * array_new = static_cast<ValType *>(malloc(sizeof(ValType)*(new_end_index-beg_index)));
				assert(array_new != NULL);
				for(IndType i = 0; i != end_index-beg_index; ++i) 
				{
					new (array_new+i) ValType(*(array+i));
					(*(array+i)).~ValType();
				}
				IndType end = new_end_index-beg_index;
				for(IndType i = end_index-beg_index; i != end; ++i) 
				{
					//std::cout << this << " " << __FILE__  << ":" << __LINE__ << " call constructor " << i << " obj " << array+i << std::endl;
					new (array_new+i) ValType();
				}
				end_index = new_end_index;
				free(array);
				array = array_new;
			}
			if( row >= last_index ) last_index = row+1;
			return array[row-beg_index];
		}
		const ValType & operator [](IndType row) const
		{
			assert(row >= beg_index );
			assert(row < end_index );
			return array[row-beg_index];
		}
		bool empty() const {return beg_index == last_index;}
		void set_interval_beg(IndType beg)
		{
			IndType shift = beg-beg_index;
			shift_interval(shift);
		}
		void set_interval_end(IndType end)
		{
			
			if( end > end_index )
			{
				ValType * array_new = static_cast<ValType *>(malloc(sizeof(ValType)*(end-beg_index)));
				assert(array_new != NULL);
				IndType cycle_end = end_index-beg_index;
				for(IndType i = 0; i != cycle_end; ++i) 
				{
					new (array_new+i) ValType(*(array+i));
					(*(array+i)).~ValType();
				}
				cycle_end = end-beg_index;
				for(IndType i = end_index-beg_index; i != cycle_end; ++i) 
				{
					//std::cout << this << " " << __FILE__  << ":" << __LINE__ << " call constructor " << i << " obj " << array+i << std::endl;
					new (array_new+i) ValType();
				}
				end_index = end;
				free(array);
				array = array_new;
			}
			last_index = end;
		}
		
		void shift_interval(IndType shift)
		{
			beg_index += shift;
			end_index += shift;
			last_index += shift;
		}
		iterator begin() {return array;}
		const_iterator begin() const {return array;}
		iterator end() {return array + (last_index-beg_index);}
		
		const_iterator end() const {return array + (last_index-beg_index);}
		IndType get_interval_beg() const { return beg_index; }
		IndType get_interval_end() const { return last_index; }
		int size() const {return last_index - beg_index;}
	};
	*/
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
				pbegin = static_cast<element *>(malloc(sizeof(element)*n));
				assert(pbegin != NULL);
				pend = pbegin+n;
				preserved = pbegin+n;
			}
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
				for(size_type i = n; i < k; i++) pbegin[i].~element();
				if( pbegin == stack )
				{
					pbegin = static_cast<element *>(malloc(sizeof(element)*n));
					assert(pbegin != NULL);
					for(size_type i = 0; i < k; i++) 
					{
						new (pbegin+i) element(stack[i]);
						stack[i].~element();
					}
				}
				else
				{
					element * pbegin_new = static_cast<element *>(malloc(sizeof(element)*n));
					assert(pbegin_new != NULL);
					for(size_type i = 0; i < k; i++) 
					{
						new (pbegin_new+i) element(pbegin[i]);
						pbegin[i].~element();
					}
					free(pbegin);
					pbegin = pbegin_new;
				}
				pend = pbegin+ (k < n ? k : n);
				preserved = pbegin + n;
			}
			/*
			else if( pbegin != stack )
			{
				memcpy(stack,pbegin,sizeof(element)*n);
				for(size_type i = n; i < k; i++) pbegin[i].~element();
				free(pbegin);
				pbegin = stack;
				pend = stack+n;
				preserved = stack+static_cast<size_type>(stacked);
			}
			*/
		}
		dynarray()
		{
			pbegin = pend = stack;
			preserved = stack+static_cast<size_type>(stacked);
		}
		dynarray(size_type n,element c = element())
		{
			preallocate(n);
			for(element * i = pbegin; i < pend; i++) new (i) element(c);
		}
		template<class InputIterator>
		dynarray(InputIterator first, InputIterator last)
		{
			size_type n = static_cast<size_type>(std::distance(first,last));
			preallocate(n);
			{
				InputIterator it = first;
				element * i = pbegin;
				while(it != last) {new (i++) element(*(it++));}
			}
		}
		dynarray(const dynarray & other)
		{
			//std::cout << __FUNCTION__ << std::endl;
			size_type n = other.size();
			preallocate(n);
			for(size_type k = 0; k < n; k++)
				new (pbegin+k) element(other.pbegin[k]);
		}
		
		~dynarray()
		{
			for(element * i = pbegin; i < pend; i++) (*i).~element();
			if( pbegin != stack ) free(pbegin);
		}
		dynarray & operator =(dynarray const & other)
		{
			if( this != &other )
			{
				size_type n = size();
				for(element * i = pbegin; i != pend; ++i) (*i).~element();
				if(pbegin != stack) free(pbegin);
				n = other.size();
				preallocate(n);
				for(size_type k = 0; k < n; k++)
					new (pbegin+k) element(other.pbegin[k]);
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
			new (pend++) element(e);	
		}
		void pop_back()
		{
			(*(pend--)).~element();
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
			while( capacity() < n ) reserve(capacity()*2);
			for(element * i = pbegin+n; i < pbegin+oldsize; i++) (*i).~element(); //delete elements, located over the size
			for(element * i = pbegin+oldsize; i < pbegin+n; i++) new (i) element(c); //initialize extra entities
			pend = pbegin + n;
		}
		__INLINE size_type size() const {return static_cast<size_type>(pend-pbegin);}
		void clear() 
		{ 
			for(element * i = pbegin; i < pend; i++) (*i).~element();
			pend = pbegin;
			//pbegin = pend = stack;
			//preserved = stack+static_cast<size_type>(stacked);
		}
		void move(dynarray<element,stacked> & other)
		{
			size_type k = size(), n = other.size();
			if( n > static_cast<size_type>(stacked) ) free(other.pbegin);
			if( k > static_cast<size_type>(stacked) )
			{
				other.pbegin = pbegin;
				other.pend = pend;
				other.preserved = preserved;
			}
			else
			{
				memcpy((void *)other.stack,(void *)stack,sizeof(element)*k);
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
			if( k <= static_cast<size_type>(stacked) && n <= static_cast<size_type>(stacked) )
			{
				char temp[stacked*sizeof(element)];
				memcpy(temp,(void *)stack,sizeof(element)*k);
				memcpy((void *)stack,(void *)other.stack,sizeof(element)*n);
				memcpy((void *)other.stack,temp,sizeof(element)*k);
				other.pend = other.pbegin+k;
				pend = pbegin+n;
			}
			else if( k <= static_cast<size_type>(stacked) && n > static_cast<size_type>(stacked) )
			{
				memcpy((void *)other.stack,(void *)stack,sizeof(element)*k);
				pbegin = other.pbegin;
				pend = other.pend;
				preserved = other.preserved;
				other.pbegin = other.stack;
				other.pend = other.stack+k;
				other.preserved = other.stack+static_cast<size_type>(stacked);
			}
			else if( k > static_cast<size_type>(stacked) && n <= static_cast<size_type>(stacked) )
			{
				memcpy((void *)stack,(void *)other.stack,sizeof(element)*n);
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
			(*pos).~element();
			if( s > 0 ) memmove((void *)(pbegin+d),(void *)(pbegin+d+1),sizeof(element)*s);
			pend--;
			return pbegin+d;
		}
		iterator erase(iterator b, iterator e)
		{
			ptrdiff_t d = b-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-e;
			ptrdiff_t n = e-b;
			for(iterator i = b; i != e; i++) (*i).~element();
			if( s > 0 ) memmove((void *)(pbegin+d),(void *)(pbegin+d+1),sizeof(element)*s);
			pend -= n;
			return pbegin+d;
		}
		iterator insert(iterator pos, const element & x)
		{
			ptrdiff_t d = pos-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-pos;
			if( pend == preserved ) reserve(capacity()*2);
			pend++;
			if( s > 0 ) memmove((void *)(pbegin+d+1),(void *)(pbegin+d),sizeof(element)*s);
			new (pbegin+d) element(x);
			return pbegin+d;
		}
		void insert(iterator pos, size_type n, const element & x)
		{
			ptrdiff_t d = pos-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-pos;
			while( capacity()-size() < n ) reserve(capacity()*2);
			if( s > 0 ) memmove((void *)(pbegin+d+n),(void *)(pbegin+d),sizeof(element)*s);
			pend+=n;
			for(size_type i = 0; i < n; i++) new (pbegin+d+i) element(x);
		}
		template <class InputIterator>
		void insert(iterator pos, InputIterator first, InputIterator last)
		{
			ptrdiff_t n = static_cast<ptrdiff_t>(std::distance(first,last));
			ptrdiff_t d = pos-iterator(pbegin);
			ptrdiff_t s = iterator(pend)-pos;
			while( capacity() < size()+n ) reserve(capacity()*2);
			if( s > 0 ) memmove((void *)(pbegin+d+n),(void *)(pbegin+d),sizeof(element)*s);
			{
				InputIterator it = first;
				element * i = pbegin+d;
				while(it != last) new (i++) element(*it++);
			}
			pend+=n;
		}
	};
	
	template<class key,class value, unsigned int stacked>
	class tiny_map
	{
	public:
		typedef dynarray<std::pair<key,value>,stacked> container;
		typedef typename dynarray<std::pair<key,value>,stacked>::size_type size_type;
	private:
		container inner_data;
	public:
		class iterator : public container::iterator 
		{
			public: 
				iterator() : container::iterator() {}
				iterator( const typename container::iterator & other) :container::iterator(other) {}
		};
		class const_iterator : public container::const_iterator 
		{
			public: 
				const_iterator() : container::const_iterator() {}
				const_iterator( const typename container::const_iterator & other) :container::const_iterator(other) {}
		};
		tiny_map() :inner_data() {}
		tiny_map(const tiny_map & other) :inner_data(other.inner_data) {}
		tiny_map & operator = (tiny_map const & other) {inner_data = other.inner_data; return *this;}
		~tiny_map() {}
		const_iterator find(const key & x) const
		{
			for(const_iterator it = inner_data.begin(); it != inner_data.end(); ++it)
				if( it->first == x ) return it;
			return end();
		}
		iterator find(const key & x)
		{
			for(iterator it = inner_data.begin(); it != inner_data.end(); ++it)
				if( it->first == x ) return it;
			return end();
		}
		iterator begin() {return inner_data.begin();}
		iterator end() {return inner_data.end();}
		const_iterator begin() const {return inner_data.begin();}
		const_iterator end() const {return inner_data.end();}
		value & operator [](const key & x)
		{
			for(typename container::iterator it = inner_data.begin(); it != inner_data.end(); ++it)
				if( it->first == x ) return it->second;
			inner_data.push_back(std::pair<key,value>(x,value()));
			return inner_data.back().second;
		}
		size_type size() {return inner_data.size();}
		void clear() {inner_data.clear();}
		bool empty() const {return inner_data.empty();}
		iterator erase(iterator pos) {return inner_data.erase(typename container::iterator(pos));}
		
	};
	
	
	template<typename IndType, typename TValue, int HSize>
	class small_hash
	{
	private:
		typedef sparse_data<IndType,TValue> inner_class;
		typedef typename inner_class::pair_t inner_type;
		typedef typename inner_class::iterator inner_iter;
		typedef inner_class inner_class_array[HSize];
		//typedef typename inner_class::reverse_iterator reverse_inner_iter;
		inner_class_array lists;
		IndType compute_pos(IndType key) {return (key*15637)%HSize;}
	public:
		template<typename dtype>
		class _iterator
		{
		private:
			inner_class_array *lists;
			int cur_list;
			inner_iter it;
		public:
			typedef dtype * pointer;
			typedef dtype & reference;
			typedef dtype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::bidirectional_iterator_tag iterator_category;
			_iterator(int cur_list, inner_iter it, inner_class_array * lists) :cur_list(cur_list), it(it), lists(lists) {}
			_iterator():cur_list(-1),it(NULL), lists(NULL){}
			_iterator(const _iterator & other){it = other.it; cur_list = other.cur_list; lists = other.lists;}
			~_iterator() {};
			_iterator & operator ++()
			{ 
				++it; 
				if(it == (*lists)[cur_list].end() ) 
				{
					do
					{
						++cur_list;
					} while ((*lists)[cur_list].empty() && cur_list < HSize);
					if( cur_list < HSize ) 
						it = (*lists)[cur_list].begin(); 
					else 
						it = NULL;
				} 
				return *this;
			}
			_iterator operator ++(int){ _iterator ret = *this; ++(*this); return ret;}
			_iterator & operator --()
			{ 
				--it; 
				if(it == (*lists)[cur_list].begin()-1 ) 
				{
					do
					{
						--cur_list;
					} while((*lists)[cur_list].empty() && cur_list >= 0 ); 
					if( cur_list >= 0 ) 
						it = (*lists)[cur_list].end()-1; 
					else 
						it = NULL;
				}  
				return *this; 
			}
			_iterator operator --(int){ _iterator ret = *this; --(*this); return ret; }
			dtype & operator *() { return *it; }
			dtype * operator ->() { return &*it; }
			_iterator & operator =(_iterator const & other) {it = other.it; cur_list = other.cur_list; lists = other.lists; return *this; }
			bool operator ==(const _iterator & other) const { return it == other.it;}
			bool operator !=(const _iterator & other) const { return it != other.it;}
			bool operator <(const _iterator & other) const { return it < other.it;}
			bool operator >(const _iterator & other) const { return it > other.it;}
			bool operator <=(const _iterator & other) const { return it <= other.it;}
			bool operator >=(const _iterator & other) const { return it >= other.it;}
		};
		typedef _iterator<inner_type> iterator;
		typedef _iterator<const inner_type> const_iterator;
		iterator begin() 
		{
			int i = 0; 
			while(lists[i].empty() && i < HSize) i++; 
			return iterator(i,i < HSize ? lists[i].begin() : NULL,&lists);
		}
		iterator end() {return iterator(HSize,NULL,&lists);}
		/*
		template<typename dtype>
		class _reverse_iterator
		{
		private:
			inner_class_array * lists;
			int cur_list;
			inner_iter it;
		public:
			typedef dtype * pointer;
			typedef dtype & reference;
			typedef dtype value_type;
			typedef ptrdiff_t difference_type;
			typedef std::bidirectional_iterator_tag iterator_category;
			_reverse_iterator(int cur_list, inner_iter it, inner_class_array * lists) :cur_list(cur_list), it(it), lists(lists) {}
			_reverse_iterator():cur_list(-1),it(NULL), lists(NULL){}
			_reverse_iterator(const _reverse_iterator & other){it = other.it; cur_list = other.cur_list; lists = other.lists;}
			~_reverse_iterator() {};
			_reverse_iterator & operator ++(){ ++it; if(it == (*lists)[cur_list]->rend() ) {--cur_list; if( cur_list >= 0 ) it = (*lists)[cur_list]->rbegin(); else it = NULL;} return *this;}
			_reverse_iterator operator ++(int){ _iterator ret = *this; ++(*this); return ret; }
			_reverse_iterator & operator --(){ --it; if(it == (*lists)[cur_list]->rbegin()-1 ) {++cur_list; if( cur_list < HSize ) it = (*lists)[cur_list]->rend()-1; else it = NULL;}  return *this; }
			_reverse_iterator operator --(int){  _iterator ret = *this; --(*this); return ret; }
			dtype & operator *() { return *it; }
			dtype * operator ->() { return &*it; }
			_reverse_iterator & operator =(_reverse_iterator const & other) {it = other.it; return *this;}
			bool operator ==(const _reverse_iterator & other) { return it == other.it;}
			bool operator !=(const _reverse_iterator & other) { return it != other.it;}
			bool operator <(const _reverse_iterator & other) { return it < other.it;}
			bool operator >(const _reverse_iterator & other) { return it > other.it;}
			bool operator <=(const _reverse_iterator & other) { return it <= other.it;}
			bool operator >=(const _reverse_iterator & other) { return it >= other.it;}
		};
		typedef _reverse_iterator<inner_type> reverse_iterator;
		typedef _reverse_iterator<const inner_type> const_reverse_iterator;
		*/
	public:
		TValue & operator [] (IndType key)
		{
			IndType pos = compute_pos(key);
			return lists[pos][key];
		}
		bool is_present(IndType key)
		{
			IndType pos = compute_pos(key);
			return lists[pos].find(key) != lists[pos].end();
		}
		IndType size()
		{
			IndType ret = 0;
			for(IndType i = 0; i < HSize; i++) ret += lists[i].size();
			return ret;
		}
		std::vector< std::pair<IndType, TValue> > serialize()
		{
			std::vector< std::pair<IndType, TValue> > ret;
			ret.resize(size());
			IndType q = 0;
			for(IndType i = 0; i < HSize; i++)
			{
				inner_iter it;
				for(it = lists[i].begin(); it != lists[i].end(); ++it)
					ret[q++] = std::make_pair(it->first,it->second);
			}
			return ret;
		}
		void clear()
		{
			for(IndType i = 0; i < HSize; i++) lists[i].clear();
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
		element ** chunks;
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
			//~ size_type oldnchunks  = m_size  /chunk;
			size_type oldnchunks2 = (static_cast<uenum>(m_size) >> block_bits) + ( (m_size & block_bits_mask) ? 1 : 0);
			size_type oldn = (static_cast<uenum>(oldnchunks2) >> fwd_alloc_chunk_bits) + ( (oldnchunks2 & fwd_alloc_chunk_bits_mask) ? 1 : 0);
			//~ size_type newnchunks = new_size/chunk;
			size_type newnchunks2 = (static_cast<uenum>(new_size) >> block_bits) + ( (new_size & block_bits_mask)? 1 : 0);
			size_type newn = (static_cast<uenum>(newnchunks2) >> fwd_alloc_chunk_bits) + ( (newnchunks2 & fwd_alloc_chunk_bits_mask) ? 1 : 0);
			
			//~ std::cout << "new " << new_size << " " << newn << " " << newnchunks2 << " old " << m_size << " " << oldn << " " << oldnchunks2 << " block " << block << " chunk " << chunk << std::endl;
			
			for(size_type q = new_size; q < m_size; q++) 
				access_element(q).~element();
			for(size_type q = newnchunks2; q < oldnchunks2; q++) 
			{
				assert(chunks[q] != NULL);
				free(chunks[q]);
				chunks[q] = NULL;
			}
			if( newn != oldn )	
			{
				if( newn > 0 )
				{
					chunks = (element **) realloc(chunks,fwd_alloc_chunk_size*newn);
					assert(chunks != NULL);
					if( newn > oldn ) memset(chunks+oldn*fwd_alloc_chunk_val,0,fwd_alloc_chunk_size*(newn-oldn));
				}
				else
				{
					free(chunks);
					chunks = NULL;
				}
			}
			for(size_type q = oldnchunks2; q < newnchunks2; q++) 
			{
				assert(chunks[q] == NULL);
				chunks[q] = (element *)malloc(block_size);
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
			for(size_type q = 0; q < m_size; q++) access_element(q).~element();
			size_type cend = (static_cast<uenum>(m_size) >> block_bits) + ((m_size & block_bits_mask)? 1 : 0);
			for(size_type q = 0; q < cend; q++) 
			{
				free(chunks[q]);
				chunks[q] = NULL;
			}
			free(chunks);
			chunks = NULL;
			m_size = 0;
		}
		chunk_array()
		{
			m_size = 0;
			chunks = NULL;
		}
		chunk_array(const chunk_array & other)
		{
			chunks = NULL;
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
		/*for future
		chunk_array & operator =(chunk_array && other)
		{
			if( this != &other )
			{
				clear();
				chunks = other.chunks;
				m_size = other.m_size;
				other.chunks = NULL;
				other.m_size = 0;
			}
			return *this;
		}
		*/
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
			//access_element(m_size-1).~element();
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
			//destruct current
			//(*pos).~element();
			iterator it = pos, jt = it++;
			while(it != end()) (*jt++) = (*it++);//std::move(*jt++);
			//destruct last
			//access_element(m_size-1).~element();
			inner_resize(m_size-1);
			m_size--;
			return pos;
		}
		
		iterator begin() {return iterator(this,0);}
		iterator end() {return iterator(this,m_size);}
		const_iterator begin() const {return const_iterator(this,0);}
		const_iterator end() const {return const_iterator(this,m_size);}
		
		//~ reverse_iterator rbegin() {return reverse_iterator(this,m_size-1);}
		//~ reverse_iterator rend() {return reverse_iterator(this,-1);}
		//~ const_reverse_iterator rbegin() const {return const_reverse_iterator(this,m_size-1);}
		//~ const_reverse_iterator rend() const {return const_reverse_iterator(this,-1);}
		
		//don't need other standard functions?
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
		char ** chunks;
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
				free(chunks[q]);
				chunks[q] = NULL;
			}
			if( newn != oldn )	
			{
				if( newn > 0 )
				{
					chunks = reinterpret_cast<char **>(realloc(chunks,fwd_alloc_chunk_size*newn));
					assert(chunks != NULL);
					if( newn > oldn ) memset(chunks+oldn*fwd_alloc_chunk_val,0,fwd_alloc_chunk_size*(newn-oldn));
				}
				else
				{
					free(chunks);
					chunks = NULL;
				}
			}
			for(size_type q = oldnchunks2; q < newnchunks2; q++) 
			{
				assert(chunks[q] == NULL);
				chunks[q] = static_cast<char *>(malloc(block_size*record_size));
				memset(chunks[q],0,block_size*record_size);
				assert(chunks[q] != NULL);
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
				free(chunks[q]);
				chunks[q] = NULL;
			}
			free(chunks);
			chunks = NULL;
			m_size = 0;
		}
		chunk_bulk_array(size_type set_record_size = 1)
		{
			record_size = set_record_size;
			m_size = 0;
			chunks = NULL;
		}
		chunk_bulk_array(const chunk_bulk_array & other)
		{
			chunks = NULL;
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
}

#endif
