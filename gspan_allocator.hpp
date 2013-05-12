#ifndef GSPAN_ALLOCATOR_H_
#define GSPAN_ALLOCATOR_H_

#include <cstdlib>
#include <vector>
#include <stack>
#include <memory>
#include <cassert>

#include <stdint.h>

namespace gSpan
{
    class FixedAllocator
    {
	std::stack<void*> free_ptrs_;
	std::size_t data_size_;
    public:
	FixedAllocator(std::size_t data_size) :data_size_(data_size) {}
	~FixedAllocator();

	void* allocate();
	void deallocate(void* p)	{ free_ptrs_.push(p); }
    };


    class MemAllocator
    {
	std::vector<FixedAllocator*> fallocs_;
    public:
	~MemAllocator();

	FixedAllocator* get_fixed_allocator(std::size_t data_size);

	void* allocate(std::size_t data_size)
	    { return get_fixed_allocator(data_size)->allocate(); }

	void deallocate(void*, std::size_t data_size);

	template<class T>
	T* alloc_array(std::size_t n)
	    {
		return static_cast<T*>(allocate(sizeof(T) * n));
	    }

	template<class T>
	void dealloc_array(T* p, std::size_t n)
	    {
		if (p)
		{
		    assert(n > 0);
		    deallocate(p, sizeof(T) * n);
		}
	    }
    };


    template<class T>
    class STL_Allocator : public std::allocator<T>
    {
    public:
	typedef T		value_type;
	typedef T*		pointer;
	typedef T&		reference;
	typedef const T*	const_pointer;
	typedef const T&	const_reference;
	typedef std::size_t	size_type;
	typedef std::ptrdiff_t	difference_type;

	template <class U> 
	struct rebind { typedef STL_Allocator<U> other; };

	explicit
	STL_Allocator(MemAllocator* m) throw() :malloc_(m) {}

	STL_Allocator(const STL_Allocator& r) throw() :malloc_(r.malloc_) {}

	template <class U>
	STL_Allocator(const STL_Allocator<U>& r) throw() :malloc_(r.get_mem_alloc()) {}

	pointer allocate(size_type n, std::allocator<void>::const_pointer hint = 0)
	    { return static_cast<pointer>(alloc_(n)); }

	void deallocate(pointer p, size_type n)
	    { dealloc_(p, n); }

	MemAllocator* get_mem_alloc() const { return malloc_; }
    private:
	MemAllocator* malloc_;	// no owner
	
	const static std::size_t OBJECT_ALIGN = sizeof(T);
	typedef size_type BlockHeader;

	static std::size_t alignup(std::size_t x)
	    {	
		std::size_t a = OBJECT_ALIGN - 1U;
		return (x + a) & ~a;
	    }

	static void* alignup_ptr(void* p)
	    {
		std::size_t a = OBJECT_ALIGN - 1U;
		::uintptr_t x = reinterpret_cast< ::uintptr_t >(p);
		return reinterpret_cast<void*>((x + a) & ~a);
	    }

	void* alloc_(size_type n);
	void dealloc_(void* p, size_type n);
    };

    template<class T>
    void* STL_Allocator<T>::alloc_(size_type n)
    {
	std::size_t data_offset = alignup(sizeof(BlockHeader));
	std::size_t alloc_size = alignup(data_offset + alignup(sizeof(T) * n));
	unsigned char* block_ptr = reinterpret_cast<unsigned char*>(malloc_->allocate(alloc_size));
	return reinterpret_cast<void*>(block_ptr + data_offset);
    }

    template<class T>
    void STL_Allocator<T>::dealloc_(void* p, size_type n)
    {
	std::size_t data_offset = alignup(sizeof(BlockHeader));
	std::size_t alloc_size = alignup(data_offset + alignup(sizeof(T) * n));
	::uintptr_t x = reinterpret_cast< ::uintptr_t >(p);
	void* block_ptr = reinterpret_cast<void*>(x - data_offset);
	malloc_->deallocate(block_ptr, alloc_size);
    }

}

#endif
