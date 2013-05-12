#include "gspan_allocator.hpp"
#include <new>
#include <cassert>

#ifndef BR
#define BR asm volatile ("int3;")
#endif


namespace gSpan
{

    // -------------- FixedAllocator -----------------------
    
    FixedAllocator::~FixedAllocator()
    {
	while (!free_ptrs_.empty())
	{
	    void* p = free_ptrs_.top();
	    free_ptrs_.pop();
	    ::operator delete(p);
	}
    }

    
    void* FixedAllocator::allocate()
    {
	if (!free_ptrs_.empty())
	{
	    void* p = free_ptrs_.top();
	    free_ptrs_.pop();
	    return p;
	}
	return ::operator new (data_size_);
    }

    // ------------- MemAllocator --------------------------

    MemAllocator::~MemAllocator()
    {
	for (std::size_t i = 0; i < fallocs_.size(); ++i)
	    delete fallocs_[i];
    }

    FixedAllocator* MemAllocator::get_fixed_allocator(std::size_t data_size)
    {
	assert(data_size > 0);
	std::size_t old_size = fallocs_.size();
	if (old_size < data_size)
	{
	    fallocs_.resize(data_size);
	    for (std::size_t i = old_size; i < data_size; ++i)
		fallocs_[i] = new FixedAllocator(i + 1);
	}
	return fallocs_[data_size - 1];
    }

    void MemAllocator::deallocate(void* p, std::size_t data_size)
    {
	assert(data_size <= fallocs_.size());
	fallocs_[data_size - 1]->deallocate(p);
    }
}
