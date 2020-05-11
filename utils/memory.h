#ifndef CADO_UTILS_MEMORY_H_
#define CADO_UTILS_MEMORY_H_

// IWYU pragma: private, include "utils.h"

#include <stdlib.h>
#ifdef __cplusplus
#include <memory>
#endif

#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void * malloc_check(const size_t x);
ATTRIBUTE((malloc)) extern void * malloc_aligned(size_t size, size_t alignment);
ATTRIBUTE((warn_unused_result)) void * realloc_aligned(void * p, 
        const size_t old_size, const size_t new_size, const size_t alignment);
extern void free_aligned(void * ptr);

extern void * aligned_alloc (size_t alignment, size_t size);

extern void * malloc_pagealigned(size_t sz) ATTR_ASSUME_ALIGNED(32);
extern void free_pagealigned(void * ptr);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
template<typename T, int align = sizeof(T)> class aligned_allocator : public std::allocator<T> {
    typedef std::allocator<T> super;
    public:
    template <typename U> struct rebind {
        typedef aligned_allocator<U, align> other;
    } ;
    typename super::pointer allocate(size_t n) const {
        return (typename super::pointer) malloc_aligned(n * sizeof(T), align);
    }
    void deallocate(typename super::pointer p, size_t) const {
        return free_aligned(p);
    }
    template <typename X>
        T * allocate(const size_t n, const X *) const {
            return allocate(n);
        }
};
#endif

#endif
