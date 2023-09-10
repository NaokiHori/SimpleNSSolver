#if !defined(MEMORY_H)
#define MEMORY_H

#include <stddef.h> // size_t

// general-purpose memory allocator
extern void * memory_calloc(
    const size_t count,
    const size_t size
);

// corresponding memory deallocator
extern void memory_free(
    void * ptr
);

#endif // MEMORY_H
