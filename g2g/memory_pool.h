#pragma once
#ifndef __MEMORY_POOL_H
#define __MEMORY_POOL_H

#include <omp.h>

static const size_t ALIGN = 64;

#define gigabytes(amount) (((size_t) (amount)) * 1024 * 1024 * 1024)

struct MemoryPool {
    MemoryPool() { }
    void initialize(size_t permanent_bytes, size_t transient_bytes)
    {
        size_t total_bytes = permanent_bytes + transient_bytes;

        char * memory = NULL;
        posix_memalign((void **) &memory, ALIGN, total_bytes);

        permanent_pool = memory;
        transient_pool = memory + permanent_bytes;

        permanent_size = permanent_remaining = permanent_bytes;
        transient_size = transient_remaining = transient_bytes;

        omp_init_lock(&permanent_lock);
        omp_init_lock(&transient_lock);
    }
    ~MemoryPool()
    {
        free(permanent_pool);
        omp_destroy_lock(&permanent_lock); 
        omp_destroy_lock(&transient_lock); 
    }
    void * get_from_permanent(size_t bytes)
    {
        bytes = (((bytes + ALIGN - 1) / ALIGN) * ALIGN);
        omp_set_lock(&permanent_lock);
        void * result = NULL;
        if (bytes > permanent_remaining) {
            goto end;
        }
        result = permanent_pool + (permanent_size - permanent_remaining);
        permanent_remaining -= bytes;
end:
        omp_unset_lock(&permanent_lock);
        return result;
    }

    void * get_from_transient(size_t bytes)
    {
        bytes = (((bytes + ALIGN - 1) / ALIGN) * ALIGN);
        omp_set_lock(&transient_lock);
        void * result = NULL;
        if (bytes > transient_remaining) {
            goto end;
        }
        result = transient_pool + (transient_size - transient_remaining);
        transient_remaining -= bytes;
end:
        omp_unset_lock(&transient_lock);
        return result;
    }

    void reset_transient()
    {
        omp_set_lock(&transient_lock);
        transient_remaining = transient_size;
        omp_unset_lock(&transient_lock);
    }

    char * permanent_pool, * transient_pool;
    size_t permanent_size, transient_size;
    size_t permanent_remaining, transient_remaining;
    omp_lock_t permanent_lock, transient_lock;
};

#endif
