#include "memory_pool.h"
#include <stdio.h>

MemoryPool pool;
bool is_initialized = false;

void * get_memory(size_t bytes, bool in_permanent)
{
    if (!is_initialized) {
        pool.initialize(gigabytes(4), gigabytes(1));
        is_initialized = true;
    }
    if (in_permanent) {
        return pool.get_from_permanent(bytes);
    } else {
        return pool.get_from_transient(bytes);
    }
}

void reset_memory()
{
    pool.reset_transient();
}

void memory_diagnostic()
{
    printf("permanent = %zu / %zu (%lf), transient = %zu / %zu (%lf)\n",
            pool.permanent_remaining, pool.permanent_size, 
            ((double) 100.0 * pool.permanent_remaining) / pool.permanent_size,
            pool.transient_remaining, pool.transient_size, 
            ((double) 100.0 * pool.transient_remaining) / pool.transient_size);
}
