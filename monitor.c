// A monitor to keep track of dynamic memory usage

#include <malloc.h>
#include "monitor.h"

size_t memory = 0;
size_t peak_memory = 0;

void *__real_malloc(size_t size);
void *__real_calloc(size_t nmemb, size_t size);
void *__real_realloc(void *ptr, size_t size);
void *__real_reallocarray(void *ptr, size_t nmemb, size_t size);
void __real_free(void *ptr);


/** Register the amount of newly allocated memory */
static void add_memory(void *ptr) {
    memory += malloc_usable_size(ptr);
    if(memory > peak_memory) {
        peak_memory = memory;
    }
}


/** Deregister the amount of memory to be freed */
static void sub_memory(void *ptr) {
    memory -= malloc_usable_size(ptr);
}


void *__wrap_malloc(size_t size) {
    void *ptr = __real_malloc(size);
    add_memory(ptr);
    return ptr;
}


void *__wrap_calloc(size_t nmemb, size_t size) {
    void *ptr = __real_calloc(nmemb, size);
    add_memory(ptr);
    return ptr;
}


void *__wrap_realloc(void *ptr, size_t size) {
    sub_memory(ptr);
    ptr = __real_realloc(ptr, size);
    add_memory(ptr);
    return ptr;
}


void *__wrap_reallocarray(void *ptr, size_t nmemb, size_t size) {
    sub_memory(ptr);
    ptr = __real_reallocarray(ptr, nmemb, size);
    add_memory(ptr);
    return ptr;
}


void __wrap_free(void *ptr) {
    sub_memory(ptr);
    __real_free(ptr);
}
