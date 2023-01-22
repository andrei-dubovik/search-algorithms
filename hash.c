// A basic hash table (using separate chaining)

#include <stdlib.h>
#include <string.h>
#include "hash.h"


/** See https://github.com/skeeto/hash-prospector */
static hsize_t hash32(hkey_t x) {
    x ^= x >> 16;
    x *= 0x21f0aaad;
    x ^= x >> 15;
    x *= 0xd35a2d97;
    x ^= x >> 15;
    return x;
}


/** Create a hash table */
Hash *hash_create() {
    Hash *hash = malloc(sizeof(Hash));
    hash->size = 0;
    hash->capacity = 1;
    hash->data  = calloc(1, sizeof(HashElement));
    hash->table = calloc(1, sizeof(hsize_t));
    return hash;
}


/** Destroy a hash table */
void hash_free(Hash *hash) {
    free(hash->data);
    free(hash->table);
    free(hash);
}


/** Reallocate memory, clear extra memory */
static void *recallocarray(void *ptr, size_t oldnmemb, size_t nmemb, size_t size) {
    ptr = reallocarray(ptr, nmemb, size);
    if(nmemb > oldnmemb) {
        memset(ptr + oldnmemb*size, 0, (nmemb - oldnmemb)*size);
    }
    return ptr;
}


/** Double hash capacity */
void hash_expand(Hash *hash) {
    // Expand capacity
    hsize_t oldcap = hash->capacity;
    hash->capacity *= 2;
    hash->data = recallocarray(
        hash->data,
        oldcap,
        hash->capacity,
        sizeof(HashElement)
    );

    // Reset the hash table
    free(hash->table);
    hash->table = calloc(hash->capacity, sizeof(hsize_t));

    // Rebuild the hash table
    for(hsize_t j = 0; j < hash->size; j++) {
        hsize_t i = hash32(hash->data[j].key) % hash->capacity;
        hash->data[j].next = hash->table[i];
        hash->table[i] = j + 1;
    }
}


/** Get a value from a hash table if exists */
hval_t hash_get(const Hash *hash, const hkey_t key) {
    hsize_t i = hash32(key) % hash->capacity;
    hsize_t j = hash->table[i];
    while(j != 0) {
        j -= 1;
        if(hash->data[j].key == key) {
            return hash->data[j].value;
        }
        j = hash->data[j].next;
    }
    return NULL;
}


/** Get an element of a hash table, create if does not exist */
HashElement *hash_entry(Hash *hash, const hkey_t key) {
    // Search for existing record
    hsize_t i = hash32(key) % hash->capacity;
    hsize_t j0 = hash->table[i];
    hsize_t j = j0;
    while(j != 0) {
        j -= 1;
        if(hash->data[j].key == key) {
            return &hash->data[j];
        }
        j = hash->data[j].next;
    }

    // Expand the hash if full capacity is reached
    if(hash->size == hash->capacity) {
        hash_expand(hash);
        i = hash32(key) % hash->capacity;
        j0 = hash->table[i];
    }

    // Create a new record
    // LIFO preference
    j = hash->size++;
    hash->data[j].key = key;
    hash->data[j].next = j0;
    hash->table[i] = j + 1;
    return &hash->data[j];
}
