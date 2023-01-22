// A basic hash table (using separate chaining)

#ifndef HASH_H
#define HASH_H

#include <stdint.h>

typedef uint32_t hsize_t;
typedef uint32_t hkey_t;
typedef void* hval_t;

/** An element of a hash table */
typedef struct {
    hsize_t key;
    hval_t value;
    hsize_t next;
} HashElement;

/** A hash table */
typedef struct {
    hsize_t size;
    hsize_t capacity;
    HashElement *data;
    hsize_t *table;
} Hash;

/** Create a hash table */
Hash *hash_create();

/** Destroy a hash table */
void hash_free(Hash *hash);

/** Get a value from a hash table if exists */
hval_t hash_get(const Hash *hash, const hkey_t key);

/** Get an element of a hash table, create if does not exist */
HashElement *hash_entry(Hash *hash, const hkey_t key);

#endif
