// A basic max-heap (priority queue)

#ifndef HEAP_H
#define HEAP_H

/** A an element of a priority queue */
typedef struct {
    double priority;
    void *value;
} Element;

/** A priority queue */
typedef struct {
    int size;
    int capacity;
    Element *tree;
    int *(*position)(void*);
} Heap;

/** Initialize a priorioty queue */
Heap *heap_create(int *(*position)(void*));

/** Add an element to a priority queue */
void heap_insert(Heap *heap, void *value, const double priority);

/** Pop the element with the highest priority */
void *heap_pop(Heap *heap);

/** Replace a given element with a new element of different priority */
void heap_replace(Heap *heap, void *oldval, void *newval, const double priority);

/** Destroy a priority queue */
void heap_free(Heap *heap);

#endif
