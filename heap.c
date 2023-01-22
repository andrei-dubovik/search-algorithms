// A basic max-heap (priority queue)
#include <stdlib.h>
#include "heap.h"


/** Initialize a priorioty queue */
Heap *heap_create(int *(*position)(void*)) {
    Heap *heap = malloc(sizeof(Heap));
    heap->size = 0;
    heap->capacity = 1;
    heap->tree = malloc(sizeof(Element));
    heap->position = position;
    return heap;
}


/** Sift up */
static void sift_up(Heap *heap, int h) {
    Element *tree = heap->tree;
    int *(*pos)(void*) = heap->position;

    int p;
    while(h > 0 && tree[p = (h-1)/2].priority < tree[h].priority) {
        Element buf = tree[p];
        tree[p] = tree[h];
        tree[h] = buf;
        *pos(tree[p].value) = p;
        *pos(tree[h].value) = h;
        h = p;
    }
}


/** Sift down */
static void sift_down(Heap *heap, int h) {
    Element *tree = heap->tree;
    int *(*pos)(void*) = heap->position;

    while(1) {
        int c1 = 2*h + 1;
        int c2 = 2*h + 2;
        c1 = c1 < heap->size ? c1 : h;
        c2 = c2 < heap->size ? c2 : h;
        int c = tree[c1].priority > tree[c2].priority ? c1 : c2;
        if(tree[c].priority <= tree[h].priority) break;
        Element buf = tree[c];
        tree[c] = tree[h];
        tree[h] = buf;
        *pos(tree[c].value) = c;
        *pos(tree[h].value) = h;
        h = c;
    }
}


/** Add an element to a priority queue */
void heap_insert(Heap *heap, void *value, const double priority) {
    // Resize if necessary
    if(++heap->size > heap->capacity) {
        heap->capacity *= 2;
        heap->tree = realloc(heap->tree, heap->capacity*sizeof(Element));
    }

    // Add new element
    Element *tree = heap->tree;
    int h = heap->size - 1;
    tree[h].priority = priority;
    tree[h].value = value;
    *heap->position(value) = h;

    // Sift up
    sift_up(heap, h);
}


/** Pop the element with the highest priority */
void *heap_pop(Heap *heap) {
    if(heap->size == 0) return NULL;
    Element *tree = heap->tree;
    void *value = tree[0].value;

    // Remove first element
    tree[0] = tree[--heap->size];

    // Sift down
    sift_down(heap, 0);

    return (void*)value;
}


/** Replace a given element with a new element of different priority */
void heap_replace(Heap *heap, void *oldval, void *newval, const double priority) {
    Element *tree = heap->tree;
    int *(*pos)(void*) = heap->position;

    // Change the value and priority
    int h = *pos(oldval);
    double oldpriority = tree[h].priority;
    tree[h].value = newval;
    tree[h].priority = priority;
    *pos(newval) = h;

    // Sift up or sift down
    if(priority > oldpriority) {
        sift_up(heap, h);
    } else {
        sift_down(heap, h);
    }
}


/** Destroy a priority queue */
void heap_free(Heap *heap) {
    free(heap->tree);
    free(heap);
}
