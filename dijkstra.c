#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"
#include "hash.h"
#include "heap.h"
#include "monitor.h"


/** Search node */
struct Node {
    int i;                // state (depth)
    int j;                // state (offset)
    struct Node *parent;
    double V;             // accumulated value
    int rc;               // number of strong references
    int position;         // priority queue position
};

typedef struct Node Node;


/** Get a reference to node's position */
int *node_position(void *node) {
    return &((Node*)node)->position;
}


/** Increase reference count */
Node *borrow(Node *node) {
    if(node == NULL) return NULL;
    node->rc++;
    return node;
}


/** Decrease reference count, free memory */
void release(Node *node) {
    if(node == NULL) return;
    if(--node->rc == 0) {
        Node *parent = node->parent;
        free(node);
        release(parent);
    }
}


/** Initialize a new node */
Node *make_node(int i, int j, Node *parent, double V) {
    Node *node = malloc(sizeof(Node));
    node->i = i;
    node->j = j;
    node->parent = borrow(parent);
    node->V = V;
    node->rc = 1;
    return node;
}


/** Compute the objective function at a given point in the tree */
double fplus(const int i, const int j, const int uj) {
    const double dt = 4.0/T;
    const double du = 1.0/M;
    double u = uj*du - 1.0;
    double x = -i*dt + j*dt*du;
    return (u*u + x + 4.0)*dt;  // x + 4.0 >= 0 (can't have negative costs)
}


/** Extract the optimal path from the solution */
int *optimal_path(Node *node) {
    int *path = calloc(T + 1, sizeof(int));
    path[T] = M;
    for(int i = T - 1; i >= 0; i--) {
        path[i] = node->j - node->parent->j;
        node = node->parent;
    }
    return path;
}


/** Solve the problem using the Dijkstra's algorithm */
void dijkstra(double *value, int **path) {
    // Initialize data structures
    Node *node = make_node(0, 0, NULL, 0.0);
    Heap *frontier = heap_create(node_position);
    heap_insert(frontier, borrow(node), 0.0);
    Hash *reached = hash_create();
    hash_entry(reached, 0)->value = borrow(node);
    release(node);
    #ifdef TRACK_VISITS
    int *visit_stats = calloc(tree_size, sizeof(int));
    #endif

    // Best-first search (Dijkstra's algorithm)
    while(frontier->size > 0) {
        node = heap_pop(frontier);
        #ifdef TRACK_VISITS
        visit_stats[idx(node->i, node->j)]++;
        #endif
        if(node->i == T) break;  // goal state
        for(int k = 0; k < N + 1; k++) {
            Node *child = make_node(
                node->i + 1,
                node->j + k,
                node,
                node->V + fplus(node->i + 1, node->j + k, k)
            );
            int s = idx(child->i, child->j);
            #ifdef TRACK_VISITS
            visit_stats[s]++;
            #endif
            Node **reached_s = (Node**)&hash_entry(reached, s)->value;
            if(*reached_s == NULL) {
                *reached_s = borrow(child);
                heap_insert(frontier, borrow(child), -child->V);
            } else if(child->V < (*reached_s)->V) {
                heap_replace(frontier, *reached_s, borrow(child), -child->V);
                release(*reached_s);  // `frontier` reference
                release(*reached_s);  // `reached` reference
                *reached_s = borrow(child);
            }
            release(child);
        }
        release(node);
    }

    // Save the solution
    *path = optimal_path(node);
    *value = node->V - 16.0;

    // Report various statistics
    printf("Frontier size: %'i\n", frontier->size);
    printf("Stored states: %'i (%.1f%%)\n", reached->size, 100.0*reached->size/tree_size);
    printf("Memory / state: %'i bytes\n", peak_memory/reached->size);
    #ifdef TRACK_VISITS
    print_visits(visit_stats);
    #endif

    // Free memory
    release(node);
    for(hsize_t i = 0; i < reached->size; i++) {
        release(reached->data[i].value);
    }
    hash_free(reached);
    while(node = heap_pop(frontier)) {
        release(node);
    }
    heap_free(frontier);
    #ifdef TRACK_VISITS
    free(visit_stats);
    #endif
}


int main() {
    setlocale(LC_NUMERIC, "");

    // Solve the problem
    double V, timer;
    int *path;
    TIMEIT(dijkstra(&V, &path), timer);

    // Save the solution
    dump(path, "dijkstra.dat");
    free(path);

    // Report the optimal objective function
    printf("V = %f\n", V);

    // Report resource usage
    printf("Running time: %.2f s\n", timer);
    printf("Peak dynamic memory usage: %'i bytes\n", peak_memory);

    return 0;
}
