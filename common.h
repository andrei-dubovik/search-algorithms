#ifndef COMMON_H
#define COMMON_H

extern const int M;
extern const int N;
extern const int T;
extern const int tree_size;

/** Compute the objective function at a given point in the tree */
double objf(const int i, const int j, const int uj);

/** Compute a flat index for a tree node */
int idx(int i, int j);

/** Traverse the optimal path */
void dump(const int *path, const char *name);

#ifdef TRACK_VISITS
void print_visits(int *stats);
#endif

// Macros

#define TIMEIT(expr, timer) \
    clock_t _clock = clock(); \
    expr; \
    _clock = clock() - _clock; \
    timer = (double)_clock/CLOCKS_PER_SEC


#endif
