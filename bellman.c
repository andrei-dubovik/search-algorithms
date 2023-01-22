#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"
#include "monitor.h"


/** Compute tree breadth for given depth */
int breadth(const int depth) {
    return 2*M*depth + 1;
}


/** Extract the optimal path from the solution */
int *optimal_path(const int *u) {
    int *path = calloc(T + 1, sizeof(int));
    int j = 0;
    for(int i = 0; i <= T; i++) {
        int uj = u[idx(i, j)];
        path[i] = uj;
        j += uj;
    }
    return path;
}


/** Solve the problem using the Bellman's optimality principle */
void bellman(double *value, int **path) {
    // Declare and allocate data structures
    int size_V = (2*T*M + 1)*sizeof(double);
    int size_u = tree_size*sizeof(int);
    double *V = malloc(size_V);  // Value function (one tree level)
    int *u = malloc(size_u);  // Control function (whole tree)
    #ifdef TRACK_VISITS
    int *visit_stats = calloc(tree_size, sizeof(int));
    #endif

    // Initialize data structures
    for (int i = 0; i < breadth(T); i++) {
        V[i] = objf(T, i, M);
        u[idx(T, i)] = M;
    }

    // Solve the problem backwards
    for (int i = T - 1; i >= 0; i--) {
        for (int j = 0; j < breadth(i); j++) {
            #ifdef TRACK_VISITS
            visit_stats[idx(i, j)]++;
            #endif
            int uj = M;
            double min_vj = objf(i, j, uj) + V[j + uj];
            for (int k = 0; k < N + 1; k++) {
                #ifdef TRACK_VISITS
                visit_stats[idx(i + 1, j + k)]++;
                #endif
                double vj = objf(i, j, k) + V[j + k];
                if (vj < min_vj) {
                    min_vj = vj;
                    uj = k;
                }
            }
            V[j] = min_vj;
            u[idx(i, j)] = uj;
        }
    }

    // Save the results
    *value = V[0];
    *path = optimal_path(u);

    // Report various statistics
    #ifdef TRACK_VISITS
    print_visits(visit_stats);
    #endif

    // Free memory
    free(u);
    free(V);
    #ifdef TRACK_VISITS
    free(visit_stats);
    #endif
}


int main() {
    setlocale(LC_NUMERIC, "");

    // Solve the problem
    double V, timer;
    int *path;
    TIMEIT(bellman(&V, &path), timer);

    // Save the solution
    dump(path, "bellman.dat");
    free(path);

    // Report the optimal objective function
    printf("V = %f\n", V);

    // Report resource usage
    printf("Running time: %.2f s\n", timer);
    printf("Peak dynamic memory usage: %'i bytes\n", peak_memory);
    printf("Memory / state: %'i bytes\n", peak_memory/tree_size);

    return 0;
}
