#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"
#include "hash.h"
#include "monitor.h"


/** Maximum number of MCTS iterations */
const int MCTS_STEPS = 100000;

/** Search node */
typedef struct {
    int i;     // state (depth)
    int j;     // state (offset)
    double V;  // cumulative value
    int N;     // number of passes
} Node;

#ifdef TRACK_VISITS
int *visit_stats;
#endif


/** Generate a normally distributed random variable (Box-Muller) */
double rnorm(double mu, double sigma) {
    double u = (double)(random() + 1)/0x80000001;
    double v = (double)(random() + 1)/0x80000001;
    double x = sqrt(-2.0*log(u))*cos(2.0*M_PI*v);
    return x*sigma + mu;
}


/** Generate a uniformly distributed random variable */
double runif() {
    return (double)random()/0x7fffffff;
}


/** Initialize a new search node */
Node *make_node(const int i, const int j) {
    Node *node = malloc(sizeof(Node));
    node->i = i;
    node->j = j;
    node->V = 0.0;
    node->N = 0;
    return node;
}


/** Sample a random delta u using a truncated Wiener process */
int rdeltau(int u, double sigma) {
    double s = sigma*sqrt(4.0/T);
    while(1) {
        int du = round(rnorm(0.0, s)*M);
        if(u + du >= 0 && u + du <= N) return du;
    }
}


/** Compute the value of a random path from a given point */
double sample(const int i0, const int j0, double sigma) {
    double V = 0.0;
    int j = j0;
    int uj = M;
    for(int i = i0; i <= T; i++) {
        #ifdef TRACK_VISITS
        if(i > i0) visit_stats[idx(i, j)]++;
        #endif
        if(i < T) {
            uj += rdeltau(uj, sigma);
        } else {
            uj = M;  // Fixed action at the terminal point
        }
        V += -objf(i, j, uj);
        j += uj;
    }
    return V;
}


/** Select a child node that maximizes the value function */
Node *argmaxv(Node *node, Hash *visited) {
    double maxv = -INFINITY;
    Node *arg = NULL;
    for(int uj = 0; uj <= N; uj++) {
        int j = node->j + uj;
        int s = idx(node->i + 1, j);
        #ifdef TRACK_VISITS
        visit_stats[s]++;
        #endif
        Node *child = (Node*)hash_get(visited, s);
        if(child != NULL) {
            double v = -objf(node->i, node->j, uj) + child->V;
            if(v > maxv) {
                maxv = v;
                arg = child;
            }
        }
    }
    return arg;
}


/** Select the next node using an e-greedy policy */
Node *traverse(Node *node, Hash *visited, int parent_u, double eps, double sigma) {
    if(runif() > eps/log(node->N)) {  // eps/0.0 = +inf
        // Select a visited node greedily
        Node *child = argmaxv(node, visited);
        if(child != NULL) {
            return child;
        }
    }

    // Select a node at random
    int uj;
    if(node->i == 0) {
        uj = random()%(N + 1);
    } else {
        uj = parent_u + rdeltau(parent_u, sigma);
    }
    int j = node->j + uj;
    int s = idx(node->i + 1, j);
    #ifdef TRACK_VISITS
    visit_stats[s]++;
    #endif
    Node **child = (Node**)&hash_entry(visited, s)->value;
    if(*child == NULL) {
        *child = make_node(node->i + 1, j);
    }
    return *child;
}


/** Grow a search tree */
void grow(Node *node, Hash *visited, double eps, double sigma) {
    Node *path[T];
    int t = 0;
    int uj = M;


    // Selection, expansion
    // (stop at the boundary or if it's a newly added node)
    while(node->i < T && node->N > 0) {
        path[t++] = node;
        node = traverse(node, visited, uj, eps, sigma);
        uj = node->j - path[t-1]->j;
    }

    // Simulation for newly created nodes
    if(node->N++ == 0) {
        node->V = sample(node->i, node->j, sigma);
    }


    // Back-propogation
    while(t-- > 0) {
        Node *parent = path[t];
        uj = node->j - parent->j;
        double V = -objf(parent->i, parent->j, uj) + node->V;
        if(V > parent->V) parent->V = V;
        parent->N++;
        node = parent;
        #ifdef TRACK_VISITS
        visit_stats[idx(node->i, node->j)]++;
        #endif
    }
}


/** Extract the optimal path from the solution */
int *optimal_path(Node *node, Hash *visited) {
    int *path = calloc(T + 1, sizeof(int));
    for(int i = 0; i < T; i++) {
        Node *child = argmaxv(node, visited);
        if(child == NULL) {
            fprintf(stderr, "Error: the optimal path is incomplete. Aborting.\n");
            exit(1);
        }
        path[i] = child->j - node->j;
        node = child;
    }
    path[T] = M;
    return path;
}


/** Solve the problem using Monte Carlo tree search */
void mcts(double eps, double sigma, double *value, int **path) {
    // Initialize RNG (initial seed has been sampled from /dev/urandom)
    srandom(0xb64a8ac3);

    // Create an empty search tree
    Node *tree = make_node(0, 0);

    // Initialize a hash for visited nodes (it's a recombining tree)
    Hash *visited = hash_create();

    #ifdef TRACK_VISITS
    if(path != NULL) visit_stats = calloc(tree_size, sizeof(int));
    #endif

    // Run Monte Carlo Tree Search
    for(int i = 0; i < MCTS_STEPS; i++) {
        grow(tree, visited, eps, sigma);
    }
    *value = -tree->V;

    if(path != NULL) {
        // Report various statistics
        printf("Stored states: %'i (%.1f%%)\n", visited->size, 100.0*visited->size/tree_size);
        printf("Memory / state: %'i bytes\n", peak_memory/visited->size);
        #ifdef TRACK_VISITS
        print_visits(visit_stats);
        #endif

        // Save the solution
        *path = optimal_path(tree, visited);
    }

    // Free memory
    for(hsize_t i = 0; i < visited->size; i++) {
        free(visited->data[i].value);
    }
    hash_free(visited);
    free(tree);
    #ifdef TRACK_VISITS
    if(path != NULL) free(visit_stats);
    #endif
}


/** Minimize a function using golden search (verbose) */
double golden_section(double (*f)(double, void *), double x1, double x4, void *args, int steps) {
    const double invphi = (sqrt(5.0) - 1.0)/2.0;
    double x2 = x4 - invphi*(x4 - x1);
    double x3 = x1 + invphi*(x4 - x1);

    double f1 = f(x1, args);
    double f2 = f(x2, args);
    double f3 = f(x3, args);
    double f4 = f(x4, args);

    printf("Initializing:\n");
    printf("x1 = %.3f, f = %.3f\n", x1, f1);
    printf("x2 = %.3f, f = %.3f\n", x2, f2);
    printf("x3 = %.3f, f = %.3f\n", x3, f3);
    printf("x4 = %.3f, f = %.3f\n", x4, f4);
    printf("\n");

    for(int i = 0; i < steps; i++) {
        // f1 < f2,f3,f4 is an heuristic for non-convex or stochastic functions
        if(f2 < f3 || (f1 < f2 && f1 < f3 && f1 < f4)) {
            x4 = x3;
            f4 = f3;
            x3 = x2;
            f3 = f2;
            x2 = x4 - invphi*(x4 - x1);
            f2 = f(x2, args);
        } else {
            x1 = x2;
            f1 = f2;
            x2 = x3;
            f2 = f3;
            x3 = x1 + invphi*(x4 - x1);
            f3 = f(x3, args);
        }
        printf("Step %i:\n", i + 1);
        printf("x1 = %.3f, f = %.3f\n", x1, f1);
        printf("x2 = %.3f, f = %.3f\n", x2, f2);
        printf("x3 = %.3f, f = %.3f\n", x3, f3);
        printf("x4 = %.3f, f = %.3f\n", x4, f4);
        printf("\n");
    }

    return f2 < f3 ? x2 : x3;
}


// Lambda functions

double mcts_eps(double x, void *args) {
    double V;
    mcts(x, *(double*)args, &V, NULL);
    return V;
}


double mcts_sigma(double x, void *args) {
    double V;
    mcts(*(double*)args, x, &V, NULL);
    return V;
}


/** Optimize MCTS over metaparameters (verbose) */
void optimize(double *eps, double *sigma) {
    const int steps1 = 4;
    const int steps2 = 10;
    double eps_ub = 2.0*(*eps);
    double sigma_ub = 2.0*(*sigma);
    printf("\e[31mStarting at eps = %.3f, sigma = %.3f\e[0m\n\n", *eps, *sigma);
    for(int i = 0; i < steps1; i++) {
        printf("Optimizing in eps...\n\n");
        *eps = golden_section(mcts_eps, 0.0, eps_ub, sigma, steps2);
        printf("Optimizing in sigma...\n\n");
        *sigma = golden_section(mcts_sigma, 0.0, sigma_ub, eps, steps2);
        printf("\e[31mStep %i: eps = %.3f, sigma = %.3f\e[0m\n\n", i + 1, *eps, *sigma);
    }
}


int main() {
    setlocale(LC_NUMERIC, "");

    // Optimize over metaparameters
    // double eps = 2.0;
    // double sigma = 1.0;
    // optimize(&eps, &sigma);

    // We start with optimized values
    double eps = 2.885;
    double sigma = 0.318;

    // Solve the problem
    double V, timer;
    int *path;
    TIMEIT(mcts(eps, sigma, &V, &path), timer);

    // Save the solution
    dump(path, "mcts.dat");
    free(path);

    // Report the optimal objective function
    printf("V = %f\n", V);

    // Report resource usage
    printf("Running time: %.2f s\n", timer);
    printf("Peak dynamic memory usage: %'i bytes\n", peak_memory);

    return 0;
}
