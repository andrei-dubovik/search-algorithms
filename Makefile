## Variables

# Compile these algorithms
algos := bellman dijkstra mcts pontryagin

# Use wrappers around malloc, etc. to track peak memory usage
lmon := -Wl,--wrap=malloc \
        -Wl,--wrap=calloc \
        -Wl,--wrap=realloc \
        -Wl,--wrap=reallocarray \
        -Wl,--wrap=free

# Track visits when compiling for profiling (uses additional memory)
profile := -DTRACK_VISITS

# Common dependencies
deps = monitor.h monitor.c common.h common.c | %

# Compilation command
gcc = gcc -O2 -lm $($*) $(lmon) $(filter %.c,$^) -o $@

## Rules

all: $(algos:%=release/%) $(algos:%=profile/%)

release profile :
	mkdir $@

%/bellman: bellman.c $(deps)
	$(gcc)

%/dijkstra: dijkstra.c heap.h heap.c hash.h hash.c $(deps)
	$(gcc)

%/mcts: mcts.c hash.c hash.h $(deps)
	$(gcc)

%/pontryagin: pontryagin.c $(deps)
	$(gcc)

clean:
	rm -rf release profile
