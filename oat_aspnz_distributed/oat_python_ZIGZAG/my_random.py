import oat_python as oat
import numpy as np
import random
import time




def random_sampled_subsets(max_vertex=0, num_edges=0, min_edge_size=0, max_edge_size=0, random_seed=0):
    subsets = []
    random.seed(random_seed)
    for _ in range(num_edges):
        # Randomly sample a subset size between 1 and max_vertex (inclusive)
        max_edge_size = min(max_vertex + 1, max_edge_size)

        subset_size = random.randint(min_edge_size, max_edge_size)

        # Randomly sample `subset_size` unique elements from [0 .. max_vertex] and sort them
        subset = sorted(random.sample(range(max_vertex + 1), subset_size))

        # Add the sorted subset to the list
        subsets.append(subset)

    return subsets


random_sampled_subsets(max_vertex=10, num_edges=3,
                       min_edge_size=3, max_edge_size=5)

random_seed = 0
max_vertex = 50
num_edges = 13
num_hypergraphs = 200
min_edge_size = 6
max_edge_size = 16

hypergraphs = [random_sampled_subsets(max_vertex=max_vertex, num_edges=num_edges, min_edge_size=min_edge_size,
                                      max_edge_size=max_edge_size, random_seed=random_seed + p) for p in range(num_hypergraphs)]



# =============================================================================
# WARM UP
# =============================================================================


start_time = time.time()
interval_decomposition = oat.rust.homology_decomposition_for_union_zigzag_of_hypergraphs_z2_coefficients(
    [
        [[0,1],[2,3]],
        [[3,4],[4,5]],        
    ], 
    max_homology_dimension=3,
    return_cycle_representatives=False,
    print_profiling_statistics=False,
)

print(f"Elapsed time: {time.time() - start_time} seconds")    
print(interval_decomposition)


# =============================================================================
# BENCHMARK
# =============================================================================


start_time = time.time()
interval_decomposition = oat.rust.homology_decomposition_for_union_zigzag_of_hypergraphs_z2_coefficients(
    hypergraphs, 
    max_homology_dimension=3,
    return_cycle_representatives=False,
    print_profiling_statistics=False,
)

print(f"Elapsed time: {time.time() - start_time} seconds")    
print(interval_decomposition)