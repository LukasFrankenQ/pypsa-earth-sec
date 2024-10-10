'''
This algorithm estimates the supply curve to meet industrial heating demands of
varying temperatures and volumes through EGS. It identifies potential clusters of
sites that are eligible for a heat-network (based on the passed `maximum_network_extent`).
These networks can achieve better economies of scale, which can outweigh the additional
cost of piping, which is also considered here.
For each cluster, the algorithm estimates some basic properties of the cost-optimal piping
network that is suitable to meet all heat demands in the network. It only considers
layouts where temperature is highest at the EGS site, and monotonically decreases as heat
is distributed through the network. The location of the EGS plant also minimizes the
required piping volume.
From the size of the demand, each cluster is assigned an investment (and operational) cost
(see `get_egs_plant_cost()`) of a hypothetical EGS plant. Demands and investment/operations
costs are then compiled into an EGS supply curve, which can be inserted into the energy
system model.
For full functionality, the algorithm should be enabled to query a method provided by
Project InnerSpace that returns investment and operational cost of an EGS plant based on
plant coordinates, size and temperature needs.

'''


import numpy as np
import pandas as pd
import networkx as nx
import geopandas as gpd
from collections import deque
from itertools import product
from itertools import combinations
from shapely.geometry import MultiPoint

import matplotlib.pyplot as plt


#####################        TEMPORARY FUNCTIONS        #####################
def get_cost_evalator_cost(lower, upper, demand):
    opex_factor = 0.01
    cost = (lower + (upper - lower)) / demand
    return cost, cost * opex_factor


def get_drilling():

    default_gradient = 250 / 10
    gradient_dev = np.random.normal(1, 0.2)

    drilling_cost_factor = 100

    depth = np.arange(0, 10, 100)
    temp = default_gradient * gradient_dev * depth

    cost = drilling_cost_factor * depth

    return depth, cost, temp


def get_plant_capex(demand):
    return 500 + 1000 / (demand * 2)


##################         CLUSTERING HELPERS          ####################
'''
def is_monotonically_decreasing(G, values):
    if isinstance(values, list):
        values = {i: values[i] for i in range(len(values))}

    start_node = max(values, key=values.get)

    prev_value = float('inf')

    for node in nx.dfs_preorder_nodes(G, start_node):
        if values[node] >= prev_value:
            return False
        prev_value = values[node]

    return True
'''


def get_fc_graph(points, min_widths):
    G = nx.Graph()

    for i, (x, y) in enumerate(points):
        G.add_node(
            i,
            pos=(x, y),
            min_width=min_widths[i],
            )

    for i in range(len(points)):
        for j in range(i + 1, len(points)):

            length = np.linalg.norm(np.array(points[i]) - np.array(points[j]))
            width = max(min_widths[i], min_widths[j])
            cost = length * width

            G.add_edge(i, j, weight=cost)
    
    return G


def measure_graph(G):
    total = 0

    for u, v in G.edges():
        p1 = np.array(G.nodes[u]['pos'])
        p2 = np.array(G.nodes[v]['pos'])

        total += np.linalg.norm(p1 - p2)
    
    return total


'''
def constrained_mst(G, temp):
    """
    Computes the Minimum Spanning Tree (MST) of a graph with an additional constraint.

    Parameters:
    - G: A NetworkX graph (assumed to be fully connected).
    - constraint_func: A function that takes (G, u, v, data) and returns True or False.

    Returns:
    - mst: A NetworkX graph representing the MST.
    """
    # Create a list of all edges with their weights
    edges = [(u, v, data) for u, v, data in G.edges(data=True)]
    # Sort edges by weight
    edges.sort(key=lambda x: x[2].get('weight', 1))

    # Initialize a new graph to store the MST
    mst = nx.Graph()
    mst.add_nodes_from(G.nodes(data=True))

    # Initialize disjoint sets for cycle detection
    parent = {node: node for node in G.nodes()}

    def find(node):
        # Find the root parent of the node
        while parent[node] != node:
            parent[node] = parent[parent[node]]  # Path compression
            node = parent[node]
        return node

    def union(u_root, v_root):
        # Merge two subsets
        parent[u_root] = v_root

    for u, v, data in edges:
        u_root = find(u)
        v_root = find(v)

        # Check if adding this edge creates a cycle
        hold = G.copy()
        hold.add_edge(u, v, **data)

        if u_root != v_root and is_monotonically_decreasing(hold, temp):
                mst.add_edge(u, v, **data)
                union(u_root, v_root)

    return mst
'''
    

def prufer_decode(prufer_sequence):
    m = len(prufer_sequence)
    n = m + 2
    degree = [1] * n

    for node in prufer_sequence:
        degree[node] += 1

    edges = []
    leaves = deque([i for i in range(n) if degree[i] == 1])

    for node in prufer_sequence:
        leaf = leaves.popleft()

        edges.append((leaf, node))

        degree[leaf] -= 1
        degree[node] -= 1

        if degree[node] == 1:
            leaves.append(node)
            leaves = deque(sorted(leaves))

    u = leaves.popleft()
    v = leaves.popleft()
    edges.append((u, v))

    return edges


def generate_all_trees(n):

    prufer_sequences = product(range(n), repeat=n-2)
    all_trees = []

    for prufer_seq in prufer_sequences:
        edges = prufer_decode(list(prufer_seq))
        all_trees.append(edges)

    return all_trees


def can_traverse_monotonically_decreasing(tree, values):
    max_value = max(values)
    start_nodes = [node for node, val in enumerate(values) if val == max_value]

    def dfs(current_node, current_value, visited):
        visited.add(current_node)
        for neighbor in tree.get(current_node, []):
            if neighbor not in visited and values[neighbor] <= current_value:
                dfs(neighbor, values[neighbor], visited)

    for start_node in start_nodes:
        visited = set()
        dfs(start_node, values[start_node], visited)
        if len(visited) == len(tree):
            return True

    return False


def get_costoptimal_network(sites, temps, caps):

    pipe_cost = 0.1 # $/GWkm

    # function to get the capacity needed for each pipe
    def compute_subtree_sums(tree, values, current_node, parent_node, edge_values):

        subtree_sum = values[current_node]
        
        for child in tree.get(current_node, []):
            if child != parent_node:

                child_subtree_sum = compute_subtree_sums(tree, values, child, current_node, edge_values)
                edge_values[(current_node, child)] = child_subtree_sum
                subtree_sum += child_subtree_sum

        return subtree_sum

    if isinstance(temps, np.ndarray):
        temps = list(temps)
    if isinstance(caps, np.ndarray):
        caps = list(caps)
    
    n = len(temps)
    start_node = temps.index(max(temps))

    all_trees = generate_all_trees(n)

    base = nx.Graph()
    for i, (x, y) in enumerate(sites):
        base.add_node(
            i,
            pos=(x, y),
            cap=caps[i],
            )

    keepers = []
    for tree in all_trees:
        hold = base.copy()

        for u, v in tree:
            hold.add_edge(u, v)

        if can_traverse_monotonically_decreasing(nx.to_dict_of_lists(hold), temps):
            keepers.append(tree)

    network_cost = []
    for i, tree in enumerate(keepers):
        # hold = G.copy()
        hold = base.copy()

        for u, v in tree:
            hold.add_edge(u, v, weight=np.linalg.norm(np.array(sites[u]) - np.array(sites[v])))

        edge_values = {}
        compute_subtree_sums(
            nx.to_dict_of_lists(hold),
            caps,
            start_node,
            None,
            edge_values)

        total_pipe_cost = 0

        for edge, value in edge_values.items():
            p1 = hold.nodes[edge[0]]['pos']
            p2 = hold.nodes[edge[1]]['pos']

            total_pipe_cost += value * np.linalg.norm(np.array(p1) - np.array(p2)) * pipe_cost

        network_cost.append(total_pipe_cost)

    best_network = keepers[np.argmin(network_cost)]    

    for (u, v) in best_network:
        base.add_edge(u, v)

    return base, min(network_cost)


def get_partitions(elements):

    def get_naive_partitions(elements):
        if not elements:
            return [[]]

        result = []
        for i in range(1, len(elements) + 1):
            for combo in combinations(elements, i):
                remaining_elements = [e for e in elements if e not in combo]
                for rest in get_naive_partitions(remaining_elements):
                    result.append([tuple(combo)] + rest)
        
        return result


    all_partitions = get_naive_partitions(elements)

    def remove_permutations(lst):
        seen = list()
        cleaned_list = []

        for element in lst:
            sorted_element = sorted(element)

            if sorted_element not in seen:
                cleaned_list.append(element)
                seen.append(sorted_element)

        return cleaned_list

    return remove_permutations(all_partitions)


def get_heat_network(xy: np.array, temps: list, caps: list, pipe_price: float):
    
    assert len(xy) == len(temps), 'Length of xy and temperatures must be the same'
    assert len(xy) == len(caps), 'Length of xy and capacities must be the same'
    if len(xy):
        assert xy.shape[1] == 2
    
    if not isinstance(temps, list):
        temps = temps.tolist()
    if not isinstance(caps, list):
        caps = caps.tolist()
    n = xy.shape[0]

    if n == 1:

        return caps[0], get_plant_capex(np.array(caps))[0]

    else:

        if n == 2:

            # well placement will always just be at larger hotter demand site; redundant approach?
            well_proposals = gpd.GeoSeries(
                gpd.points_from_xy(
                np.linspace(xy[0, 1], xy[1, 1], 2),
                np.linspace(xy[0, 0], xy[1, 0], 2)
                )
            )

        else:

            num_well_test_sqrt = 10

            xrange = np.linspace(xy[:, 1].min(), xy[:, 1].max(), num_well_test_sqrt)
            yrange = np.linspace(xy[:, 0].min(), xy[:, 0].max(), num_well_test_sqrt)

            well_proposals = gpd.GeoSeries(
                gpd.points_from_xy(
                np.stack([xrange for _ in range(num_well_test_sqrt)], axis=1).flatten(),
                np.stack([yrange for _ in range(num_well_test_sqrt)], axis=1).T.flatten(),
                )
            )

            hull = MultiPoint(xy[:,::-1]).convex_hull
            well_proposals = (
                well_proposals
                .loc[well_proposals.within(hull)]
            )

        pipe_lengths, well_xs, well_ys = [], [], []

        for well in well_proposals:

            well_coords = np.array(well.coords)[0][::-1]

            layout = (
                np.vstack([xy, well_coords])
            )

            G, pipe_volume = get_costoptimal_network(layout, temps + [max(temps) + 1], caps + [0])

            '''
            _, ax = plt.subplots(figsize=(3, 3))
            for u, v in G.edges():
                p1 = G.nodes[u]['pos']
                p2 = G.nodes[v]['pos']
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k-', zorder=1)
            ax.scatter(layout[:, 0], layout[:, 1], c=temps + [max(temps) + 1], cmap='magma', edgecolor='black', linewidth=1)
            plt.show()
            '''

            well_xs.append(well_coords[1])
            well_ys.append(well_coords[0])
            pipe_lengths.append(pipe_volume)

        if not pipe_lengths:
            return np.nan

        results = pd.DataFrame({
            'pipe_lengths': pipe_lengths,
            'well_xs': well_xs,
            'well_ys': well_ys
        })

        results = results.sort_values('pipe_lengths').iloc[0]
        pipe_cost = results.loc['pipe_lengths'] / 360 * 6371 * pipe_price

        return sum(caps), pipe_cost + get_plant_capex(sum(caps))


def plot_network(G, sites, caps, temps):
    """
    Plots a network graph with nodes and edges, and overlays site data with temperature-based coloring.
    Parameters:
    G (networkx.Graph): The network graph containing nodes and edges.
    sites (numpy.ndarray): A 2D array of site coordinates.
    caps (list or numpy.ndarray): A list or array of capacities for each site, used to scale the marker sizes.
    temps (list or numpy.ndarray): A list or array of temperatures for each site, used to color the markers.
    Returns:
    None
    """

    _, ax = plt.subplots(figsize=(3, 3))

    for u, v in G.edges():
        p1 = G.nodes[u]['pos']
        p2 = G.nodes[v]['pos']

        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k-', zorder=1)

    ax.scatter(
        sites[:, 0],
        sites[:, 1],
        c=temps,
        cmap='coolwarm',
        edgecolor='black',
        linewidth=1,
        s=np.array(caps) * 30,
        )

    # ax.set_title(f'Total Pipe Cost: {total_pipe_cost:.2f}')
    plt.show()