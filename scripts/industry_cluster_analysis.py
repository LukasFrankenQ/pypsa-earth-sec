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
import matplotlib.pyplot as plt

from typing import Iterable
from collections import deque
from itertools import product
from itertools import combinations
from sklearn.cluster import KMeans
from pyproj import Proj, Transformer
from shapely.geometry import MultiPoint
from scipy.spatial import ConvexHull
from scipy.spatial.distance import euclidean



import sys
from pathlib import Path
sys.path.append(str(Path.cwd().parent / 'notebooks'))
from vincenty import V_inv


#####################        TEMPORARY FUNCTIONS        #####################
# def get_cost_evalator_cost(lower, upper, demand):
#     opex_factor = 0.01
#     cost = (lower + (upper - lower)) / demand
#     return cost, cost * opex_factor


'''
def get_plant_capex(total_demand):

    # expects demand value in MWhth for annual conumption
    # return a (ballpark figure) of the capex in $/MWth
    # based on typical values for EGS plants in $/kW
    # assuming that the plant is running at full capacity

    # further assumed economies of scale:
    # at 1 MWth, capex is $7000/kWth
    # at 100 MWth, capex is $2000/kWth
    # Define the known data points

    # returns cost in $/MWth

    def get_cost(cap):
        if cap > 100:
            return 2000
        else:
            # linear interpolation
            return 7000 - 50 * cap

    avg_demand = total_demand / 8760
    

    if isinstance(avg_demand, Iterable):
        return [get_cost(d) * 1000 for d in avg_demand] 
    else:
        return get_cost(avg_demand) * 1000
'''


def round_borehole_capex(cap, one_drill_cap, one_drill_capex):

    fraction = cap / one_drill_cap
    n_wells = round(fraction)

    return one_drill_capex * max(n_wells, fraction) / min(n_wells, fraction)

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


def convex_hull_diameter(points):
    """
    Computes the diameter (maximum distance between any two points) of the convex hull.
    
    Parameters:
    - points: numpy array of points on the convex hull

    Returns:
    - max_distance: maximum distance between any two points on the convex hull
    """
    n = len(points)
    if n == 1:
        return 0.0
    elif n == 2:
        return euclidean(points[0], points[1])
    else:
        # Rotating calipers algorithm
        max_distance = 0.0
        k = 1  # Initialize antipodal point
        for i in range(n):
            j = (i + 1) % n
            while True:
                next_k = (k + 1) % n
                area = triangle_area(points[i], points[j], points[next_k])
                next_area = triangle_area(points[i], points[j], points[k])
                if area > next_area:
                    k = next_k
                else:
                    break
            distance = euclidean(points[i], points[k])
            if distance > max_distance:
                max_distance = distance
        return max_distance


def triangle_area(a, b, c):
    """
    Computes the area of a triangle given its vertices a, b, c.

    Parameters:
    - a, b, c: numpy arrays representing the vertices of the triangle

    Returns:
    - area: the area of the triangle
    """
    return 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - 
                     (c[0] - a[0]) * (b[1] - a[1]))


def cluster_points(points, indices, threshold):
    """
    Recursively clusters points based on maximum pairwise distance.

    Parameters:
    - points: numpy array of shape (n_samples, n_features)
    - indices: list of indices corresponding to the points
    - threshold: maximum allowable distance within a cluster

    Returns:
    - clusters: list of tuples, each containing indices of a cluster
    """
    if len(points) <= 1:
        return [tuple(indices)]

    # Compute the convex hull
    # if len(points) >= 3:
    #     points += 1e-6 * np.random.randn(*points.shape)  # Add noise to avoid colinearity
    #     hull = ConvexHull(points)
    #     hull_points = points[hull.vertices]
    # else:
    #     hull_points = points  # For 2 points, the convex hull is the line segment itself

    # Compute the diameter of the convex hull using rotating calipers algorithm
    # max_distance = convex_hull_diameter(hull_points)

    from scipy.spatial.distance import pdist

    # Compute pairwise distances
    pairwise_distances = pdist(points)
    max_distance = pairwise_distances.max()

    if max_distance <= threshold:
        return [tuple(indices)]
    else:
        kmeans = KMeans(n_clusters=2)
        labels = kmeans.fit_predict(points)

        cluster1_points = points[labels == 0]
        cluster1_indices = [indices[i] for i in range(len(indices)) if labels[i] == 0]

        cluster2_points = points[labels == 1]
        cluster2_indices = [indices[i] for i in range(len(indices)) if labels[i] == 1]

        clusters = []
        clusters.extend(cluster_points(cluster1_points, cluster1_indices, threshold))
        clusters.extend(cluster_points(cluster2_points, cluster2_indices, threshold))

        return clusters


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


# pipe_capex = 500000 # $/MWkm
# pipe_capex = 1000 # $/MWkm
# pipe_capex = 1700 # $/m

from networkx.algorithms.approximation.traveling_salesman import traveling_salesman_problem

def get_simple_costoptimal_network(sites, pipe_capex=1700):

    G = nx.Graph()

    for i, (x, y) in enumerate(sites):
        G.add_node(
            i,
            pos=(x, y),
            )
    hold = G.copy()
    
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):

            length = np.linalg.norm(np.array(sites[i]) - np.array(sites[j]))
            edge_weight = length * pipe_capex * 2 # 2 because there is a hot and cold cycle

            G.add_edge(i, j, weight=edge_weight)

    T = nx.minimum_spanning_tree(G, weight='weight', algorithm='kruskal')
    # T = nx.minimum_spanning_tree(G, weight='weight', algorithm='prim')
    total_cost = 0
    for _, _, data in T.edges(data=True):
        total_cost += data.get('weight', 0)

    # G.remove_edges_from(nx.selfloop_edges(G))
    # node_list = traveling_salesman_problem(G, weight='weight')
    # total_cost = 0

    # for u, v in zip(node_list[:-1], node_list[1:]):
        # hold.add_edge(u, v)
        # total_cost += G.get_edge_data(u, v)['weight']

    return hold, total_cost



def get_costoptimal_network(sites, temps, caps):

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

            dist = V_inv(p1, p2)[0]
            # total_pipe_cost += value / 8760 * np.linalg.norm(np.array(p1) - np.array(p2))
            total_pipe_cost += value / 8760 * dist

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


def get_heat_network(
        xy: np.array,
        temps: list,
        caps: list,
        *args,
        pipe_capex=1700):
    '''
    args are arguments passed to fct round_bolehole_capex
    '''

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

        # return caps[0], get_plant_capex(np.array(caps))[0]
        return caps[0], round_borehole_capex(caps[0], *args)

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

            num_well_test_sqrt = 3

            # xrange = np.linspace(xy[:, 1].min(), xy[:, 1].max(), num_well_test_sqrt)
            # yrange = np.linspace(xy[:, 0].min(), xy[:, 0].max(), num_well_test_sqrt)
            xrange = np.linspace(xy[:, 1].min(), xy[:, 1].max(), num_well_test_sqrt)[[1]]
            yrange = np.linspace(xy[:, 0].min(), xy[:, 0].max(), num_well_test_sqrt)[[1]]

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

            # G, pipe_volume = get_costoptimal_network(layout, temps + [max(temps) + 1], caps + [0])
            G, pipe_volume = get_simple_costoptimal_network(layout, pipe_capex)

            well_xs.append(well_coords[1])
            well_ys.append(well_coords[0])
            pipe_lengths.append(pipe_volume)

        if not pipe_lengths:
            return sum(caps), np.nan

        results = pd.DataFrame({
            'pipe_lengths': pipe_lengths,
            'well_xs': well_xs,
            'well_ys': well_ys
        })

        results = results.sort_values('pipe_lengths').iloc[0]
        pipe_cost = results.loc['pipe_lengths'] * pipe_capex

        # return sum(caps), pipe_cost#  + get_plant_capex(sum(caps))
        return sum(caps), pipe_cost / sum(caps) + round_borehole_capex(sum(caps), *args)


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
    
    # for i, (x, y) in enumerate(sites):
    #     ax.text(x, y, f'{i:.0f}', ha='center', va='center')

    # ax.set_title(f'Total Pipe Cost: {total_pipe_cost:.2f}')
    plt.show()


def coords_to_relative_utm(coords):
    """
    Transforms a list of longitude and latitude coordinates to UTM coordinates relative to the centroid.

    Parameters:
    - coords: list of tuples
        List containing (latitude, longitude) tuples.

    Returns:
    - relative_coords_km: numpy.ndarray
        Array of transformed coordinates in kilometers relative to the centroid.
    """
    coords_array = np.array(coords)
    lons = coords_array[:, 1]
    lats = coords_array[:, 0]

    centroid_lon = np.mean(lons)
    centroid_lat = np.mean(lats)

    utm_zone = int(np.floor((centroid_lon + 180) / 6) % 60) + 1
    hemisphere = 'north' if centroid_lat >= 0 else 'south'

    utm_proj = Proj(proj='utm', zone=utm_zone, hemisphere=hemisphere)

    transformer = Transformer.from_proj(
        proj_from='epsg:4326',  # WGS84 Latitude and Longitude
        proj_to=utm_proj,
        always_xy=True
    )

    eastings, northings = transformer.transform(lons, lats)
    centroid_easting, centroid_northing = transformer.transform(centroid_lon, centroid_lat)

    relative_eastings = eastings - centroid_easting
    relative_northings = northings - centroid_northing

    relative_coords = np.column_stack((relative_eastings, relative_northings))

    return relative_coords
