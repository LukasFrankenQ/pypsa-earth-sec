
#### Building an EGS Supply Curve for Industrial Heating Demands

- __define__ ``piping_cost``        [ __$/MWth/km__ ] <span style="color: red;">InnerSpace</span>

- __define__ ``maximum_network_extent``        [ __km__ ] <span style="color: red;">InnerSpace</span>
- __define__ ``get_egs_plant_cost()`` [ __$/kWth__ ] __#__ evaluates cost of geothermal site based on ``x``, ``y``, ``demand``, ``temperature`` <span style="color: red;">Data and Method InnerSpace, Implementation Lukas</span> <br><br>

#### Algorithm Summary

This algorithm estimates the supply curve to meet industrial heating demands of varying temperatures and volumes through EGS. It identifies potential clusters of sites that are eligible for a heat-network (based on the passed `maximum_network_extent`). These networks can achieve better economies of scale, which can outweigh the additional cost of piping, which is also considered here.
For each cluster, the algorithm estimates some basic properties of the cost-optimal piping network that is suitable to meet all heat demands in the network. It only considers layouts where temperature is highest at the EGS site, and monotonically decreases as heat is distributed through the network. The location of the EGS plant also minimizes the required piping volume.
From the size of the demand, each cluster is assigned an investment (and operational) cost (see `get_egs_plant_cost()`) of a hypothetical EGS plant. Demands and investment/operations costs are then compiled into an EGS supply curve, which can be inserted into the energy system model.
For full functionality, the algorithm should be enabled to query a method provided by Project InnerSpace that returns investment and operational cost of an EGS plant based on plant coordinates, size and temperature needs.


#### Pseudocode


- __define__ ``region_supply_curves`` ← empty dict <br><br>
- __for__ ``region`` in network_regions:
    - ``demand_sites`` ← _load_region_sites_(``region``)

    - __define__ ``supply_potential`` ← [ ]
    - __define__ ``supply_capex`` ← [ ]
    - `clusters` ← _DBScan_(``demand_sites``, ``max_network_extent``)

    - __for__ ``cluster`` in ``clusters`` where size of cluster == 1:
        
        - __define__ ``plant_capex`` ← ``get_egs_plant_cost(cluster)``
        - ``supply_potential``._append_(``cluster_peak_demand``)         __# skipped decision__
        - ``supply_capex``._append_(``plant_capex``)


    - __for__ ``cluster`` in ``clusters`` where size of cluster > 1:

        - __define__ ``partitions`` ← _get_cluster_partitions_(``cluster``)
        - __define__ ``partition_qualities`` ← [ ]  

        - __for__ ``partition`` in ``partitions``:

            - __define__ ``network_costs`` ← [ ]
            - __define__ ``network_capacities``  ← [ ]

            - __for__ ``heat_network`` in ``partition``:

                - __if__ number of sites in heat_network == 1:
                    - ``network_capacities``._append_(``heat_network_peak_demand``)         __# skipped decision__
                    - ``network_costs``._append_(``plant_capex``)

                - __if__ number of sites in heat_network > 1:
                    - // greedy approach: tests all potential site location within convex hull spanned by sites in network
                    - ``heat_network`` ← _get_cost_optimal_heating_network_(``site_locations``, ``site_demands``, ``piping_cost``, ``temperatures``)  __# returns lowest cost piping configuration while ensuring monotonically decreasing temperature__

                    - ``network_capacities``._append_(``heat_network_peak_demand``) 
                    - ``network_costs``._append_(_get_plant_cost_(``heat_network``) + ``network_piping_costs``)
                - ``partition_qualities``._append_(_sum_(``network_costs``))           

        - __define__ ``partition`` ← _get_best_partition_(``partition_qualities``)
        - _update_supply_(``supply_potential``, ``supply_capex``, ``data_of_best_partition``)

    ``region_supply_curves``[``region``] ← _build_supply_curve_(``supply_potential``, ``supply_capex``)


## New Version

#### Building an EGS Supply Curve for Industrial Heating Demands

- __define__ ``piping_cost``        [ __$/MWth/km__ ] <span style="color: red;">InnerSpace</span>

- __define__ ``maximum_network_extent``        [ __km__ ] <span style="color: red;">InnerSpace</span>
- __define__ ``get_egs_plant_cost()`` [ __$/kWth__ ] __#__ evaluates cost of geothermal site based on ``x``, ``y``, ``demand``, ``temperature`` <span style="color: red;">Data and Method InnerSpace, Implementation Lukas</span> <br><br>

#### Algorithm Summary

This algorithm estimates the supply curve to meet industrial heating demands of varying temperatures and volumes through EGS. It identifies potential clusters of sites that are close enough to each other to be eligible for a heat-network. Supplying demand sites jointly through a heat-network can make EGS, which is only economic from a demand of around 10MWth, viable.
The algorithm first scans each site for potential clusters based on spatial proximity using a _Density-Based Scan_. Each network is split into multiple smaller networks, until the maximum feasible total network extent is undercut.
For each resulting network, the algorithm then computes the most efficient piping layout to distribute heat from an EGS site to the demands. The cost of the resulting piping is added to the cost factors related to the plant to make up the overall investment and operation cost of the heat network.


#### Pseudocode


Algorithm EstimateSupplyCurveEGS

Input:
- DemandSites: List of industrial heating demand sites with their locations, temperature requirements, and volume demands
- MaxNetworkExtent: Maximum feasible total network extent
- MinEconomicDemand: Minimum demand (e.g., 10 MWth) for EGS viability

Output:
- SupplyCurve: Estimated cost to meet industrial heating demands through EGS

Steps:

1. Initialize PotentialClusters as an empty list
2. Filter out demands that requires temperatures above 250C 

3. // Identify potential clusters based on spatial proximity
4. For each site in DemandSites:
    a. Perform Density-Based Scan to find neighboring sites within a proximity threshold
    b. Form clusters of sites that are close enough for a heat network
    c. Add these clusters to PotentialClusters

5. // Split clusters exceeding maximum network extent
6. For each cluster in PotentialClusters:
    a. While TotalExtent(cluster) > MaxNetworkExtent:
        i. Split cluster into smaller networks
       ii. Recalculate TotalExtent for each new network

7. Initialize NetworkCosts as an empty list

8. // Compute efficient piping layout and costs for each network
9. For each network in PotentialClusters:
    a. ComputeMostEfficientPipingLayout(network)
    b. Calculate PipingCost for the network
    c. If TotalDemand(network) ≥ MinEconomicDemand:
        i. Calculate EGSPlantCost based on TotalDemand
       ii. TotalCost = PipingCost + EGSPlantCost
      iii. Add (network, TotalCost) to NetworkCosts

10. // Aggregate costs to form the supply curve
11. SupplyCurve = AggregateCosts(NetworkCosts)

12. Return SupplyCurve