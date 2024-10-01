
#### Building an EGS Supply Curve for Industrial Heating Demands

- __define__ ``piping_cost``        [ __$/MWth/km__ ] <span style="color: red;">InnerSpace</span>

- __define__ ``maximum_network_extent``        [ __km__ ] <span style="color: red;">InnerSpace</span>
- __define__ ``get_egs_plant_cost()`` [ __$/kWth__ ] __#__ evaluates cost of geothermal site based on ``x``, ``y``, ``demand``, ``temperature`` <span style="color: red;">Data and Method InnerSpace, Implementation Lukas</span> <br><br>

#### Algorithm Summary

This algorithm estimates the supply-curve to meet industrial heating demands of varying temperatures and volumes through EGS. It identifies potential clusters of sites that are eligible for a heat-network (based on the passed `maximum_network_extent`). These networks can achieve better economies of scale, which can outwigh the additional cost of piping which is also considered here.
For each cluster, the algorithm estimates some basic properties of the cost-optimal piping network that is suitable to meet all heat demands in the network. It only considers layouts where temperature is highest at the EGS site, and monotonically decreases as heat is distributed through the network. The location of the EGS plant also minimizes the required piping volume.
From the size of the demand, each cluster is assigned an investment (and operational) cost (see `get_egs_plant_cost()`) of a hypothetical EGS plant. Demands and investment/operations costs are then compiled into an EGS supply curve, that can be inserted into the energy system model.
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
