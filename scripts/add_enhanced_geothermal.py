import pypsa
import pandas as pd

from helpers import override_component_attrs


def add_enhanched_geothermal(n):

    egs_potential = pd.read_csv(snakemake.input["egs_potential"], index_col=[0,1])

    idx = pd.IndexSlice

    n.add(
        "Bus",
        "EGS",
        carrier="geothermal heat",
        unit="MWh_th",
    )

    n.add(
        "Generator",
        "EGS",
        bus="EGS",
        carrier="geothermal heat",
        p_nom_extendable=True,
    )

    for bus in egs_potential.index.get_level_values(0).unique():

        ss = egs_potential.loc[idx[bus,:]]
        nodes = f"{bus} " + pd.Index(range(len(ss)), dtype=str)

        p_nom_max = ss["potential"].values
        capex = ss.index.values
        opex = ss["opex"].values
        # eta = ss["efficiency"].values
        eta = 0.12

        n.madd(
            "Bus",
            nodes,
            suffix=" EGS surface",
            carrier="geothermal heat",
        )

        n.madd(
            "Link",
            nodes,
            suffix=" EGS well",
            bus0="EGS",
            bus1=nodes + " EGS surface",
            p_nom_max=p_nom_max / eta,
            capital_cost=capex * eta,
            p_nom_extendable=True,
        )

        n.madd(
            "StorageUnit",
            nodes,
            suffix=" EGS reservoir",
            bus=nodes + " EGS surface",
            max_hours=100, # should be agreed on, constraint to be implemented
        )

        n.madd(
            "Link",
            nodes,
            suffix=" EGS surface",
            bus0=nodes + " EGS surface",
            bus1=bus,
            carrier="orc",
            efficiency=eta,
            marginal_cost=opex,
            p_nom_extendable=True,
        )


if __name__ == "__main__":

    # Load input network
    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    add_enhanched_geothermal(n)

    n.export_to_netcdf(snakemake.output["network"])