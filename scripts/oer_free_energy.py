"""
oer_free_energy.py

Analyze a conventional four-step oxygen evolution reaction (OER)
using the computational hydrogen electrode (CHE) framework and
generate a publication-style free-energy diagram.

The script reports:
    - Free energy of each elementary OER step at U = 0 V
    - Potential-determining step (PDS)
    - Limiting potential
    - OER overpotential

The script generates:
    - oer_free_energy.png
    - oer_summary.csv

Usage:
    python3 oer_free_energy.py DG1 DG2 DG3 DG4

Example:
    python3 oer_free_energy.py 1.10 1.45 1.80 0.57

Author: Radha Somaiya
"""

import argparse
import csv

import matplotlib

# Non-interactive backend for HPC/remote systems
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ============================================================
# OER THERMODYNAMIC ANALYSIS
# ============================================================

def calculate_oer_metrics(delta_g):
    """
    Calculate basic OER thermodynamic quantities.

    Parameters
    ----------
    delta_g : list of float
        Free energies of the four elementary OER steps
        at U = 0 V, in eV.

    Returns
    -------
    dict
        Potential-determining step, limiting potential,
        and OER overpotential.
    """

    if len(delta_g) != 4:
        raise ValueError(
            "Exactly four OER step free energies are required."
        )

    maximum_delta_g = max(delta_g)

    pds_index = delta_g.index(
        maximum_delta_g
    )

    # Each conventional OER step transfers one electron.
    # When Delta G is expressed in eV, the numerical value
    # of the limiting potential is obtained in volts.
    limiting_potential = maximum_delta_g

    equilibrium_potential = 1.23

    overpotential = (
        limiting_potential
        - equilibrium_potential
    )

    return {
        "pds_step": pds_index + 1,
        "maximum_delta_g": maximum_delta_g,
        "limiting_potential": limiting_potential,
        "overpotential": overpotential,
    }


# ============================================================
# CUMULATIVE FREE ENERGIES
# ============================================================

def cumulative_free_energies(delta_g):
    """
    Convert individual OER step free energies into
    cumulative free energies.
    """

    energies = [0.0]

    total = 0.0

    for value in delta_g:

        total += value
        energies.append(total)

    return energies


# ============================================================
# POTENTIAL DEPENDENCE
# ============================================================

def apply_potential(energies, potential):
    """
    Apply CHE potential dependence to cumulative OER energies.

    For the conventional four-step OER:

        G_n(U) = G_n(0) - nU

    where n is the number of transferred proton-electron pairs.
    """

    shifted_energies = []

    for n, energy in enumerate(energies):

        shifted_energy = (
            energy
            - n * potential
        )

        shifted_energies.append(
            shifted_energy
        )

    return shifted_energies


# ============================================================
# DRAW FREE-ENERGY PROFILE
# ============================================================

def draw_profile(
    ax,
    energies,
    color="black",
    connector_style="--",
    line_width=2.5,
    connector_width=1.5,
):
    """
    Draw a publication-style free-energy profile.
    """

    number_of_states = len(energies)

    level_half_width = 0.28

    for i, energy in enumerate(energies):

        # Solid horizontal intermediate level
        ax.plot(
            [
                i - level_half_width,
                i + level_half_width,
            ],
            [
                energy,
                energy,
            ],
            color=color,
            linewidth=line_width,
            linestyle="-",
            solid_capstyle="butt",
        )

        # Connector to next intermediate
        if i < number_of_states - 1:

            next_energy = energies[i + 1]

            ax.plot(
                [
                    i + level_half_width,
                    i + 1 - level_half_width,
                ],
                [
                    energy,
                    next_energy,
                ],
                color=color,
                linewidth=connector_width,
                linestyle=connector_style,
            )


# ============================================================
# FREE-ENERGY DIAGRAM
# ============================================================

def plot_oer_free_energy(
    delta_g,
    filename="oer_free_energy.png",
):
    """
    Generate publication-style OER free-energy diagrams.

    Profiles:
        U = 0 V
            black solid intermediate levels
            black dashed connectors

        U = 1.23 V
            red solid intermediate levels
            red dotted connectors
    """

    energies_0 = cumulative_free_energies(
        delta_g
    )

    energies_123 = apply_potential(
        energies_0,
        1.23,
    )

    states = [
        r"$*$",
        r"$*\mathrm{OH}$",
        r"$*\mathrm{O}$",
        r"$*\mathrm{OOH}$",
        r"$*+\mathrm{O_2}$",
    ]

    fig, ax = plt.subplots(
        figsize=(7.5, 5.5)
    )

    # U = 0 V profile
    draw_profile(
        ax,
        energies_0,
        color="black",
        connector_style="--",
        line_width=2.5,
        connector_width=1.5,
    )

    # U = 1.23 V profile
    draw_profile(
        ax,
        energies_123,
        color="red",
        connector_style=":",
        line_width=2.5,
        connector_width=1.8,
    )

    # Energy labels for U = 0 V
    for i, energy in enumerate(
        energies_0
    ):

        ax.text(
            i,
            energy + 0.12,
            f"{energy:.2f}",
            ha="center",
            va="bottom",
            fontsize=9,
            color="black",
        )

    # Legend
    legend_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linewidth=2.5,
            linestyle="-",
            label=r"$U = 0$ V",
        ),
        Line2D(
            [0],
            [0],
            color="red",
            linewidth=2.5,
            linestyle="-",
            label=r"$U = 1.23$ V",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        loc="upper left",
        frameon=False,
        fontsize=10,
        handlelength=3.0,
    )

    # Axis limits
    ax.set_xlim(
        -0.55,
        4.55,
    )

    all_energies = (
        energies_0
        + energies_123
    )

    minimum_energy = min(
        all_energies
    )

    maximum_energy = max(
        all_energies
    )

    ax.set_ylim(
        minimum_energy - 0.45,
        maximum_energy + 0.55,
    )

    # X axis
    ax.set_xticks(
        range(len(states))
    )

    ax.set_xticklabels(
        states,
        fontsize=11,
    )

    ax.set_xlabel(
        "Reaction pathway",
        fontsize=12,
    )

    # Y axis
    ax.set_ylabel(
        r"$\Delta G$ (eV)",
        fontsize=12,
    )

    # Tick formatting
    ax.tick_params(
        axis="x",
        which="both",
        bottom=True,
        top=False,
        labelbottom=True,
        labeltop=False,
        direction="out",
        width=1.2,
        length=4,
        labelsize=10,
    )

    ax.tick_params(
        axis="y",
        which="both",
        left=True,
        right=False,
        labelleft=True,
        labelright=False,
        direction="out",
        width=1.2,
        length=4,
        labelsize=10,
    )

    # Border thickness
    for spine in ax.spines.values():

        spine.set_linewidth(
            1.3
        )

    # Save figure
    fig.tight_layout()

    fig.savefig(
        filename,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)


# ============================================================
# CSV EXPORT
# ============================================================

def save_csv(
    delta_g,
    results,
    filename="oer_summary.csv",
):
    """
    Save OER step free energies and thermodynamic metrics
    to a CSV file.
    """

    step_names = [
        "* -> *OH",
        "*OH -> *O",
        "*O -> *OOH",
        "*OOH -> * + O2",
    ]

    with open(
        filename,
        "w",
        newline="",
    ) as csvfile:

        writer = csv.writer(
            csvfile
        )

        writer.writerow(
            [
                "Step",
                "Reaction",
                "DeltaG_eV",
            ]
        )

        for i, (
            step_name,
            value,
        ) in enumerate(
            zip(
                step_names,
                delta_g,
            ),
            start=1,
        ):

            writer.writerow(
                [
                    i,
                    step_name,
                    f"{value:.6f}",
                ]
            )

        writer.writerow([])

        writer.writerow(
            [
                "Metric",
                "Value",
                "Unit",
            ]
        )

        writer.writerow(
            [
                "Potential-determining step",
                results["pds_step"],
                "",
            ]
        )

        writer.writerow(
            [
                "Maximum Delta G",
                f"{results['maximum_delta_g']:.6f}",
                "eV",
            ]
        )

        writer.writerow(
            [
                "Limiting potential",
                f"{results['limiting_potential']:.6f}",
                "V",
            ]
        )

        writer.writerow(
            [
                "OER overpotential",
                f"{results['overpotential']:.6f}",
                "V",
            ]
        )


# ============================================================
# PRINT RESULTS
# ============================================================

def print_results(
    delta_g,
    results,
):
    """
    Print OER thermodynamic results.
    """

    step_names = [
        "* -> *OH",
        "*OH -> *O",
        "*O -> *OOH",
        "*OOH -> * + O2",
    ]

    print()

    print(
        "OER Free-Energy Analysis"
    )

    print(
        "------------------------"
    )

    for i, (
        step_name,
        value,
    ) in enumerate(
        zip(
            step_names,
            delta_g,
        ),
        start=1,
    ):

        print(
            f"Step {i} "
            f"({step_name}): "
            f"{value:.3f} eV"
        )

    print()

    print(
        "Potential-determining step: "
        f"Step {results['pds_step']}"
    )

    print(
        "Maximum Delta G:             "
        f"{results['maximum_delta_g']:.3f} eV"
    )

    print(
        "Limiting potential:          "
        f"{results['limiting_potential']:.3f} V"
    )

    print(
        "OER overpotential:           "
        f"{results['overpotential']:.3f} V"
    )

    print()

    print(
        "Free-energy diagram saved to: "
        "oer_free_energy.png"
    )

    print(
        "Results saved to: "
        "oer_summary.csv"
    )


# ============================================================
# COMMAND-LINE INPUT
# ============================================================

def parse_arguments():
    """
    Read the four OER step free energies from the command line.

    Example:

        python3 oer_free_energy.py 1.10 1.45 1.80 0.57
    """

    parser = argparse.ArgumentParser(
        description=(
            "Analyze four OER free-energy steps and generate "
            "a publication-style free-energy diagram."
        )
    )

    parser.add_argument(
        "delta_g1",
        type=float,
        help="Delta G for * -> *OH in eV",
    )

    parser.add_argument(
        "delta_g2",
        type=float,
        help="Delta G for *OH -> *O in eV",
    )

    parser.add_argument(
        "delta_g3",
        type=float,
        help="Delta G for *O -> *OOH in eV",
    )

    parser.add_argument(
        "delta_g4",
        type=float,
        help="Delta G for *OOH -> * + O2 in eV",
    )

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():
    """
    Main program.
    """

    # Read command-line input
    args = parse_arguments()

    delta_g = [
        args.delta_g1,
        args.delta_g2,
        args.delta_g3,
        args.delta_g4,
    ]

    # Thermodynamic analysis
    results = calculate_oer_metrics(
        delta_g
    )

    # Generate free-energy figure
    plot_oer_free_energy(
        delta_g,
        "oer_free_energy.png",
    )

    # Save CSV results
    save_csv(
        delta_g,
        results,
        "oer_summary.csv",
    )

    # Print terminal summary
    print_results(
        delta_g,
        results,
    )


if __name__ == "__main__":
    main()
