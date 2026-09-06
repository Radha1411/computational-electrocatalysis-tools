"""
oer_limiting_potential.py

Analyze a conventional four-step oxygen evolution reaction (OER)
using the computational hydrogen electrode (CHE) framework and
compare the free-energy landscape at:

    - U = 0 V
    - U = 1.23 V
    - U = U_L

where U_L is the limiting potential.

The script reports:
    - Free energy of each elementary OER step at U = 0 V
    - Potential-determining step (PDS)
    - Limiting potential
    - OER overpotential

The script generates:
    - oer_limiting_potential.png

Usage:
    python3 oer_limiting_potential.py DG1 DG2 DG3 DG4

Example:
    python3 oer_limiting_potential.py 1.10 1.45 1.80 0.57

Author: Radha Somaiya
"""

import argparse

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

    # For the conventional OER mechanism,
    # each step transfers one electron.
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

    Example:

        Delta G = [1.10, 1.45, 1.80, 0.57]

    gives:

        [0.00, 1.10, 2.55, 4.35, 4.92]
    """

    energies = [0.0]

    total = 0.0

    for value in delta_g:

        total += value

        energies.append(
            total
        )

    return energies


# ============================================================
# POTENTIAL DEPENDENCE
# ============================================================

def apply_potential(
    energies,
    potential,
):
    """
    Apply CHE potential dependence to cumulative OER energies.

    For the conventional four-step OER:

        G_n(U) = G_n(0) - nU

    where n is the number of transferred proton-electron pairs.

    Parameters
    ----------
    energies : list of float
        Cumulative free energies at U = 0 V.

    potential : float
        Applied potential in V vs RHE.

    Returns
    -------
    list of float
        Potential-dependent cumulative free energies.
    """

    shifted_energies = []

    for n, energy in enumerate(
        energies
    ):

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
    color,
    connector_style,
    line_width=2.5,
    connector_width=1.5,
):
    """
    Draw a publication-style free-energy profile.
    """

    number_of_states = len(
        energies
    )

    level_half_width = 0.28

    for i, energy in enumerate(
        energies
    ):

        # Horizontal intermediate level
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

        # Connector to the next intermediate
        if i < number_of_states - 1:

            next_energy = (
                energies[i + 1]
            )

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

def plot_limiting_potential_diagram(
    delta_g,
    limiting_potential,
    filename="oer_limiting_potential.png",
):
    """
    Plot OER free-energy profiles at:

        U = 0 V
        U = 1.23 V
        U = U_L
    """

    # --------------------------------------------------------
    # U = 0 V
    # --------------------------------------------------------

    energies_0 = (
        cumulative_free_energies(
            delta_g
        )
    )

    # --------------------------------------------------------
    # U = 1.23 V
    # --------------------------------------------------------

    energies_123 = (
        apply_potential(
            energies_0,
            1.23,
        )
    )

    # --------------------------------------------------------
    # U = limiting potential
    # --------------------------------------------------------

    energies_ul = (
        apply_potential(
            energies_0,
            limiting_potential,
        )
    )

    # --------------------------------------------------------
    # Reaction states
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # U = 0 V profile
    # --------------------------------------------------------

    draw_profile(
        ax,
        energies_0,
        color="black",
        connector_style="--",
        line_width=2.5,
        connector_width=1.5,
    )

    # --------------------------------------------------------
    # U = 1.23 V profile
    # --------------------------------------------------------

    draw_profile(
        ax,
        energies_123,
        color="red",
        connector_style=":",
        line_width=2.5,
        connector_width=1.8,
    )

    # --------------------------------------------------------
    # U = U_L profile
    # --------------------------------------------------------

    draw_profile(
        ax,
        energies_ul,
        color="blue",
        connector_style="-.",
        line_width=2.5,
        connector_width=1.5,
    )

    # --------------------------------------------------------
    # Energy labels for U = 0 V
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

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
        Line2D(
            [0],
            [0],
            color="blue",
            linewidth=2.5,
            linestyle="-",
            label=(
                rf"$U = U_L = "
                rf"{limiting_potential:.2f}$ V"
            ),
        ),
    ]

    ax.legend(
        handles=legend_handles,
        loc="upper left",
        frameon=False,
        fontsize=10,
        handlelength=3.0,
    )

    # --------------------------------------------------------
    # Axis limits
    # --------------------------------------------------------

    ax.set_xlim(
        -0.55,
        4.55,
    )

    all_energies = (
        energies_0
        + energies_123
        + energies_ul
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

    # --------------------------------------------------------
    # X axis
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Y axis
    # --------------------------------------------------------

    ax.set_ylabel(
        r"$\Delta G$ (eV)",
        fontsize=12,
    )

    # --------------------------------------------------------
    # Tick formatting
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Border thickness
    # --------------------------------------------------------

    for spine in ax.spines.values():

        spine.set_linewidth(
            1.3
        )

    # --------------------------------------------------------
    # Save figure
    # --------------------------------------------------------

    fig.tight_layout()

    fig.savefig(
        filename,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)


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
        "OER Limiting-Potential Analysis"
    )

    print(
        "-------------------------------"
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
        "At the limiting potential, "
        "the PDS becomes thermoneutral."
    )

    print()

    print(
        "Free-energy diagram saved to: "
        "oer_limiting_potential.png"
    )


# ============================================================
# COMMAND-LINE INPUT
# ============================================================

def parse_arguments():
    """
    Read four OER step free energies from the command line.

    Example:

        python3 oer_limiting_potential.py \
        1.10 1.45 1.80 0.57
    """

    parser = argparse.ArgumentParser(
        description=(
            "Analyze OER energetics at U = 0 V, "
            "U = 1.23 V, and the limiting potential."
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
    results = (
        calculate_oer_metrics(
            delta_g
        )
    )

    # Generate limiting-potential diagram
    plot_limiting_potential_diagram(
        delta_g,
        results[
            "limiting_potential"
        ],
        "oer_limiting_potential.png",
    )

    # Print terminal results
    print_results(
        delta_g,
        results,
    )


if __name__ == "__main__":
    main()
