"""
oer_dft_thermodynamics.py

Calculate conventional four-step OER thermodynamics directly
from DFT total energies of the clean surface and OER
intermediates.

Required DFT energies:
    E(*)
    E(*OH)
    E(*O)
    E(*OOH)
    E(H2)
    E(H2O)

Thermodynamic corrections for H2, H2O, OH*, O*, and OOH*
are based on the supplementary information of:

Nature Communications (2020) 11, 2522
DOI: 10.1038/s41467-020-16237-1

Free energy is evaluated as:

    G = E_DFT + ZPE + dH - TS + E_solvation

The script calculates:
    - Corrected free energies
    - Four elementary OER free-energy steps
    - Potential-determining step (PDS)
    - Limiting potential
    - OER overpotential
    - Free-energy diagram at U = 0 and 1.23 V
    - CSV summary

Usage:
    python3 oer_dft_thermodynamics.py \
        E_clean E_OH E_O E_OOH E_H2 E_H2O

Example using synthetic energies:
    python3 oer_dft_thermodynamics.py \
        -100.000 -110.000 -105.500 -114.000 \
        -7.000 -14.000

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
# CONSTANTS
# ============================================================

OER_TOTAL_FREE_ENERGY = 4.92
OER_EQUILIBRIUM_POTENTIAL = 1.23


# ============================================================
# THERMODYNAMIC CORRECTIONS
# ============================================================
#
# Nature Communications (2020) 11, 2522
#
# Gcorr = ZPE + dH - TS + E_solvation
#
# Values below follow the rounded values reported in the
# supplementary information.
# ============================================================


# ------------------------------------------------------------
# H2
# ------------------------------------------------------------

ZPE_H2 = 0.27
DH_H2 = 0.09
TS_H2 = 0.41
ESOLV_H2 = 0.00

GCORR_H2 = (
    ZPE_H2
    + DH_H2
    - TS_H2
    + ESOLV_H2
)


# ------------------------------------------------------------
# H2O
# ------------------------------------------------------------

ZPE_H2O = 0.56
DH_H2O = 0.10
TS_H2O = 0.68
ESOLV_H2O = 0.00

GCORR_H2O = (
    ZPE_H2O
    + DH_H2O
    - TS_H2O
    + ESOLV_H2O
)


# ------------------------------------------------------------
# OH*
# ------------------------------------------------------------

ZPE_OH = 0.39
DH_OH = 0.03
TS_OH = 0.03
ESOLV_OH = -0.35

GCORR_OH = (
    ZPE_OH
    + DH_OH
    - TS_OH
    + ESOLV_OH
)


# ------------------------------------------------------------
# O*
# ------------------------------------------------------------

ZPE_O = 0.07
DH_O = 0.03
TS_O = 0.05
ESOLV_O = 0.00

GCORR_O = (
    ZPE_O
    + DH_O
    - TS_O
    + ESOLV_O
)


# ------------------------------------------------------------
# OOH*
# ------------------------------------------------------------

ZPE_OOH = 0.47
DH_OOH = 0.05
TS_OOH = 0.08
ESOLV_OOH = -0.40

GCORR_OOH = (
    ZPE_OOH
    + DH_OOH
    - TS_OOH
    + ESOLV_OOH
)


# ------------------------------------------------------------
# Clean surface
# ------------------------------------------------------------

GCORR_CLEAN = 0.00


# ============================================================
# CORRECT DFT ENERGIES
# ============================================================

def calculate_corrected_energies(
    e_clean,
    e_oh,
    e_o,
    e_ooh,
    e_h2,
    e_h2o,
):
    """
    Apply thermodynamic corrections to DFT total energies.
    """

    g_clean = (
        e_clean
        + GCORR_CLEAN
    )

    g_oh = (
        e_oh
        + GCORR_OH
    )

    g_o = (
        e_o
        + GCORR_O
    )

    g_ooh = (
        e_ooh
        + GCORR_OOH
    )

    g_h2 = (
        e_h2
        + GCORR_H2
    )

    g_h2o = (
        e_h2o
        + GCORR_H2O
    )

    return {
        "G_clean": g_clean,
        "G_OH": g_oh,
        "G_O": g_o,
        "G_OOH": g_ooh,
        "G_H2": g_h2,
        "G_H2O": g_h2o,
    }


# ============================================================
# OER STEP FREE ENERGIES
# ============================================================

def calculate_oer_steps(corrected):
    """
    Calculate conventional four-step OER free energies
    at U = 0 V using the computational hydrogen electrode.

    Step 1:
        * + H2O -> *OH + H+ + e-

    Step 2:
        *OH -> *O + H+ + e-

    Step 3:
        *O + H2O -> *OOH + H+ + e-

    Step 4:
        *OOH -> * + O2 + H+ + e-

    CHE relation:

        G(H+ + e-) = 1/2 G(H2)

    The fourth step is evaluated from overall OER closure:

        DeltaG1 + DeltaG2 + DeltaG3 + DeltaG4 = 4.92 eV

    This avoids requiring a calculated O2 total energy.
    """

    g_clean = corrected["G_clean"]
    g_oh = corrected["G_OH"]
    g_o = corrected["G_O"]
    g_ooh = corrected["G_OOH"]
    g_h2 = corrected["G_H2"]
    g_h2o = corrected["G_H2O"]

    # --------------------------------------------------------
    # Step 1
    # * + H2O -> *OH + H+ + e-
    # --------------------------------------------------------

    delta_g1 = (
        g_oh
        + 0.5 * g_h2
        - g_clean
        - g_h2o
    )

    # --------------------------------------------------------
    # Step 2
    # *OH -> *O + H+ + e-
    # --------------------------------------------------------

    delta_g2 = (
        g_o
        + 0.5 * g_h2
        - g_oh
    )

    # --------------------------------------------------------
    # Step 3
    # *O + H2O -> *OOH + H+ + e-
    # --------------------------------------------------------

    delta_g3 = (
        g_ooh
        + 0.5 * g_h2
        - g_o
        - g_h2o
    )

    # --------------------------------------------------------
    # Step 4
    #
    # Obtained from the total OER free energy.
    # --------------------------------------------------------

    delta_g4 = (
        OER_TOTAL_FREE_ENERGY
        - delta_g1
        - delta_g2
        - delta_g3
    )

    return [
        delta_g1,
        delta_g2,
        delta_g3,
        delta_g4,
    ]


# ============================================================
# OER METRICS
# ============================================================

def calculate_oer_metrics(delta_g):
    """
    Determine PDS, limiting potential, and overpotential.
    """

    maximum_delta_g = max(
        delta_g
    )

    pds_index = delta_g.index(
        maximum_delta_g
    )

    # One electron is transferred in each conventional OER step.
    limiting_potential = (
        maximum_delta_g
    )

    overpotential = (
        limiting_potential
        - OER_EQUILIBRIUM_POTENTIAL
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
    Convert elementary OER step energies into cumulative
    free-energy levels.
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
    Apply CHE potential dependence.

        G_n(U) = G_n(0) - nU

    Each conventional OER step transfers one electron.
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
    Draw horizontal intermediate levels and connecting lines.
    """

    level_half_width = 0.28

    for i, energy in enumerate(
        energies
    ):

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

        if i < len(energies) - 1:

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

def plot_oer_diagram(
    delta_g,
    filename="oer_dft_thermodynamics.png",
):
    """
    Generate OER free-energy profiles at:

        U = 0 V
        U = 1.23 V
    """

    energies_0 = (
        cumulative_free_energies(
            delta_g
        )
    )

    energies_123 = (
        apply_potential(
            energies_0,
            OER_EQUILIBRIUM_POTENTIAL,
        )
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

    # --------------------------------------------------------
    # U = 0 V
    # --------------------------------------------------------

    draw_profile(
        ax,
        energies_0,
        color="black",
        connector_style="--",
    )

    # --------------------------------------------------------
    # U = 1.23 V
    # --------------------------------------------------------

    draw_profile(
        ax,
        energies_123,
        color="red",
        connector_style=":",
        connector_width=1.8,
    )

    # --------------------------------------------------------
    # U = 0 V labels
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
    ]

    ax.legend(
        handles=legend_handles,
        loc="upper left",
        frameon=False,
        fontsize=10,
        handlelength=3.0,
    )

    # --------------------------------------------------------
    # Axes
    # --------------------------------------------------------

    ax.set_xlim(
        -0.55,
        4.55,
    )

    all_energies = (
        energies_0
        + energies_123
    )

    ax.set_ylim(
        min(all_energies) - 0.45,
        max(all_energies) + 0.55,
    )

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

    ax.set_ylabel(
        r"$\Delta G$ (eV)",
        fontsize=12,
    )

    ax.tick_params(
        axis="x",
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
        left=True,
        right=False,
        labelleft=True,
        labelright=False,
        direction="out",
        width=1.2,
        length=4,
        labelsize=10,
    )

    for spine in ax.spines.values():

        spine.set_linewidth(
            1.3
        )

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
    corrected,
    delta_g,
    results,
    filename="oer_dft_thermodynamics.csv",
):
    """
    Save corrected energies and OER thermodynamic results.
    """

    step_names = [
        "* + H2O -> *OH + H+ + e-",
        "*OH -> *O + H+ + e-",
        "*O + H2O -> *OOH + H+ + e-",
        "*OOH -> * + O2 + H+ + e-",
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
                "Corrected species",
                "Free_energy_eV",
            ]
        )

        for name, value in corrected.items():

            writer.writerow(
                [
                    name,
                    f"{value:.6f}",
                ]
            )

        writer.writerow([])

        writer.writerow(
            [
                "Step",
                "Reaction",
                "DeltaG_eV",
            ]
        )

        for i, (
            reaction,
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
                    reaction,
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
    corrected,
    delta_g,
    results,
):
    """
    Print detailed OER thermodynamic results.
    """

    step_names = [
        "* + H2O -> *OH + H+ + e-",
        "*OH -> *O + H+ + e-",
        "*O + H2O -> *OOH + H+ + e-",
        "*OOH -> * + O2 + H+ + e-",
    ]

    print()

    print(
        "DFT-Based OER Thermodynamic Analysis"
    )

    print(
        "------------------------------------"
    )

    print()

    print(
        "Thermodynamic corrections:"
    )

    print(
        f"  H2   : {GCORR_H2:+.3f} eV"
    )

    print(
        f"  H2O  : {GCORR_H2O:+.3f} eV"
    )

    print(
        f"  OH*  : {GCORR_OH:+.3f} eV"
    )

    print(
        f"  O*   : {GCORR_O:+.3f} eV"
    )

    print(
        f"  OOH* : {GCORR_OOH:+.3f} eV"
    )

    print()

    print(
        "Corrected free energies:"
    )

    for name, value in corrected.items():

        print(
            f"  {name:<10s} = "
            f"{value:12.6f} eV"
        )

    print()

    print(
        "OER elementary steps:"
    )

    for i, (
        reaction,
        value,
    ) in enumerate(
        zip(
            step_names,
            delta_g,
        ),
        start=1,
    ):

        print(
            f"Step {i}: "
            f"{reaction}"
        )

        print(
            f"        Delta G = "
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
        "oer_dft_thermodynamics.png"
    )

    print(
        "Results saved to: "
        "oer_dft_thermodynamics.csv"
    )


# ============================================================
# COMMAND-LINE INPUT
# ============================================================

def parse_arguments():
    """
    Read DFT total energies from the command line.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Calculate conventional OER thermodynamics "
            "directly from DFT total energies."
        )
    )

    parser.add_argument(
        "e_clean",
        type=float,
        help="DFT total energy of clean surface * in eV",
    )

    parser.add_argument(
        "e_oh",
        type=float,
        help="DFT total energy of *OH surface in eV",
    )

    parser.add_argument(
        "e_o",
        type=float,
        help="DFT total energy of *O surface in eV",
    )

    parser.add_argument(
        "e_ooh",
        type=float,
        help="DFT total energy of *OOH surface in eV",
    )

    parser.add_argument(
        "e_h2",
        type=float,
        help="DFT total energy of H2 in eV",
    )

    parser.add_argument(
        "e_h2o",
        type=float,
        help="DFT total energy of H2O in eV",
    )

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():
    """
    Main program.
    """

    args = parse_arguments()

    # --------------------------------------------------------
    # Correct DFT energies
    # --------------------------------------------------------

    corrected = (
        calculate_corrected_energies(
            args.e_clean,
            args.e_oh,
            args.e_o,
            args.e_ooh,
            args.e_h2,
            args.e_h2o,
        )
    )

    # --------------------------------------------------------
    # Calculate OER step energies
    # --------------------------------------------------------

    delta_g = (
        calculate_oer_steps(
            corrected
        )
    )

    # --------------------------------------------------------
    # Thermodynamic metrics
    # --------------------------------------------------------

    results = (
        calculate_oer_metrics(
            delta_g
        )
    )

    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    plot_oer_diagram(
        delta_g,
        "oer_dft_thermodynamics.png",
    )

    # --------------------------------------------------------
    # CSV
    # --------------------------------------------------------

    save_csv(
        corrected,
        delta_g,
        results,
        "oer_dft_thermodynamics.csv",
    )

    # --------------------------------------------------------
    # Terminal output
    # --------------------------------------------------------

    print_results(
        corrected,
        delta_g,
        results,
    )


if __name__ == "__main__":
    main()
