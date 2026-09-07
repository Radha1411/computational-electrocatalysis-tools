"""
surface_phase_diagram.py

Potential-dependent surface phase diagram for catalytic
surfaces containing combinations of OH and O species.

The script is designed for a surface with four active sites,
where each site can be represented as either:

    OH
    O

Examples:

    OH_OH_OH_OH
    OH_O_OH_OH
    O_O_OH_OH
    O_O_O_O

No OH/O surface configuration is shifted to zero.

The surface formation free energy is calculated as:

    DeltaG(U)
        = E_surface
        - E_pristine
        - 4 G(H2O)
        + n/2 G(H2)
        - n U

where:

    E_surface
        DFT total energy of a surface configuration

    E_pristine
        DFT total energy of the pristine surface

    G(H2)
        Gibbs free energy of H2

    G(H2O)
        Gibbs free energy of H2O

    n
        Number of O species in the surface configuration

    U
        Applied potential in V vs RHE

The thermodynamically stable surface phase is the configuration
with the lowest DeltaG at a given potential.

Input CSV format:

    configuration,total_energy
    OH_OH_OH_OH,-100.000
    OH_O_OH_OH,-95.000
    O_O_OH_OH,-90.000

Use synthetic/example values when sharing publicly.

Outputs:

    surface_phase_diagram.png
    surface_phase_diagram.pdf
    surface_phase_data.csv
    stable_phase_regions.csv

Author: Radha Somaiya
"""

import argparse
import csv

import matplotlib

# Non-interactive backend for HPC systems
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


# ============================================================
# DEFAULT SETTINGS
# ============================================================

DEFAULT_U_MIN = 0.60
DEFAULT_U_MAX = 2.00
DEFAULT_POINTS = 5000

DEFAULT_H2_CORRECTION = -0.04
DEFAULT_H2O_CORRECTION = -0.01

NUMBER_OF_SURFACE_SITES = 4


# ============================================================
# READ SURFACE ENERGIES
# ============================================================

def read_surface_energies(filename):
    """
    Read surface configurations and total energies from CSV.

    Expected format:

        configuration,total_energy
        OH_OH_OH_OH,-100.000
        OH_O_OH_OH,-95.000

    Returns
    -------
    dict
        configuration -> DFT total energy
    """

    energies = {}

    with open(filename, "r", newline="") as csvfile:

        reader = csv.DictReader(csvfile)

        required_columns = {
            "configuration",
            "total_energy",
        }

        if not required_columns.issubset(
            reader.fieldnames or []
        ):
            raise ValueError(
                "Input CSV must contain columns:\n"
                "configuration,total_energy"
            )

        for row in reader:

            configuration = (
                row["configuration"].strip()
            )

            if configuration in energies:
                raise ValueError(
                    f"Duplicate configuration found: "
                    f"{configuration}"
                )

            try:
                energy = float(
                    row["total_energy"]
                )

            except ValueError as error:
                raise ValueError(
                    f"Invalid total energy for "
                    f"{configuration}"
                ) from error

            validate_configuration(
                configuration
            )

            energies[
                configuration
            ] = energy

    if not energies:
        raise ValueError(
            "No surface configurations were found "
            "in the input file."
        )

    return energies


# ============================================================
# VALIDATE CONFIGURATION
# ============================================================

def validate_configuration(configuration):
    """
    Check that a configuration contains exactly four sites
    and that every site is either OH or O.
    """

    species = configuration.split("_")

    if len(species) != NUMBER_OF_SURFACE_SITES:

        raise ValueError(
            f"Configuration '{configuration}' "
            f"does not contain exactly "
            f"{NUMBER_OF_SURFACE_SITES} sites."
        )

    for item in species:

        if item not in ("OH", "O"):

            raise ValueError(
                f"Invalid species '{item}' in "
                f"configuration '{configuration}'. "
                "Allowed species are OH and O."
            )


# ============================================================
# COUNT SURFACE SPECIES
# ============================================================

def count_surface_species(configuration):
    """
    Count OH and O species.

    Example:

        OH_O_OH_OH

    gives:

        OH = 3
        O  = 1
    """

    species = configuration.split("_")

    number_oh = species.count("OH")
    number_o = species.count("O")

    return (
        number_oh,
        number_o,
    )


# ============================================================
# GAS-PHASE GIBBS ENERGIES
# ============================================================

def calculate_gas_free_energy(
    electronic_energy,
    correction,
):
    """
    Calculate Gibbs free energy:

        G = E_DFT + correction
    """

    return (
        electronic_energy
        + correction
    )


# ============================================================
# SURFACE FORMATION FREE ENERGY
# ============================================================

def calculate_surface_states(
    energies,
    e_pristine,
    g_h2,
    g_h2o,
    potentials,
):
    """
    Calculate DeltaG(U) for every surface configuration.

    Equation:

        DeltaG(U)
            = E_surface
            - E_pristine
            - 4 G(H2O)
            + n/2 G(H2)
            - n U

    where n is the number of O species.
    """

    states = {}
    intercepts = {}
    n_electrons = {}

    for configuration, e_surface in energies.items():

        number_oh, number_o = (
            count_surface_species(
                configuration
            )
        )

        n = number_o

        delta_g0 = (
            e_surface
            - e_pristine
            - NUMBER_OF_SURFACE_SITES * g_h2o
            + 0.5 * n * g_h2
        )

        delta_g_u = (
            delta_g0
            - n * potentials
        )

        states[
            configuration
        ] = delta_g_u

        intercepts[
            configuration
        ] = delta_g0

        n_electrons[
            configuration
        ] = n

    return (
        states,
        intercepts,
        n_electrons,
    )


# ============================================================
# FIND STABLE PHASE
# ============================================================

def determine_stable_phases(
    states,
):
    """
    Determine the lowest-free-energy phase
    at every potential.
    """

    state_names = list(
        states.keys()
    )

    all_g = np.vstack(
        [
            states[name]
            for name in state_names
        ]
    )

    stable_indices = np.argmin(
        all_g,
        axis=0,
    )

    stable_names = np.array(
        state_names
    )[stable_indices]

    minimum_energies = np.min(
        all_g,
        axis=0,
    )

    return (
        state_names,
        stable_names,
        minimum_energies,
    )


# ============================================================
# INITIAL STABLE REGIONS
# ============================================================

def find_stable_regions(
    potentials,
    stable_names,
):
    """
    Find approximate stable-phase regions
    using the potential grid.
    """

    stable_regions = []

    current_state = stable_names[0]

    start_index = 0

    for i in range(
        1,
        len(potentials),
    ):

        if (
            stable_names[i]
            != current_state
        ):

            stable_regions.append(
                {
                    "state": current_state,
                    "U_start": potentials[
                        start_index
                    ],
                    "U_end": potentials[
                        i - 1
                    ],
                }
            )

            current_state = (
                stable_names[i]
            )

            start_index = i

    stable_regions.append(
        {
            "state": current_state,
            "U_start": potentials[
                start_index
            ],
            "U_end": potentials[-1],
        }
    )

    return stable_regions


# ============================================================
# ANALYTICAL TRANSITION POTENTIALS
# ============================================================

def refine_transition_potentials(
    stable_regions,
    intercepts,
    n_electrons,
):
    """
    Refine transition potentials analytically.

    For two states:

        G1 = b1 - n1 U
        G2 = b2 - n2 U

    At crossover:

        b1 - n1 U
            =
        b2 - n2 U

    therefore:

        U =
        (b2 - b1)
        /
        (n2 - n1)
    """

    for i in range(
        len(stable_regions) - 1
    ):

        state_left = (
            stable_regions[i]["state"]
        )

        state_right = (
            stable_regions[i + 1]["state"]
        )

        b_left = intercepts[
            state_left
        ]

        b_right = intercepts[
            state_right
        ]

        n_left = n_electrons[
            state_left
        ]

        n_right = n_electrons[
            state_right
        ]

        denominator = (
            n_right
            - n_left
        )

        if denominator == 0:
            continue

        transition_u = (
            b_right
            - b_left
        ) / denominator

        stable_regions[i][
            "U_end"
        ] = transition_u

        stable_regions[i + 1][
            "U_start"
        ] = transition_u

    return stable_regions


# ============================================================
# REPORT TARGET POTENTIAL
# ============================================================

def stable_phase_at_potential(
    target_u,
    state_names,
    intercepts,
    n_electrons,
):
    """
    Determine thermodynamically stable phase
    at a selected potential.
    """

    free_energies = {}

    for name in state_names:

        free_energies[name] = (
            intercepts[name]
            - n_electrons[name]
            * target_u
        )

    stable_state = min(
        free_energies,
        key=free_energies.get,
    )

    return (
        stable_state,
        free_energies[
            stable_state
        ],
    )


# ============================================================
# PRINT RESULTS
# ============================================================

def print_results(
    energies,
    intercepts,
    n_electrons,
    stable_regions,
    target_u,
    target_state,
    target_energy,
    g_h2,
    g_h2o,
):
    """
    Print thermodynamic results.
    """

    print()

    print("=" * 80)
    print(
        "SURFACE PHASE DIAGRAM"
    )
    print("=" * 80)

    print()

    print(
        f"G(H2)  = {g_h2:.6f} eV"
    )

    print(
        f"G(H2O) = {g_h2o:.6f} eV"
    )

    print()

    print(
        "Surface free-energy equations"
    )

    print("-" * 80)

    for configuration in energies:

        number_oh, number_o = (
            count_surface_species(
                configuration
            )
        )

        intercept = intercepts[
            configuration
        ]

        n = n_electrons[
            configuration
        ]

        print(
            f"{configuration:18s} | "
            f"OH = {number_oh} | "
            f"O = {number_o} | "
            f"DeltaG(U) = "
            f"{intercept:8.3f} "
            f"- {n}U"
        )

    print()

    print("=" * 80)
    print(
        "THERMODYNAMICALLY STABLE SURFACE PHASES"
    )
    print("=" * 80)

    for region in stable_regions:

        print(
            f"{region['state']:18s}: "
            f"{region['U_start']:.3f} "
            f"to "
            f"{region['U_end']:.3f} V"
        )

    print()

    print("=" * 80)

    print(
        f"STABLE PHASE AT U = "
        f"{target_u:.2f} V"
    )

    print("=" * 80)

    print(
        f"Stable phase = "
        f"{target_state}"
    )

    print(
        f"Surface formation free energy = "
        f"{target_energy:.3f} eV"
    )


# ============================================================
# SAVE POTENTIAL-DEPENDENT DATA
# ============================================================

def save_phase_data(
    filename,
    potentials,
    state_names,
    states,
    stable_names,
    minimum_energies,
):
    """
    Save DeltaG(U) for all configurations.
    """

    with open(
        filename,
        "w",
        newline="",
    ) as csvfile:

        writer = csv.writer(
            csvfile
        )

        header = [
            "Potential_V"
        ]

        header.extend(
            [
                f"DG_{name}_eV"
                for name in state_names
            ]
        )

        header.extend(
            [
                "Stable_phase",
                "Minimum_DG_eV",
            ]
        )

        writer.writerow(
            header
        )

        for i, potential in enumerate(
            potentials
        ):

            row = [
                f"{potential:.6f}"
            ]

            for name in state_names:

                row.append(
                    f"{states[name][i]:.6f}"
                )

            row.extend(
                [
                    stable_names[i],
                    f"{minimum_energies[i]:.6f}",
                ]
            )

            writer.writerow(
                row
            )


# ============================================================
# SAVE STABLE REGIONS
# ============================================================

def save_stable_regions(
    filename,
    stable_regions,
):
    """
    Save phase-transition regions.
    """

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
                "Stable_phase",
                "U_start_V",
                "U_end_V",
            ]
        )

        for region in stable_regions:

            writer.writerow(
                [
                    region["state"],
                    f"{region['U_start']:.6f}",
                    f"{region['U_end']:.6f}",
                ]
            )


# ============================================================
# PLOT PHASE DIAGRAM
# ============================================================

def plot_phase_diagram(
    potentials,
    states,
    stable_regions,
    output_png,
    output_pdf,
):
    """
    Generate free-energy diagram and stable-phase strip.
    """

    fig = plt.figure(
        figsize=(13.5, 8.5)
    )

    grid = GridSpec(
        nrows=2,
        ncols=1,
        height_ratios=[
            5.0,
            1.0,
        ],
        hspace=0.16,
    )

    ax = fig.add_subplot(
        grid[0]
    )

    ax_phase = fig.add_subplot(
        grid[1],
        sharex=ax,
    )

    # --------------------------------------------------------
    # Plot all surface configurations
    # --------------------------------------------------------

    line_colors = {}

    for state_name, free_energy in states.items():

        line = ax.plot(
            potentials,
            free_energy,
            linewidth=1.7,
            alpha=0.60,
            label=state_name,
        )[0]

        line_colors[
            state_name
        ] = line.get_color()

    # --------------------------------------------------------
    # Highlight stable line segments
    # --------------------------------------------------------

    for region in stable_regions:

        state = region[
            "state"
        ]

        u_start = region[
            "U_start"
        ]

        u_end = region[
            "U_end"
        ]

        mask = (
            (potentials >= u_start)
            &
            (potentials <= u_end)
        )

        ax.plot(
            potentials[mask],
            states[state][mask],
            linewidth=4.0,
            color=line_colors[state],
            solid_capstyle="round",
            zorder=20,
        )

    # --------------------------------------------------------
    # Transition lines
    # --------------------------------------------------------

    transition_potentials = [
        region["U_end"]
        for region
        in stable_regions[:-1]
    ]

    for transition_u in transition_potentials:

        ax.axvline(
            transition_u,
            linestyle="--",
            linewidth=1.5,
        )

        ax_phase.axvline(
            transition_u,
            linestyle="--",
            linewidth=1.3,
        )

        ax.text(
            transition_u + 0.01,
            0.96,
            f"{transition_u:.2f} V",
            transform=(
                ax.get_xaxis_transform()
            ),
            fontsize=12,
            verticalalignment="top",
        )

    # --------------------------------------------------------
    # Upper panel formatting
    # --------------------------------------------------------

    ax.set_ylabel(
        r"$\Delta G$ (eV)",
        fontsize=16,
    )

    ax.set_xlim(
        potentials[0],
        potentials[-1],
    )

    ax.tick_params(
        axis="both",
        direction="in",
        top=True,
        right=True,
        width=1.5,
        length=5,
        labelsize=11,
    )

    ax.tick_params(
        labelbottom=False
    )

    for spine in ax.spines.values():

        spine.set_linewidth(
            1.5
        )

    ax.legend(
        ncol=1,
        loc="center left",
        bbox_to_anchor=(
            1.01,
            0.50,
        ),
        frameon=False,
        fontsize=8.5,
    )

    # --------------------------------------------------------
    # Stable-phase strip
    # --------------------------------------------------------

    for region in stable_regions:

        state = region[
            "state"
        ]

        u_start = region[
            "U_start"
        ]

        u_end = region[
            "U_end"
        ]

        width = (
            u_end
            - u_start
        )

        rectangle = Rectangle(
            (
                u_start,
                0.0,
            ),
            width,
            1.0,
            facecolor=(
                line_colors[state]
            ),
            edgecolor="black",
            linewidth=1.0,
        )

        ax_phase.add_patch(
            rectangle
        )

        ax_phase.text(
            u_start
            + width / 2.0,
            0.50,
            state,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=10,
            color="white",
        )

    # --------------------------------------------------------
    # Lower panel formatting
    # --------------------------------------------------------

    ax_phase.set_ylim(
        0.0,
        1.0,
    )

    ax_phase.set_yticks(
        []
    )

    ax_phase.set_xlabel(
        "Potential (V vs. RHE)",
        fontsize=15,
        labelpad=10,
    )

    ax_phase.tick_params(
        axis="x",
        direction="out",
        width=1.5,
        length=5,
        labelsize=11,
    )

    for spine in ax_phase.spines.values():

        spine.set_linewidth(
            1.5
        )

    ax_phase.text(
        potentials[0],
        1.15,
        "Stable phase",
        fontsize=12,
        horizontalalignment="left",
        verticalalignment="bottom",
        clip_on=False,
    )

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    plt.subplots_adjust(
        left=0.10,
        right=0.76,
        top=0.96,
        bottom=0.12,
    )

    plt.savefig(
        output_png,
        dpi=600,
        bbox_inches="tight",
    )

    plt.savefig(
        output_pdf,
        bbox_inches="tight",
    )

    plt.close(fig)


# ============================================================
# COMMAND-LINE INPUT
# ============================================================

def parse_arguments():
    """
    Parse command-line options.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Generate a potential-dependent OH/O "
            "surface phase diagram."
        )
    )

    parser.add_argument(
        "input_csv",
        help=(
            "CSV file containing configuration "
            "and total_energy columns"
        ),
    )

    parser.add_argument(
        "--pristine",
        type=float,
        required=True,
        help=(
            "DFT total energy of pristine surface "
            "in eV"
        ),
    )

    parser.add_argument(
        "--h2",
        type=float,
        required=True,
        help=(
            "DFT total energy of H2 in eV"
        ),
    )

    parser.add_argument(
        "--h2o",
        type=float,
        required=True,
        help=(
            "DFT total energy of H2O in eV"
        ),
    )

    parser.add_argument(
        "--h2-correction",
        type=float,
        default=DEFAULT_H2_CORRECTION,
        help=(
            "Free-energy correction for H2 "
            f"(default: "
            f"{DEFAULT_H2_CORRECTION} eV)"
        ),
    )

    parser.add_argument(
        "--h2o-correction",
        type=float,
        default=DEFAULT_H2O_CORRECTION,
        help=(
            "Free-energy correction for H2O "
            f"(default: "
            f"{DEFAULT_H2O_CORRECTION} eV)"
        ),
    )

    parser.add_argument(
        "--umin",
        type=float,
        default=DEFAULT_U_MIN,
        help=(
            "Minimum potential in V vs RHE"
        ),
    )

    parser.add_argument(
        "--umax",
        type=float,
        default=DEFAULT_U_MAX,
        help=(
            "Maximum potential in V vs RHE"
        ),
    )

    parser.add_argument(
        "--points",
        type=int,
        default=DEFAULT_POINTS,
        help=(
            "Number of potential grid points"
        ),
    )

    parser.add_argument(
        "--target",
        type=float,
        default=1.23,
        help=(
            "Potential at which stable phase "
            "is reported in terminal"
        ),
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

    if args.umax <= args.umin:

        raise ValueError(
            "Umax must be greater than Umin."
        )

    if args.points < 2:

        raise ValueError(
            "At least two potential points "
            "are required."
        )

    # --------------------------------------------------------
    # Read surface configurations
    # --------------------------------------------------------

    energies = read_surface_energies(
        args.input_csv
    )

    # --------------------------------------------------------
    # Gas-phase Gibbs energies
    # --------------------------------------------------------

    g_h2 = calculate_gas_free_energy(
        args.h2,
        args.h2_correction,
    )

    g_h2o = calculate_gas_free_energy(
        args.h2o,
        args.h2o_correction,
    )

    # --------------------------------------------------------
    # Potential range
    # --------------------------------------------------------

    potentials = np.linspace(
        args.umin,
        args.umax,
        args.points,
    )

    # --------------------------------------------------------
    # Calculate surface states
    # --------------------------------------------------------

    (
        states,
        intercepts,
        n_electrons,
    ) = calculate_surface_states(
        energies,
        args.pristine,
        g_h2,
        g_h2o,
        potentials,
    )

    # --------------------------------------------------------
    # Stable phases
    # --------------------------------------------------------

    (
        state_names,
        stable_names,
        minimum_energies,
    ) = determine_stable_phases(
        states
    )

    stable_regions = (
        find_stable_regions(
            potentials,
            stable_names,
        )
    )

    stable_regions = (
        refine_transition_potentials(
            stable_regions,
            intercepts,
            n_electrons,
        )
    )

    # --------------------------------------------------------
    # Stable phase at selected potential
    # --------------------------------------------------------

    (
        target_state,
        target_energy,
    ) = stable_phase_at_potential(
        args.target,
        state_names,
        intercepts,
        n_electrons,
    )

    # --------------------------------------------------------
    # Print
    # --------------------------------------------------------

    print_results(
        energies,
        intercepts,
        n_electrons,
        stable_regions,
        args.target,
        target_state,
        target_energy,
        g_h2,
        g_h2o,
    )

    # --------------------------------------------------------
    # Save CSV files
    # --------------------------------------------------------

    save_phase_data(
        "surface_phase_data.csv",
        potentials,
        state_names,
        states,
        stable_names,
        minimum_energies,
    )

    save_stable_regions(
        "stable_phase_regions.csv",
        stable_regions,
    )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plot_phase_diagram(
        potentials,
        states,
        stable_regions,
        "surface_phase_diagram.png",
        "surface_phase_diagram.pdf",
    )

    print()

    print(
        "Generated files:"
    )

    print(
        "  surface_phase_diagram.png"
    )

    print(
        "  surface_phase_diagram.pdf"
    )

    print(
        "  surface_phase_data.csv"
    )

    print(
        "  stable_phase_regions.csv"
    )


if __name__ == "__main__":
    main()
