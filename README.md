# Computational Electrocatalysis Tools

Python workflows for computational electrocatalysis, including computational hydrogen electrode (CHE) analysis, oxygen evolution reaction (OER) thermodynamics, limiting-potential calculations, and potential-dependent surface phase stability.

The repository is intended to provide lightweight and reusable tools for converting first-principles energetics into electrochemical thermodynamic quantities and publication-ready visualizations.

## Features

- OER free-energy diagrams from reaction-step free energies
- Potential-dependent OER free-energy profiles
- Limiting-potential and overpotential analysis
- OER thermodynamics directly from DFT total energies
- Potential-dependent OH/O surface phase diagrams
- Automatic identification of thermodynamically stable surface phases
- Analytical determination of surface-phase transition potentials
- CSV export of calculated thermodynamic data
- Publication-quality PNG and PDF figures

## Repository Structure

```text
computational-electrocatalysis-tools/
├── scripts/
│   ├── oer_free_energy.py
│   ├── oer_limiting_potential.py
│   ├── oer_dft_thermodynamics.py
│   └── surface_phase_diagram.py
├── examples/
│   └── surface_energies_example.csv
├── README.md
├── requirements.txt
├── LICENSE
└── .gitignore
```

## Installation

Clone the repository:

```bash
git clone https://github.com/Radha1411/computational-electrocatalysis-tools.git
cd computational-electrocatalysis-tools
```

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Current dependencies:

- NumPy
- Matplotlib

---

## 1. OER Free-Energy Analysis

### `oer_free_energy.py`

Generates a conventional four-step OER free-energy diagram from user-supplied reaction free energies.

For the four OER steps,

```text
* + H2O  -> *OH + H+ + e-
*OH      -> *O  + H+ + e-
*O + H2O -> *OOH + H+ + e-
*OOH     -> * + O2 + H+ + e-
```

the script accepts the four reaction free energies at \(U=0\) V.

### Usage

```bash
python3 scripts/oer_free_energy.py 1.10 1.45 1.80 0.57
```

The four numerical arguments correspond to:

```text
ΔG1  ΔG2  ΔG3  ΔG4
```

The script determines:

- potential-determining step (PDS)
- maximum uphill free-energy step
- limiting potential
- theoretical OER overpotential

For a conventional four-electron OER pathway,

```text
U_L = max(ΔG1, ΔG2, ΔG3, ΔG4)
```

and

```text
η_OER = U_L - 1.23 V
```

### Outputs

```text
oer_free_energy.png
oer_summary.csv
```

---

## 2. OER Limiting-Potential Analysis

### `oer_limiting_potential.py`

Extends the conventional OER free-energy analysis by explicitly comparing the reaction profile at:

- \(U=0\) V
- \(U=1.23\) V
- the calculated limiting potential \(U_L\)

### Usage

```bash
python3 scripts/oer_limiting_potential.py 1.10 1.45 1.80 0.57
```

The limiting potential corresponds to the potential at which the largest uphill one-electron OER step becomes thermoneutral.

### Output

```text
oer_limiting_potential.png
```

---

## 3. DFT-Based OER Thermodynamics

### `oer_dft_thermodynamics.py`

Calculates conventional OER thermodynamics directly from DFT total energies of the clean surface and the key adsorbed intermediates:

```text
*
*OH
*O
*OOH
```

together with calculated gas-phase energies of H2 and H2O.

The script applies thermodynamic corrections and uses the computational hydrogen electrode framework.

The first three OER reaction free energies are evaluated from the corrected DFT energies. The fourth step is obtained using the total four-electron OER free-energy requirement:

```text
ΔG1 + ΔG2 + ΔG3 + ΔG4 = 4.92 eV
```

This avoids relying directly on the DFT total energy of gas-phase O2.

### Usage

```bash
python3 scripts/oer_dft_thermodynamics.py \
E_clean E_OH E_O E_OOH E_H2 E_H2O
```

The arguments are supplied in the following order:

```text
1. E_clean
2. E_OH
3. E_O
4. E_OOH
5. E_H2
6. E_H2O
```

All electronic energies should be supplied in eV and obtained using a consistent computational setup.

### Outputs

```text
oer_dft_thermodynamics.png
oer_dft_thermodynamics.csv
```

The script reports:

- four OER reaction free energies
- potential-determining step
- limiting potential
- OER overpotential

---

## 4. Potential-Dependent Surface Phase Diagram

### `surface_phase_diagram.py`

Analyzes the potential-dependent thermodynamic stability of catalytic surface configurations containing combinations of OH and O species.

The current implementation considers four explicitly defined surface sites, with each site represented as either:

```text
OH
```

or

```text
O
```

Possible configurations therefore include, for example:

```text
OH_OH_OH_OH
OH_O_OH_OH
O_O_OH_OH
O_O_O_OH
O_O_O_O
```

Different arrangements at the same overall OH/O coverage can be supplied independently. This allows inequivalent active sites and different adsorbate arrangements to be compared explicitly.

The surface free energies are evaluated as a function of applied potential according to the thermodynamic formulation implemented in the script:

```text
ΔG(U) = ΔG(0) - nU
```

where `n` is determined from the surface configuration.

No OH/O surface configuration is shifted to zero. All supplied surface phases are evaluated on the same thermodynamic scale.

The script then determines the configuration with the lowest free energy at every potential and constructs the corresponding stable-phase diagram.

### Input CSV

Surface total energies are supplied through a CSV file:

```csv
configuration,total_energy
OH_OH_OH_OH,-120.000
OH_O_OH_OH,-116.500
O_OH_OH_OH,-116.420
O_O_OH_OH,-112.900
O_O_O_OH,-109.100
O_O_O_O,-105.400
```

The values above are synthetic and are provided only to demonstrate the required input format.

A complete example input is available at:

```text
examples/surface_energies_example.csv
```

The user is not required to provide all possible configurations. Only explicitly calculated configurations are compared.

### Usage

```bash
python3 scripts/surface_phase_diagram.py \
examples/surface_energies_example.csv \
--pristine -60.000 \
--h2 -7.000 \
--h2o -13.500 \
--h2-correction -0.04 \
--h2o-correction -0.01 \
--umin 0.6 \
--umax 2.0 \
--target 1.23
```

All numerical values in this example are synthetic.

The main command-line options are:

```text
--pristine        DFT total energy of the pristine surface
--h2              DFT total energy of H2
--h2o             DFT total energy of H2O
--h2-correction   Free-energy correction for H2
--h2o-correction  Free-energy correction for H2O
--umin             Minimum applied potential
--umax             Maximum applied potential
--points           Number of potential-grid points
--target           Potential for reporting the stable phase
```

### Outputs

```text
surface_phase_diagram.png
surface_phase_diagram.pdf
surface_phase_data.csv
stable_phase_regions.csv
```

The workflow:

1. reads the calculated surface configurations,
2. counts the OH and O species,
3. constructs the potential-dependent surface free energies,
4. compares all supplied configurations,
5. identifies the thermodynamically stable phase,
6. determines phase-transition potentials analytically,
7. reports the stable phase at a selected potential,
8. exports the complete potential-dependent dataset, and
9. generates a surface phase diagram with a stable-phase map.

The phase-transition potential between two stable configurations is obtained from the intersection of their free-energy lines.

For

```text
ΔG1(U) = b1 - n1 U
ΔG2(U) = b2 - n2 U
```

the transition occurs at

```text
Utransition = (b2 - b1) / (n2 - n1)
```

when \(n_1 \ne n_2\).

---

## Computational Hydrogen Electrode

The workflows in this repository use the computational hydrogen electrode (CHE) framework to introduce electrochemical potential dependence into DFT-derived reaction and surface thermodynamics.

For a proton-electron transfer,

```text
H+ + e-
```

the electrochemical potential is related to the hydrogen reference through the CHE formalism.

For an electrochemical state involving `n` potential-dependent electron transfers, the corresponding free energy varies linearly with applied potential according to the reaction convention implemented in the relevant script.

This enables DFT energetics to be used for constructing OER free-energy diagrams and potential-dependent surface stability diagrams.

---

## Notes

- All energies should be supplied in eV.
- DFT energies entering the same thermodynamic analysis should be calculated using a consistent computational setup.
- Free-energy corrections should be chosen consistently with the thermodynamic model being used.
- The potential-determining step is a thermodynamic quantity and should not automatically be interpreted as the kinetic rate-determining step.
- Surface phase diagrams describe relative thermodynamic stability among the configurations included in the calculation; they do not by themselves establish kinetic accessibility.
- Example energies distributed with this repository are synthetic and are intended only for demonstrating the workflows.

## License

This project is licensed under the MIT License.

## Author

**Radha Somaiya**

Computational materials science and electrocatalysis  
DFT | Computational electrocatalysis | Electrochemical thermodynamics | Surface stability
