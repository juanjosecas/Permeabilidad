# Membrane Permeability Analysis using Molecular Dynamics

This repository contains a Python script for analyzing membrane permeability through molecular dynamics simulations. The script utilizes the MDAnalysis library to process trajectory data and calculate the permeation of water molecules through a biological membrane. The analysis provides insights into how water molecules traverse the membrane over time.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [Important Considerations](#important-considerations)
- [Contact](#contact)

## Introduction

Understanding the permeation of molecules through biological membranes is crucial for various fields, including drug development and biomolecular research. This script aims to analyze the permeability of water molecules across a membrane using molecular dynamics simulations. The analysis involves tracking the movement of water molecules and detecting their passage through the membrane.

## Requirements

- Python (3.6 or higher)
- MDAnalysis library
- NumPy library

You can install the required libraries using the following commands:

```bash
pip install MDAnalysis numpy
```

## Usage

1. Clone this repository to your local machine:

```bash
git clone https://github.com/juanjosecas/Permeabilidad.git
```

2. Navigate to the repository directory:

```bash
cd membrane-permeability-analysis
```

3. Edit the script `permeability_analysis.py` to specify the paths to your topology (GRO) and trajectory (XTC) files. You can also adjust other parameters according to your simulation setup.

4. Run the script:

```bash
python permeability_analysis.py
```

5. After the analysis is complete, the results will be saved in the `permeation_results.csv` file, and a log of the analysis will be saved in the `permeation_log.txt` file.

## Important Considerations

- **Simulation Data**: Make sure you have valid topology and trajectory files (GRO and XTC) for your molecular dynamics simulation.
  
- **Python Knowledge**: Familiarity with Python programming and molecular dynamics concepts is essential to understand and modify the script as needed.

- **Memory Usage**: Depending on the size of your simulation data, the script may require significant memory. Ensure your system has sufficient resources to handle the analysis.

- **Contact**: For any questions or assistance, feel free to reach out to [jcasal@fmed.uba.ar](mailto:jcasal@fmed.uba.ar).

## Contact

For any inquiries or assistance, please contact:

**J. Casal**
Email: [jcasal@fmed.uba.ar](mailto:jcasal@fmed.uba.ar)
