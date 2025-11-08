# CATH Mainly Alpha - Protein Structure Analysis

Automated pipeline for downloading, cleaning, and analyzing mainly alpha-helix protein structures from the CATH database.

## Table of Contents

- [Overview](#overview)
- [Notebook Structure](#notebook-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output Structure](#output-structure)
- [Section Details](#section-details)
- [Main Functions](#main-functions)
- [Contributing](#contributing)
- [License](#license)

## Overview

This project implements a complete pipeline for analyzing proteins with mainly alpha-helix structure:

- Automated download of PDB structures from CATH database
- Cleaning and preprocessing of PDB files
- Amino acid frequency analysis
- Advanced analysis with DSSP (Define Secondary Structure of Proteins)
- Alpha-helix type classification
- Publication-quality visualizations

### Project Statistics

- 6 main processing cells
- ~1,962 lines of Python code
- 42 specialized functions
- 9 visualization functions
- Optimized for macOS with parallel processing

## Notebook Structure

The `cath-protocol.ipynb` notebook is organized into 6 main cells:

### Cell 0: Environment Setup (50 lines)
- Library imports
- Global variable configuration
- Initial environment setup

### Cell 1: PDB Structure Download (279 lines, 7 functions)
- CATH domain downloading
- Parallel download processing
- File validation

### Cell 2: Structure Cleaning (322 lines, 8 functions)
- Heteroatom removal
- Chain filtering
- PDB file standardization

### Cell 3: Amino Acid Frequency Analysis (156 lines, 3 functions)
- Amino acid counting
- Statistical analysis
- Distribution visualizations

### Cell 4: Advanced DSSP Analysis (558 lines, 13 functions)
- Secondary structure analysis
- Structural metrics calculation
- Alpha-helix detection

### Cell 5: Helix Type Classification (597 lines, 11 functions)
- Alpha-helix type identification
- Automatic classification
- Comparative analysis

## Requirements

### Python Dependencies

```bash
# Scientific core
numpy
pandas
matplotlib
seaborn

# Bioinformatics
biopython
biotite

# Structural analysis
dssp (mkdssp)

# Utilities
requests
tqdm
```

### External Tools

- **DSSP**: For secondary structure analysis
  ```bash
  # macOS
  brew install dssp

  # Linux (Ubuntu/Debian)
  sudo apt-get install dssp
  ```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/madsondeluna/cath-mainly-alpha.git
cd cath-mainly-alpha
```

2. Create a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install numpy pandas matplotlib seaborn biopython biotite requests tqdm
```

4. Install DSSP (if not already installed):
```bash
brew install dssp  # macOS
```

## Usage

### Basic Execution

Open the notebook in Jupyter:

```bash
jupyter notebook cath-protocol.ipynb
```

Execute cells in sequential order (0 → 1 → 2 → 3 → 4 → 5).

### Configuration

Edit variables in **Cell 0** to customize:

```python
# Output directories
OUTPUT_DIR = "output"
PDB_DIR = "pdb_files"
CLEANED_DIR = "cleaned_pdb"

# Download parameters
MAX_WORKERS = 4  # Number of parallel downloads
TIMEOUT = 30     # Timeout in seconds

# Analysis filters
MIN_HELIX_LENGTH = 4  # Minimum helix length
```

## Output Structure

After execution, the following directory structure will be created:

```
cath-mainly-alpha/
├── pdb_files/              # Downloaded PDB files
│   ├── domain1.pdb
│   ├── domain2.pdb
│   └── ...
├── cleaned_pdb/            # Cleaned PDB files
│   ├── domain1_clean.pdb
│   ├── domain2_clean.pdb
│   └── ...
├── output/                 # Analysis results
│   ├── amino_acid_freq.csv
│   ├── helix_analysis.csv
│   ├── helix_types.csv
│   └── figures/           # Visualizations
│       ├── aa_distribution.png
│       ├── helix_length_dist.png
│       └── ...
└── logs/                  # Execution logs
    └── processing.log
```

## Section Details

### Cell 0: Environment Setup

**Purpose**: Prepare the execution environment

**What it does**:
- Imports all necessary libraries (NumPy, Pandas, BioPython, etc.)
- Defines global variables and directory paths
- Configures matplotlib for high-quality visualizations
- Initializes loggers for tracking

**Main imports**:
```python
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
```

---

### Cell 1: PDB Structure Download

**Purpose**: Download PDB files from the CATH database

**Main functions**:

1. **`download_cath_domain(domain_id, output_dir)`**
   - Downloads a single CATH domain
   - Parameters: domain ID, output directory
   - Returns: Path to downloaded file or None

2. **`download_cath_domains_parallel(domain_list, output_dir, max_workers=4)`**
   - Parallel download of multiple domains
   - Uses ThreadPoolExecutor for parallelization
   - Progress bar with tqdm

3. **`validate_pdb_file(filepath)`**
   - Validates PDB file integrity
   - Checks format and structure
   - Returns: True/False

4. **`get_cath_domain_list(cath_version='current')`**
   - Gets list of CATH domains
   - Filters only "Mainly Alpha" class
   - Returns: List of IDs

**Process**:
1. Gets list of mainly alpha-helix domains from CATH
2. Creates output directory if it doesn't exist
3. Downloads PDB files in parallel (4 threads by default)
4. Validates each downloaded file
5. Records download statistics (success/failure)

**Output**:
- PDB files in `pdb_files/`
- Download log in `logs/download.log`

---

### Cell 2: Structure Cleaning

**Purpose**: Clean and standardize PDB files

**Main functions**:

1. **`remove_heteroatoms(structure)`**
   - Removes heteroatoms (water, ligands)
   - Keeps only protein atoms
   - Returns: Clean structure

2. **`filter_chains(structure, chain_ids=None)`**
   - Filters specific chains
   - If None, keeps all
   - Returns: Filtered structure

3. **`standardize_residues(structure)`**
   - Standardizes residue names
   - Removes non-standard residues
   - Returns: Standardized structure

4. **`clean_pdb_file(input_path, output_path)`**
   - Main cleaning function
   - Applies all filters
   - Saves clean file

5. **`batch_clean_pdb(input_dir, output_dir)`**
   - Cleans multiple files
   - Batch processing
   - Progress bar

**Process**:
1. Reads original PDB file
2. Removes heteroatoms (water, ions, ligands)
3. Filters chains (if specified)
4. Standardizes residue names
5. Removes modified residues
6. Saves clean file
7. Validates output file

**Cleaning criteria**:
- Removes HETATM
- Keeps only standard amino acids (20 types)
- Removes residues with occupancy < 0.5
- Removes atoms with B-factor > 100

**Output**:
- Clean PDB files in `cleaned_pdb/`
- Cleaning report in `logs/cleaning.log`

---

### Cell 3: Amino Acid Frequency Analysis

**Purpose**: Analyze amino acid composition in structures

**Main functions**:

1. **`count_amino_acids(pdb_file)`**
   - Counts frequency of each amino acid
   - Returns: Dictionary {aa: count}

2. **`calculate_aa_statistics(pdb_dir)`**
   - Aggregate statistics from all files
   - Mean, standard deviation, percentiles
   - Returns: DataFrame with statistics

3. **`plot_aa_distribution(stats_df, output_path)`**
   - Bar plot of distribution
   - Custom colors by property (hydrophobic, polar, etc.)
   - Saves high-resolution figure

**Analyses performed**:
- Absolute frequency of each amino acid
- Relative frequency (%)
- Comparison between structures
- Identification of most/least common amino acids
- Analysis by physicochemical property:
  - Hydrophobic (A, V, L, I, M, F, W, P)
  - Polar (S, T, N, Q, C, Y)
  - Positively charged (K, R, H)
  - Negatively charged (D, E)
  - Special (G, P)

**Visualizations**:
- Bar plot: general distribution
- Heatmap: comparison between structures
- Box plot: variability

**Output**:
- `output/amino_acid_freq.csv`: Frequency table
- `output/aa_statistics.csv`: Aggregate statistics
- `output/figures/aa_distribution.png`: Main plot
- `output/figures/aa_heatmap.png`: Comparative heatmap

---

### Cell 4: Advanced DSSP Analysis

**Purpose**: Detailed secondary structure analysis using DSSP

**Main functions**:

1. **`run_dssp(pdb_file)`**
   - Runs DSSP on PDB file
   - Returns: DSSP object with annotations

2. **`parse_dssp_output(dssp_result)`**
   - Parses DSSP output
   - Extracts secondary structure
   - Returns: Structured DataFrame

3. **`identify_helices(dssp_df)`**
   - Identifies alpha-helix segments
   - Filters by minimum length
   - Returns: List of helices [(start, end, length)]

4. **`calculate_helix_metrics(helix_segment, structure)`**
   - Calculates geometric metrics
   - Length, angles, twist
   - Returns: Dict with metrics

5. **`analyze_helix_geometry(pdb_file)`**
   - Complete geometric analysis
   - Helix axis, radius, pitch
   - Returns: DataFrame with geometry

**Calculated metrics**:

For each detected alpha-helix:
- **Length**: Number of residues
- **Phi/psi angles**: Dihedral angles
- **Rise per residue**: Advance per residue (ideal 3.6 Å)
- **Twist**: Rotation per residue (~100° ideal)
- **Radius**: Helix radius (~2.3 Å ideal)
- **Pitch**: Helix pitch (~5.4 Å ideal)
- **RMSD**: Deviation from ideal helix
- **Regularity**: Structural regularity score

**Identified secondary structures**:
- H: alpha-helix
- G: 3₁₀-helix
- I: pi-helix
- E: beta-sheet
- B: beta-bridge
- T: Turn
- S: Bend
- C: Coil/loop

**Visualizations**:
- Helix length distribution
- Ramachandran plots (phi/psi)
- Geometric metrics distribution
- Comparison with ideal values

**Output**:
- `output/dssp_analysis.csv`: Complete DSSP analysis
- `output/helix_metrics.csv`: Metrics for each helix
- `output/helix_geometry.csv`: Detailed geometry
- `output/figures/helix_length_dist.png`: Length distribution
- `output/figures/ramachandran.png`: Ramachandran plot
- `output/figures/geometry_metrics.png`: Geometric metrics

---

### Cell 5: Helix Type Classification

**Purpose**: Classify alpha-helices into specific types

**Main functions**:

1. **`classify_helix_type(helix_metrics)`**
   - Classifies helix type based on metrics
   - alpha-helix, 3₁₀-helix, pi-helix, irregular
   - Returns: Helix type

2. **`identify_helix_kinks(helix_segment, threshold=20)`**
   - Detects kinks (bends) in helices
   - Threshold: deviation angle in degrees
   - Returns: List of kink positions

3. **`calculate_helix_stability(helix_segment)`**
   - Estimates helix stability
   - Based on hydrogen bonds
   - Returns: Stability score (0-1)

4. **`find_helix_capping_residues(helix_segment)`**
   - Identifies capping residues (N-cap and C-cap)
   - Important for stability
   - Returns: {n_cap: residue, c_cap: residue}

5. **`analyze_helix_surface(helix_segment)`**
   - Analyzes solvent exposure
   - Identifies hydrophobic/hydrophilic faces
   - Returns: Dict with surface analysis

6. **`detect_helix_helix_interactions(structure)`**
   - Detects interactions between helices
   - Packing, crossing angles
   - Returns: List of interacting pairs

**Classified helix types**:

1. **Canonical alpha-helix**
   - Rise: ~1.5 Å/residue
   - Twist: ~100°/residue
   - 3.6 residues/turn

2. **3₁₀-helix**
   - Rise: ~2.0 Å/residue
   - Twist: ~120°/residue
   - 3.0 residues/turn
   - Tighter

3. **Pi-helix**
   - Rise: ~1.2 Å/residue
   - Twist: ~87°/residue
   - 4.4 residues/turn
   - Wider

4. **Irregular helix**
   - Significant deviations from patterns
   - May contain kinks
   - RMSD > 1.0 Å

**Special analyses**:

- **Helix dipole**: Helix dipole moment
- **Capping motifs**: Terminal residue patterns
- **Hydrophobic moment**: Hydrophobic moment
- **Helix-helix packing**: Packing geometry

**Visualizations**:
- Helix type distribution
- Kink map along sequence
- Capping residue analysis
- Solvent exposure profile
- Helix-helix interaction network

**Output**:
- `output/helix_classification.csv`: Classification of each helix
- `output/helix_types_summary.csv`: Summary by type
- `output/kinks_analysis.csv`: Kink analysis
- `output/capping_residues.csv`: Capping residues
- `output/helix_interactions.csv`: Helix interactions
- `output/figures/helix_types_pie.png`: Type distribution
- `output/figures/kinks_heatmap.png`: Kink map
- `output/figures/interaction_network.png`: Interaction network

---

## Main Functions

### By Category

#### Download and I/O
- `download_cath_domain()`: Individual download
- `download_cath_domains_parallel()`: Parallel download
- `validate_pdb_file()`: File validation
- `get_cath_domain_list()`: Domain list

#### Cleaning
- `remove_heteroatoms()`: Remove non-protein atoms
- `filter_chains()`: Filter chains
- `standardize_residues()`: Standardize residues
- `clean_pdb_file()`: Complete cleaning

#### Sequence Analysis
- `count_amino_acids()`: Count amino acids
- `calculate_aa_statistics()`: Aggregate statistics

#### Structural Analysis (DSSP)
- `run_dssp()`: Run DSSP
- `parse_dssp_output()`: Parse results
- `identify_helices()`: Identify helices
- `calculate_helix_metrics()`: Geometric metrics
- `analyze_helix_geometry()`: Detailed geometry

#### Classification
- `classify_helix_type()`: Classify types
- `identify_helix_kinks()`: Detect kinks
- `calculate_helix_stability()`: Stability
- `find_helix_capping_residues()`: Capping
- `detect_helix_helix_interactions()`: Interactions

#### Visualization
- `plot_aa_distribution()`: Amino acid distribution
- `plot_helix_length_distribution()`: Lengths
- `plot_ramachandran()`: Ramachandran plot
- `plot_helix_types()`: Helix types
- `plot_interaction_network()`: Interaction network

---

## Scientific Applications

This pipeline can be used for:

1. **Structural research**: Characterization of alpha-helices in different contexts
2. **Comparative analysis**: Compare properties between protein families
3. **Structure prediction**: Validate computational models
4. **Protein design**: Inform rational design
5. **Education**: Demonstrate protein structure concepts

---

## Performance

- **Download**: ~100 structures/minute (4 threads)
- **Cleaning**: ~50 structures/minute
- **DSSP analysis**: ~20 structures/minute
- **Classification**: ~30 structures/minute

For large datasets (>1000 structures), it is recommended to:
- Increase `MAX_WORKERS` to 8-16 (if you have enough CPU)
- Run on machine with SSD
- Use at least 8GB RAM

---

## Contributing

Contributions are welcome! Please:

1. Fork the project
2. Create a feature branch (`git checkout -b feature/NewAnalysis`)
3. Commit your changes (`git commit -m 'Add new analysis of X'`)
4. Push to the branch (`git push origin feature/NewAnalysis`)
5. Open a Pull Request

---

## License

This project is under the MIT License. See the `LICENSE` file for more details.

---

## Author

**Madson de Luna**
- GitHub: [@madsondeluna](https://github.com/madsondeluna)

---

## References

- CATH Database: http://www.cathdb.info/
- DSSP: Kabsch W, Sander C (1983). "Dictionary of protein secondary structure"
- BioPython: Cock et al. (2009). "Biopython: freely available Python tools"

---

## Updates

### Current Version
- Complete functional pipeline
- 42 implemented functions
- 9 visualization types
- Support for parallel processing

### Planned Features
- AlphaFold integration
- Molecular dynamics analysis
- Web interface
- REST API
- Export to additional formats (PyMOL, Chimera)
