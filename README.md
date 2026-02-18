# UF-LGP: University of Florida - Livestock Genomics Pipelines

Production-ready genomics analysis pipelines for cattle research, developed at the University of Florida.

## Overview

UF-LGP provides automated workflows for:
- **Illumina SNP Array Processing** - Convert FinalReport files to PLINK format
- **Neogen Skim-Seek VCF Merging** - Merge per-sample sequencing data
- **Cross-Platform Integration** - Combine array and sequence genotypes
- **Breed Composition Analysis** - ADMIXTURE-based ancestry inference
- **Parentage Verification** - KING-based relatedness and parent identification

All pipelines are designed for SLURM-managed HPC environments.

## Quick Start

```bash
# Clone and setup
git clone https://github.com/gzayasPR/UF-LGP.git
cd UF-LGP
bash Setup.sh

# Run a pipeline (example: breed composition)
sbatch bin/pipelines/breed_composition.sh
```

## Repository Structure

```
UF-LGP/
├── bin/
│   ├── pipelines/          # Main pipeline entry points
│   └── scripts/            # Core processing scripts
├── metadata/               # Sample IDs and reference files
├── docs/                   # Detailed pipeline documentation
├── results/                # Output files (not tracked)
└── bash_out/               # SLURM logs (not tracked)
```

## Pipelines

### 1. Illumina Array Processing
Convert Illumina Bovine 250K FinalReport files to PLINK format with ARS-UCD1.2 coordinates.

```bash
sbatch bin/pipelines/illumina_2_plink.sh
```
📖 **[Full Documentation](docs/illumina_processing.md)**

### 2. Skim-Seek VCF Merging
Merge per-sample Skim-Seek VCFs with confidence-aware filtering.

```bash
sbatch bin/pipelines/skim_seek.sh
```
📖 **[Full Documentation](docs/skim_seek.md)**

### 3. Cross-Platform Integration
Merge Skim-Seek and Illumina 250K on common SNP sites.

```bash
sbatch bin/pipelines/SSK_250K.sh
```
📖 **[Full Documentation](docs/cross_platform.md)**

### 4. Breed Composition
Supervised ADMIXTURE analysis with PCA visualization.

```bash
sbatch bin/pipelines/breed_composition.sh
```
📖 **[Full Documentation](docs/breed_composition.md)**

### 5. Parentage Verification
KING-based parent-offspring identification.

```bash
sbatch bin/pipelines/check_parentage.sh
```
📖 **[Full Documentation](docs/parentage.md)**

## Requirements

- **Core Tools:** PLINK 1.9+, PLINK2, ADMIXTURE, KING, bcftools, htslib
- **R Packages:** SNPRelate, GENESIS, ggplot2
- **System:** SLURM workload manager, 8-32 GB RAM per job

See [Installation Guide](docs/installation.md) for detailed setup instructions.

## Example Workflow

```bash
# 1. Process Illumina genotypes
sbatch bin/pipelines/illumina_2_plink.sh

# 2. Analyze breed composition
# Edit metadata/IDs/admix.ID and metadata/IDs/MAB.ID first
sbatch bin/pipelines/breed_composition.sh

# 3. Check results
ls results/GBC_*/Figures/admixture.png
cat results/GBC_*/ADMIXTURE/Supervised_Unrelated/*.combined.ancestry.csv
```

## Documentation

- **[Installation & Setup](docs/installation.md)** - Dependencies and environment setup
- **[Pipeline Guides](docs/)** - Detailed documentation for each pipeline
- **[Troubleshooting](docs/troubleshooting.md)** - Common issues and solutions
- **[File Formats](docs/file_formats.md)** - Input/output specifications

## Citation

```
Zayas, G. (2026). UF-LGP: University of Florida Livestock Genomics Pipelines.
GitHub: https://github.com/gzayasPR/UF-LGP
```

## Contact

**Author:** Gabriel Zayas  
**Institution:** University of Wisconsin-Madison, Department of Animal & Dairy Sciences  
**Email:** zayas@wisc.edu | gzayas97@ufl.edu (closes May 2026) | gabboz12@gmail.com

For questions or issues, please open a GitHub issue.

## License

Provided for academic and research use. Contact University of Florida for commercial licensing.

---

**Developed in the Mateescu Lab | University of Florida**
