# Installation & Setup Guide

Complete guide for setting up UF-LGP on an HPC system.

## System Requirements

### Minimum Requirements

- **Operating System:** Linux (Ubuntu 18.04+, CentOS 7+, or similar)
- **Scheduler:** SLURM workload manager
- **Storage:** 100GB+ available space
- **Memory:** 32GB+ RAM for typical analyses
- **CPUs:** Multi-core system recommended (8+ cores)

### Tested Environments

- University of Florida HiPerGator 3.0 (SLURM 20.11)
- Ubuntu 24 LTS
- CentOS 7

---

## Software Dependencies

### Core Tools

| Tool | Minimum Version | Purpose |
|------|----------------|---------|
| PLINK | 1.90b3.39+ | Genotype QC and format conversion |
| PLINK2 | 2.0+ | Modern PLINK formats |
| ADMIXTURE | 1.3.0+ | Ancestry inference |
| KING | 2.3.0+ | Relatedness estimation |
| bcftools | 1.22+ | VCF manipulation |
| htslib | 1.22+ | VCF indexing (tabix, bgzip) |
| R | 3.6+ | Statistical analysis and visualization |
| samtools | 1.10+ | (Optional) Additional utilities |

### R Packages

#### Bioconductor Packages
- `SNPRelate` (1.26.0+) - Genomic data structure and PCA
- `GENESIS` (2.22.0+) - PC-AiR implementation
- `gdsfmt` (1.28.0+) - GDS format support

#### CRAN Packages
- `ggplot2` (3.3.0+) - Plotting
- `dplyr` (1.0.0+) - Data manipulation
- `tidyr` (1.1.0+) - Data reshaping

---

## Installation Steps

### 1. Clone Repository

```bash
# Navigate to desired location
cd /path/to/your/projects

# Clone repository
git clone https://github.com/gzayasPR/UF-LGP.git
cd UF-LGP

# Verify structure
ls -la
```

Expected contents:
```
UF-LGP/
├── bin/
├── docs/
├── metadata/
├── README.md
└── Setup.sh
```

### 2. Run Setup Script

```bash
# Execute setup
bash Setup.sh

# This creates:
# - bin/project_env.sh (environment variables)
# - results/, docs/, bash_out/ directories
```

Verify setup:
```bash
source bin/project_env.sh
echo $my_bin
echo $my_results
```

### 3. Install Software Dependencies

#### Option A: Using Environment Modules (Recommended for HPC)

```bash
# Check available modules
module avail

# Load required modules
module load plink/1.90b3.39
module load plink2/2.0
module load admixture/1.3.0
module load king/2.3.0
module load bcftools/1.22
module load htslib/1.22
module load R/4.1.0
module load samtools/1.10

# Save module list for future use
module save uf-lgp
```

To load saved modules in future sessions:
```bash
module restore uf-lgp
```

#### Option B: Manual Installation

**PLINK 1.9:**
```bash
cd /opt/software  # or your preferred location
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
unzip plink_linux_x86_64_20231211.zip
sudo mv plink /usr/local/bin/
```

**PLINK2:**
```bash
wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip
unzip plink2_linux_x86_64_20240105.zip
sudo mv plink2 /usr/local/bin/
```

**ADMIXTURE:**
```bash
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar -xzf admixture_linux-1.3.0.tar.gz
cd admixture_linux-1.3.0
sudo cp admixture /usr/local/bin/
```

**KING:**
```bash
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzf Linux-king.tar.gz
sudo mv king /usr/local/bin/
```

**bcftools & htslib:**
```bash
# Install dependencies
sudo apt-get install libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

# Download and compile
wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
tar -xjf bcftools-1.22.tar.bz2
cd bcftools-1.22
./configure --prefix=/usr/local
make
sudo make install
```

**R:**
```bash
# Ubuntu/Debian
sudo apt-get install r-base r-base-dev

# CentOS/RHEL
sudo yum install R
```

### 4. Install R Packages

```bash
# Start R
R

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SNPRelate",
    "GENESIS",
    "gdsfmt"
))

# Install CRAN packages
install.packages(c(
    "ggplot2",
    "dplyr",
    "tidyr"
))

# Verify installation
library(SNPRelate)
library(GENESIS)
library(ggplot2)

# Exit R
quit(save="no")
```

### 5. Verify Installation

```bash
# Test command-line tools
plink --version
plink2 --version
admixture --help
king -h
bcftools --version
tabix --version

# Test R packages
R -e "library(SNPRelate); library(GENESIS); library(ggplot2)"
```

All commands should run without errors.

---

## Configuration

### Environment Variables

Edit `bin/project_env.sh` if needed:

```bash
#!/bin/bash
export proj_dir="/path/to/UF-LGP"
export my_bin="${proj_dir}/bin"
export my_metadata="${proj_dir}/metadata"
export my_results="${proj_dir}/results"
export my_docs="${proj_dir}/docs"
export my_bash="${proj_dir}/bash_out"
```

Source this file in your scripts:
```bash
source /path/to/UF-LGP/bin/project_env.sh
```

### SLURM Configuration

Default SLURM settings in pipelines:
```bash
#SBATCH --account=mateescu        # Change to your account
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your@email.edu  # Change to your email
```

Edit these in each pipeline script (`bin/pipelines/*.sh`) as needed.

---

## Directory Structure Setup

After installation, your directory should look like:

```
UF-LGP/
├── bin/
│   ├── pipelines/              # Main entry points
│   │   ├── breed_composition.sh
│   │   ├── check_parentage.sh
│   │   ├── illumina_2_plink.sh
│   │   ├── skim_seek.sh
│   │   └── SSK_250K.sh
│   ├── scripts/                # Core processing scripts
│   │   ├── combine_Q_FAM.R
│   │   ├── gbc_admixture.sh
│   │   ├── gbc_setup_unrelated.sh
│   │   ├── illumina_prep.sh
│   │   ├── king.sh
│   │   ├── make_ill_plink.sh
│   │   ├── merge_skim_seek.sh
│   │   ├── merge_SSK_250.sh
│   │   ├── PCA.plots.R
│   │   ├── Q.PLOT.R
│   │   ├── Selection_of_Unrelated_script.R
│   │   └── shell.functions.sh
│   └── project_env.sh          # Environment variables
├── metadata/                   # Your sample IDs and reference files
│   ├── IDs/
│   └── Illumina_to_Plink/
├── docs/                       # Documentation
├── results/                    # Pipeline outputs (created as needed)
├── bash_out/                   # SLURM logs (created as needed)
├── README.md
└── Setup.sh
```

---

## Prepare Reference Files

### For Illumina Processing

Place these in `metadata/Illumina_to_Plink/`:

1. **SNP_Map.txt** - Maps SNP names to genome coordinates
2. **Update IDs file** - Maps old to new sample IDs
3. **Chr update file** - Chromosome assignments
4. **BP update file** - Base pair positions (ARS-UCD1.2)
5. **INDELS file** - Allele reference

### For Breed Composition

Create sample ID lists in `metadata/IDs/`:

```bash
# Reference (purebred) population
cat > metadata/IDs/MAB.ID << EOF
REF Angus001
REF Angus002
REF Brahman001
REF Brahman002
EOF

# Mixed population
cat > metadata/IDs/admix.ID << EOF
UF Animal001
UF Animal002
EOF
```

---

## Test Installation

### Quick Test: PLINK

```bash
# Create test data
plink --dummy 100 10000 --make-bed --out test

# Run basic command
plink --bfile test --freq --out test_freq

# Check output
head test_freq.frq

# Cleanup
rm test.*
```

### Quick Test: Breed Composition Pipeline

```bash
# 1. Prepare test data
plink --dummy 100 50000 --make-bed --out test_geno

# 2. Create test ID lists
head -30 test_geno.fam | awk '{print $1,$2}' > metadata/IDs/ref_test.ID
tail -70 test_geno.fam | awk '{print $1,$2}' > metadata/IDs/mixed_test.ID

# 3. Edit breed_composition.sh with test paths
vim bin/pipelines/breed_composition.sh
# Set:
#   IN_GENO="test_geno"
#   reflist="metadata/IDs/ref_test.ID"
#   mixedlist="metadata/IDs/mixed_test.ID"
#   NAME="TEST_RUN"

# 4. Submit test job
sbatch bin/pipelines/breed_composition.sh

# 5. Monitor
squeue -u $USER

# 6. Check results
ls results/TEST_RUN/
```

---

## Updating

To update UF-LGP:

```bash
cd /path/to/UF-LGP

# Save any local changes
git stash

# Pull updates
git pull origin main

# Reapply local changes if needed
git stash pop
```

---

## Uninstallation

To remove UF-LGP:

```bash
# Remove repository
rm -rf /path/to/UF-LGP

# Optionally remove software if manually installed
sudo rm /usr/local/bin/plink
sudo rm /usr/local/bin/plink2
sudo rm /usr/local/bin/admixture
sudo rm /usr/local/bin/king
sudo rm /usr/local/bin/bcftools
sudo rm /usr/local/bin/tabix

# Remove R packages if desired
R -e "remove.packages(c('SNPRelate', 'GENESIS', 'gdsfmt', 'ggplot2'))"
```

---

## Troubleshooting Installation

### Module Not Found

```bash
# List all available modules
module spider plink

# Load specific version
module load plink/1.90b3.39
```

### R Package Installation Fails

```bash
# Missing system dependencies
sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev

# Try installing one package at a time
R -e "BiocManager::install('SNPRelate')"
```

### Permission Denied

```bash
# Use sudo for system-wide install
sudo make install

# OR install to user directory
./configure --prefix=$HOME/software
make
make install
```

### SLURM Not Available

If running without SLURM:
- Remove `#SBATCH` lines from scripts
- Run scripts directly: `bash bin/pipelines/breed_composition.sh`
- Adjust resource usage manually

---

## Getting Help

### Check Versions

```bash
# Software versions
plink --version
bcftools --version
R --version

# R packages
R -e "packageVersion('SNPRelate')"
```

### Documentation

- **UF-LGP docs:** `docs/` folder in repository
- **PLINK:** https://www.cog-genomics.org/plink/
- **ADMIXTURE:** https://dalexander.github.io/admixture/
- **KING:** https://www.kingrelatedness.com/
- **bcftools:** https://samtools.github.io/bcftools/

### Support

- **Issues:** Open issue on GitHub repository
- **Email:** gzayas97@ufl.edu
- **HPC Support:** Contact your institutional HPC support team

---

## Next Steps

After installation:

1. **Read pipeline documentation** in `docs/` folder
2. **Prepare your metadata files** (sample IDs, reference files)
3. **Test with small dataset** before running on full data
4. **Review example workflows** in README and individual docs
5. **Run your first analysis!**
