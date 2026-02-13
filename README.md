# ScarMapper

A Python package that uses an iterative break-associated alignment strategy to classify individual double-strand DNA break repair products based on deletion size, microhomology usage, and insertions.

## Getting Started

### Prerequisites

ScarMapper runs on Linux (tested on RHEL 7.x, Scientific Linux 7.x, CentOS 7.x) and may work on macOS, though it has not been tested. It will not run on Windows because ScarMapper uses Pysam, which requires a POSIX environment.

**Minimum System Requirements:**
- 4 CPUs or threads
- 20 GB RAM
- ~5x the size of the compressed FASTQ file for disk space

**Software Requirements:**
- Python >= 3.8
- [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) (for paired-end read merging)

All Python dependencies are installed automatically via pip.

### Installation

Clone or download this repository, then install with pip:

```bash
git clone https://github.com/AllanJWLui/ScarMapper.git
cd ScarMapper
pip install .
```

For development (editable install):

```bash
pip install -e ".[dev]"
```

Verify installation:

```bash
scarmapper --help
```

### Installation with Conda (recommended for HPC)

This installs ScarMapper and the PEAR binary together:

```bash
git clone https://github.com/AllanJWLui/ScarMapper.git
cd ScarMapper
conda env create -f environment.yml
conda activate scarmapper
```

Alternatively, if PEAR is already available (e.g. via `module load`), you can install ScarMapper with pip alone.

### Running with Docker / Singularity

For HPC clusters where Docker is unavailable, use Singularity to run containerized ScarMapper.

#### Pull the image

```bash
# Singularity (on HPC)
singularity pull docker://ghcr.io/allanjwlui/scarmapper:3.0.1
# Creates: scarmapper_3.0.1.sif

# Docker (local machine)
docker pull ghcr.io/allanjwlui/scarmapper:3.0.1
```

#### Quick start with Singularity

```bash
# Show help
singularity exec scarmapper_3.0.1.sif scarmapper --help

# Run with bind mounts
singularity exec \
  --bind /path/to/your/data:/data \
  --bind /path/to/your/output:/output \
  scarmapper_3.0.1.sif \
  scarmapper --options_file /data/ScarMapper_IndelProcessing.cfg
```

#### Full workflow example

```bash
# Setup
cd /scratch/your_username/scarmapper_run
singularity pull docker://ghcr.io/allanjwlui/scarmapper:3.0.1

# Organize files
mkdir -p data output
cp your_R1.fastq.gz your_R2.fastq.gz data/
cp reference.fasta reference.fasta.fai data/
cp *.tsv data/  # manifest, targets, index files
cp ScarMapper_IndelProcessing.cfg data/

# IMPORTANT: Edit the .cfg file to use container paths
# All paths in the config must use /data/ or /output/ prefixes:
#   Reference_Genome    /data/reference.fasta
#   FASTQ1              /data/your_R1.fastq.gz
#   Target_File         /data/targets.tsv
#   Working_Folder      /output/

# Run ScarMapper
singularity exec \
  --bind $(pwd)/data:/data \
  --bind $(pwd)/output:/output \
  scarmapper_3.0.1.sif \
  scarmapper --options_file /data/ScarMapper_IndelProcessing.cfg

# Results will be in ./output/
```

#### HPC job script example

```bash
#!/bin/bash
#SBATCH --job-name=scarmapper
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4:00:00

module load singularity  # if needed

cd /scratch/$USER/scarmapper_project

singularity exec \
  --bind ./data:/data \
  --bind ./output:/output \
  scarmapper_3.0.1.sif \
  scarmapper --options_file /data/config.cfg
```

#### Running with Docker

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/output \
  ghcr.io/allanjwlui/scarmapper:3.0.1 \
  --options_file /data/config.cfg
```

#### Build locally

```bash
docker build -t scarmapper .
```

#### Important notes

- **Bind mounts**: Use `--bind local_path:container_path` to make your files accessible inside the container
- **Config file paths**: All file paths in `.cfg` options files must use container-side paths (`/data/`, `/output/`), not your local filesystem paths
- **Read-only mounts**: Add `:ro` suffix for read-only access (e.g., `--bind /genomes:/genomes:ro`)
- **Multiple bind mounts**: You can bind multiple directories as needed

### Input Files

ScarMapper requires several input files. Templates are provided in the `docs/` directory.

| File | Description |
|------|-------------|
| **FASTQ files** | Paired-end R1 and R2 reads (gzipped or uncompressed) |
| **Reference genome** | FASTA file with an index (`.fai`) in the same directory |
| **Master index file** | Tab-delimited file mapping index names to forward/reverse barcode sequences |
| **Sample manifest** | Tab-delimited file mapping index names to sample names, replicates, loci, and optionally per-sample FASTQ paths and HR donor sequences |
| **Target file** | Tab-delimited file defining each genomic target: name, chromosome, start, stop, sgRNA sequence, and strand orientation |

See `docs/ScarMapper_Manifest.tsv` and `docs/ScarMapper_Targets.tsv` for example formats.

## Usage

ScarMapper has three processing modes, each driven by an options file:

```bash
scarmapper --options_file <config_file>.cfg
```

### Indel Processing (multiplexed FASTQ)

For standard runs with multiplexed FASTQ files containing multiple samples:

```bash
scarmapper --options_file ScarMapper_IndelProcessing.cfg
```

This mode demultiplexes reads by index, merges paired-end reads using PEAR, then identifies and classifies repair scars at each target locus.

### Batch Processing (pre-demultiplexed samples)

For runs where each sample already has its own FASTQ file(s):

```bash
scarmapper --options_file ScarMapper_Batch.cfg
```

The sample manifest must include a `FASTQ1_Path` column (and optionally `FASTQ2_Path`) pointing to each sample's files. Per-sample HR donor sequences can also be specified in the manifest.

### Combine Replicates

To merge frequency files from replicate experiments and produce a combined plot:

```bash
scarmapper --options_file ScarMapper_Combine.cfg
```

### Options Files

Options files are tab-delimited configuration files (`.cfg`). Each option must be separated from its value by exactly **one tab character**. Template options files for all three modes are provided in the `docs/` directory:

- `docs/ScarMapper_IndelProcessing.cfg`
- `docs/ScarMapper_Batch.cfg`
- `docs/ScarMapper_Combine.cfg`

### Platform Support

The `--Platform` option controls how ScarMapper identifies sample indices on each read:

| Platform | Index Location | Mismatch Tolerance |
|----------|---------------|-------------------|
| **Illumina** | Parsed from FASTQ header after last `:`, split on `+` | 1 |
| **TruSeq** | First and last 6 nucleotides of read sequence | 0 (exact match) |
| **Ramsden** | Start of R1 and start of R2 (or end of R1 if PEAR merged) | 3 |

## Output Files

ScarMapper produces several output files in the working folder:

| File | Description |
|------|-------------|
| `*_ScarMapper_Summary.txt` | Per-library summary with scar counts and repair pathway fractions |
| `*_ScarMapper_Frequency.txt` | Per-library scar pattern frequencies with junction details |
| `*_ScarMapper_Raw_Data.txt` | Per-read scar data (if `--OutputRawData True`) |
| `*_SNV_Frequency.txt` | SNV positional frequency data (if SNVs detected) |
| `*.<FigureType>` | Scar pattern plots for each library |
| `*_Batch_Summary.txt` | Combined summary for batch mode |

## Scar Classification

ScarMapper classifies each repair scar into one of six categories:

| Category | Criteria |
|----------|----------|
| **TMEJ** | Deletion >= 4 nt and microhomology >= 2 nt |
| **NHEJ** | Deletion < 4 nt and insertion < 5 nt |
| **Non-MH Deletion** | Deletion >= 4 nt, microhomology < 2 nt, and insertion < 5 nt |
| **Insertion** | Insertion >= 5 nt, with or without deletions |
| **SNV** | Deletion > 1 nt, deletion size equals insertion size, within kmer bounds |
| **HR** | Homologous recombination donor sequence detected (requires `--HR_Donor`) |

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Authors

* **Allan Lui**
* **Dennis Simpson** — *Initial work*

## License

This project is licensed under the MIT License. See [LICENSE.md](LICENSE.md) for details.
