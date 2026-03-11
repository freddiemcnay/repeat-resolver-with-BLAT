# Repeat Structure Classifier with BLAT Integration

A bioinformatics tool for analysing long sequencing reads containing complex repeat expansions. Identifies and quantifies tandem repeat structures by aligning sliding windows against known repeat unit sequences using minimap2. Undetermined (U) regions — internal stretches of sequence that cannot be assigned to any known unit — are automatically extracted and optionally mapped to a reference genome via a local BLAT server.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Parameters](#parameters)
- [Setting Up BLAT on farm22](#setting-up-blat-on-farm22)
- [Available Genomes](#available-genomes)
- [Adding a New Genome](#adding-a-new-genome)
- [Running on farm22 (LSF)](#running-on-farm22-lsf)
- [Understanding Results](#understanding-results)
- [Troubleshooting](#troubleshooting)

---

## Features

- **Strand-aware analysis** — correctly interprets reverse complement alignments
- **Tandem repeat detection** — automatically identifies and counts tandem copies using uppercase/lowercase logic
- **Quality metrics** — per-window alignment scores for troubleshooting
- **Population summary** — frequency table showing dominant repeat configurations across a dataset
- **Flexible parameters** — adjustable window size, identity threshold, and minimap2 presets
- **U-region extraction** — reads with internal undetermined stretches are automatically extracted into their own FASTA files, with U windows highlighted in uppercase
- **Local BLAT mapping** — U regions are queried against a reference genome via a local gfServer, with top genomic hits reported per region and annotated directly in the structures output

---

## Requirements

- Python 3.6+
- minimap2 (must be in PATH)
- BLAT tools: `gfServer`, `gfClient`, `faToTwoBit`, `twoBitInfo` (for BLAT functionality)

All BLAT binaries are installed at:
```
/lustre/scratch125/casm/staging/team267_murchison/fm15/programmes/blat/
```

---

## Quick Start

```bash
# Basic usage — repeat structure classification only
./repeat_resolver_with_BLAT.py \
    -i reads.fasta \
    -u units.fasta \
    -o output_prefix

# With BLAT against canFam3
./repeat_resolver_with_BLAT.py \
    -i reads.fasta \
    -u units.fasta \
    -o output_prefix \
    --blat \
    --blat-host $(cat ~/blat_server_canfam3_host.txt) \
    --blat-port $(cat ~/blat_server_canfam3_port.txt) \
    --blat-genome-dir /lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canFam3

# With BLAT against custom canFam4
./repeat_resolver_with_BLAT.py \
    -i reads.fasta \
    -u units.fasta \
    -o output_prefix \
    --blat \
    --blat-host $(cat ~/blat_server_canfam4_host.txt) \
    --blat-port $(cat ~/blat_server_canfam4_port.txt) \
    --blat-genome-dir /lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canfam4
```

---

## How It Works

### Repeat Classification

The classifier uses a sliding window approach combined with minimap2 alignment to identify repeat structures:

1. **Divide read into windows** — each input read is split into non-overlapping windows (default 1000 bp)
2. **Align each window to all units** — for every window, minimap2 is run separately against each unit sequence (A, B, C, D, etc.)
3. **Select best match** — percent identity is calculated for all alignments; the unit with the highest score is assigned
   - If best identity ≥ threshold (default 40%) → assign that unit
   - If best identity < threshold → assign `U` (undetermined)
4. **Determine uppercase/lowercase** — check where the window aligns on the reference unit
   - Aligns to first 30% of unit → Uppercase (A, B, C) = start of tandem repeat
   - Aligns to last 70% of unit → Lowercase (a, b, c) = continuation of tandem repeat
5. **Track strand orientation** — overall read orientation determined by majority vote across all assigned windows
6. **Build and simplify structure** — windows are concatenated into a structure string, then collapsed into a simplified notation with tandem copy counts

### U-Region Analysis

After classification, reads are checked for **internal** U windows — stretches of undetermined sequence that are not simply a short trailing window at the end of the read (which are ignored as they often represent incomplete final windows with too few bases to align confidently).

For reads with internal U regions:
- A **context FASTA** is written with the full read sequence, U windows in uppercase and all other bases in lowercase
- A **U-regions FASTA** is written with just the undetermined sequence extracted, one entry per contiguous U stretch, with coordinates encoded in the FASTA header

### BLAT Mapping

If `--blat` is enabled, all U-region sequences are submitted as a single batch to a local `gfServer` running on a farm22 compute node. The server holds the reference genome in RAM and responds to queries from `gfClient`. Everything runs within Sanger's internal network — nothing is sent to the internet. The top genomic hits are recorded per region and the best hit per U region is annotated directly in the structures output file.

---

## Input Files

### Reads file (`-i`)

FASTA file containing the long reads to classify. Supported read types: PacBio HiFi, PacBio CLR, Oxford Nanopore.

```
>read001
ATCGATCGATCG...
>read002
GCTAGCTAGCTA...
```

### Units file (`-u`)

FASTA file containing the known repeat unit sequences to classify against. Use simple single-letter names.

```
>A
ATCGATCGATCGATCG...
>B
GCTAGCTAGCTAGCTA...
>C
TTAATTAATTAATTAA...
```

---

## Output Files

### `*_structures.txt`

Tab-delimited file with one line per read. The `u_blat_top_hits` column is only present when `--blat` is used.

| Column | Description |
|--------|-------------|
| `read_name` | Unique read identifier |
| `strand` | Overall orientation (+, -, mixed, unknown) |
| `full_structure` | Window-by-window assignments (e.g. `A-a-a-B-C-U-U-D`) |
| `simplified_structure` | Collapsed with tandem counts (e.g. `A-B-C-U-D`) |
| `u_blat_top_hits` | Top BLAT hit per U region, semicolon-separated (e.g. `chr13:16682663-16689663`). `.` if no hit or BLAT not run. |

### `*_detailed_scores.txt`

One line per window with alignment metrics for each unit. Useful for understanding assignment decisions, identifying ambiguous regions, and quality control.

### `*_structure_summary.txt`

Population-level frequency table showing the most common simplified structures, strand distribution, and count of reads with internal U regions.

### `*_u_context.fasta`

One entry per read that has internal U regions. The full read sequence is written with **U windows in uppercase** and all other bases in lowercase.

```
>read001 strand=+ structure=A-U-B
acgtacgtacgt...ACGTACGTacgtacgt...
               ^^^^^^^^
               U region (uppercase)
```

### `*_u_regions.fasta`

Just the undetermined sequences, one entry per contiguous U stretch. The FASTA header encodes the read name, base-pair coordinates within the read, and window numbers:

```
>read001_Uregion1_bp3000-5000_windows4-5
ACGTACGTACGT...
```

### `*_u_blat_hits.tsv`

Detailed BLAT results for all U regions. Up to `--blat-top-hits` rows per region, ranked by BLAT score.

| Column | Description |
|--------|-------------|
| `region_id` | U region identifier (from `*_u_regions.fasta` header) |
| `rank` | Hit rank ordered by score (1 = best) |
| `chrom` | Chromosome |
| `strand` | Alignment strand (+/-) |
| `t_start` | Genomic start position (0-based) |
| `t_end` | Genomic end position |
| `matches` | Number of matching bases |
| `mismatches` | Number of mismatching bases |
| `percent_identity` | matches / (matches + mismatches) × 100 |
| `q_size` | Query sequence length (bp) |
| `span` | Genomic span of the alignment: `t_end - t_start` |
| `score` | BLAT score: matches - mismatches |
| `block_count` | Number of alignment blocks |

### `*_u_regions.psl`

Raw PSL file output from gfClient. Contains all BLAT hits before filtering. Useful for custom downstream analysis.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-i, --input` | Required | Input FASTA with reads to classify |
| `-u, --units` | Required | FASTA with known repeat unit sequences |
| `-o, --output` | Required | Output prefix for all results files |
| `-w, --window-size` | 1000 | Window size in base pairs |
| `-t, --identity-threshold` | 40.0 | Min % identity for unit assignment (below = U) |
| `-c, --uppercase-threshold` | 0.3 | Fraction of unit length for uppercase assignment |
| `-p, --preset` | map-hifi | Minimap2 preset: map-hifi, map-ont, map-pb, asm20, sr |
| `--blat` | False | Enable local BLAT mapping of U regions |
| `--blat-host` | None | Hostname where gfServer is running (required with `--blat`) |
| `--blat-port` | 7777 | gfServer port number |
| `--blat-genome-dir` | None | Directory containing the .2bit file (required with `--blat`) |
| `--blat-min-identity` | 70.0 | Min % identity for reported BLAT hits |
| `--blat-min-score` | 20 | Min score for reported BLAT hits |
| `--blat-top-hits` | 3 | Number of top BLAT hits to report per U region |

---

## Setting Up BLAT on farm22

BLAT requires a `gfServer` to be running with the reference genome loaded into RAM before running the script with `--blat`. The server takes around 3 minutes to start and uses approximately 4-8GB of RAM per genome.

The server runs as an LSF job on a farm22 compute node. When it starts, it writes the node hostname and port to files in your home directory, which your analysis jobs read at submission time to know where to connect.

---

## Available Genomes

Two genomes are currently set up, each with its own LSF startup script and port.

### canFam3

| | |
|-|---------|
| **LSF script** | `start_blat_server_canfam3.lsf` |
| **Port** | 7777 |
| **Host file** | `~/blat_server_canfam3_host.txt` |
| **Port file** | `~/blat_server_canfam3_port.txt` |
| **2bit file** | `/lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canFam3/canFam3.2bit` |
| **Genome dir** | `/lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canFam3` |

```bash
bsub < start_blat_server_canfam3.lsf
```

### Custom canFam4 (plusY plusMT)

| | |
|-|---------|
| **LSF script** | `start_blat_server_canfam4_custom.lsf` |
| **Port** | 7778 |
| **Host file** | `~/blat_server_canfam4_host.txt` |
| **Port file** | `~/blat_server_canfam4_port.txt` |
| **2bit file** | `/lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canfam4/canFam4_plusYplusMT.2bit` |
| **Genome dir** | `/lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canfam4` |

```bash
bsub < start_blat_server_canfam4_custom.lsf
```

Both servers can run simultaneously on different ports.

---

## Adding a New Genome

### Step 1: Convert FASTA to .2bit

gfServer requires the genome in `.2bit` format. If you only have a FASTA:

```bash
faToTwoBit /path/to/genome.fa /path/to/genome.2bit

# Verify
twoBitInfo /path/to/genome.2bit stdout | head -10
```

If a `.2bit` already exists, skip this step.

### Step 2: Choose a port

Pick an unused port. Current assignments are:

| Genome | Port |
|--------|------|
| canFam3 | 7777 |
| canFam4 custom | 7778 |

Use the next available number (e.g. 7779).

### Step 3: Create a new LSF startup script

Copy an existing script and edit the four lines that change:

```bash
cp start_blat_server_canfam3.lsf start_blat_server_NEWGENOME.lsf
```

Edit these lines in the new file:

```bash
#BSUB -J blat_gfserver_NEWGENOME
#BSUB -o .../logs/blat_server_NEWGENOME_%J.log
#BSUB -e .../logs/blat_server_NEWGENOME_%J.err

GENOME_2BIT=/path/to/new/genome.2bit
PORT=7779

echo $(hostname) > ~/blat_server_NEWGENOME_host.txt
echo ${PORT}     > ~/blat_server_NEWGENOME_port.txt

${BLAT_BIN}/gfServer start $(hostname) ${PORT} ${GENOME_2BIT} \
    -stepSize=5 \
    -log=.../logs/blat_server_NEWGENOME_gf.log
```

### Step 4: Start and verify

```bash
bsub < start_blat_server_NEWGENOME.lsf

# Wait ~3 minutes, then check
gfServer status $(cat ~/blat_server_NEWGENOME_host.txt) 7779
```

### Step 5: Run the script against the new genome

```bash
./repeat_resolver_with_BLAT.py \
    -i reads.fasta \
    -u units.fasta \
    -o output_prefix \
    --blat \
    --blat-host $(cat ~/blat_server_NEWGENOME_host.txt) \
    --blat-port $(cat ~/blat_server_NEWGENOME_port.txt) \
    --blat-genome-dir /path/to/new/genome/directory
```

---

## Running on farm22 (LSF)

### Interactive test (small dataset)

```bash
./repeat_resolver_with_BLAT.py \
    -i reads.fasta \
    -u units.fasta \
    -o test_output \
    --blat \
    --blat-host $(cat ~/blat_server_canfam3_host.txt) \
    --blat-port $(cat ~/blat_server_canfam3_port.txt) \
    --blat-genome-dir /lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canFam3
```

### Submitting as an LSF job

```bash
bsub \
    -G team267-grp \
    -q normal \
    -o repeat_resolver.%J.log \
    -R "select[mem>4000] rusage[mem=4000]" \
    -R "span[hosts=1]" \
    -M 4000 \
    "./repeat_resolver_with_BLAT.py \
        -i reads.fasta \
        -u units.fasta \
        -o results \
        --blat \
        --blat-host $(cat ~/blat_server_canfam3_host.txt) \
        --blat-port $(cat ~/blat_server_canfam3_port.txt) \
        --blat-genome-dir /lustre/scratch125/casm/staging/team267_murchison/fm15/reference/canFam3"
```

**Important**: The `$(cat ~/blat_server_canfam3_host.txt)` substitution happens at submission time on the head node, so the correct hostname is baked into the job. Make sure the gfServer is running and those files exist before submitting.

### Checking the server before submitting

```bash
gfServer status $(cat ~/blat_server_canfam3_host.txt) $(cat ~/blat_server_canfam3_port.txt) \
    && echo "Server ready" \
    || echo "Server not responding — resubmit the LSF script"
```

---

## Understanding Results

### Structure Notation

| Symbol | Meaning |
|--------|---------|
| `A`, `B`, `C` | Start of a unit or single-copy unit (window aligns to first 30% of reference) |
| `a`, `b`, `c` | Continuation of a tandem copy (window aligns to last 70% of reference) |
| `U` | Undetermined — window identity below threshold |
| `×n` | Tandem copies in simplified structures, e.g. `C×3` = 3 tandem copies of C |

### Examples

| Full Structure | Simplified | Interpretation |
|----------------|------------|----------------|
| `C-c-c-D-d` | `C-D` | One copy of C, then one copy of D |
| `C-c-c-C-c-c-C-c-c` | `C×3` | Three tandem copies of C |
| `A-a-U-U-U-B-b` | `A-U-B` | Unit A, undetermined region, then unit B |
| `C-c-D-d-C-c-D-d` | `C-D-C-D` | Alternating C and D |

### U-Region Notation

**Trailing U** (e.g. `A-a-a-U`) — ignored. Almost always a short final window with too few bases to align. Not extracted, not BLATted.

**Internal U** (e.g. `A-U-U-U-B`) — extracted and analysed. Represents a genuine stretch of sequence that does not match any known unit. The BLAT hit in `u_blat_top_hits` tells you where in the reference genome this sequence originates.

### BLAT SPAN vs Query Size

- **`q_size`** — the length of the U-region sequence that was BLATted (bp)
- **`span`** — the genomic distance covered on the reference: `t_end - t_start`

These will differ when the alignment contains gaps relative to the reference, for example if the U region spans a repetitive element that is slightly expanded or contracted relative to the reference genome.

---

## Troubleshooting

### Too many U assignments

- Lower identity threshold: `-t 35` or `-t 30`
- Verify all expected units are in the units file
- Review `*_detailed_scores.txt` to see actual identity values per window

### BLAT returns no hits

- Check the gfServer is still running: `gfServer status $(cat ~/blat_server_canfam3_host.txt) $(cat ~/blat_server_canfam3_port.txt)`
- Try lowering `--blat-min-identity` (default 70) or `--blat-min-score` (default 20)
- Check `*_u_regions.psl` — if empty, gfClient ran but found nothing; if absent, gfClient failed
- Check the LSF log for gfClient error messages

### gfServer job gets killed

Resubmit the relevant LSF script and wait ~3 minutes for the genome to reload before rerunning your analysis.

### minimap2 errors

- Verify minimap2 is available: `which minimap2`
- Try a different preset: `-p map-ont` or `-p map-pb`
- Check FASTA formatting — no blank lines, proper `>header` lines
