# GenomeX — Advanced DNA Sequence Analysis System

![GenomeX Overview](https://img.shields.io/badge/Status-Development-blue.svg) ![Language](https://img.shields.io/badge/Core-C11-green) ![Backend](https://img.shields.io/badge/API-Flask-black) ![Frontend](https://img.shields.io/badge/UI-Vanilla_JS-orange)

GenomeX is a high-performance, web-based genomic analysis workstation. It bridges highly optimized bioinformatics algorithms written natively in C with a Python Flask REST API and an interactive, beautifully designed "Cyber Biotech Command Center" dashboard.

This system is built from the ground-up utilizing highly efficient data structures to analyze, visualize, and compare large-scale DNA data ranging from clinical mutated genes (BRCA1, HBB) to mitochondrial tracking for human evolution.

## 🧬 System Architecture

The project consists of three main operational tiers:
1. **C-Core Database (`c-core/`)**: Highly optimized string matching arrays and graphs implemented natively in C11. Outputs a shared library (`libdna.so` / `libdna.dll`).
2. **Python Application Layer (`api/`)**: Built on Flask, it dynamically connects to the C binary via Python's `ctypes` (`c_bridge.py`) for O(1) latency integration. 
3. **Frontend Dashboard (`frontend/`)**: A rich, responsive UI utilizing pure Vanilla JS and D3.js to execute dynamic `fetch()` requests and map out sequence discrepancies.

### Core Data Structures
The heavy lifting is spread across multiple advanced data structures assigned across the team:
* **Needleman-Wunsch DP Matrix (`needleman_wunsch.c`)**: Global DNA alignment and mismatch gap detection.
* **Trie with DFS Fuzzy Search (`trie_errors.c`)**: Clinical SNP tracking allowing k=0, 1, or 2 mismatches.
* **Suffix Tree / Ukkonen's (`suffix_tree.c`)**: Exact sequence identification, used for fast Pathogen Identification.
* **Segment Tree (`segment_tree.c`)**: High-speed position and coordinate compressed range queries across mutation loci.
* **Augmented Interval Tree (`interval_tree.c`)**: Mapping genomic intervals to overlap disease loci efficiently. 
* **Disjoint Set Union / Union-Find (`union_find.c`)**: Species evolutionary classification & clustered phylogenetic tree generation via Path Compression.

## 🚀 Installation & Setup

### 1. Fetching the Datasets
GenomeX relies heavily on real biological NCBI FASTA datasets and custom ClinVar mappings. First, download all the necessary subsets:
```bash
python fetch_data.py
```
*(This command queries NCBI E-Utilities to establish exactly 12 required genomic datasets into `data/human`, `data/mitochondrial`, and `data/pathogens`.)*

### 2. Building the C-Core Shared Library
The `Makefile` inside the `c-core/` directory binds all implemented C graphs into a generic library that the Python middleware can load.
```bash
cd c-core
make
```

### 3. Starting the Flask Server
Make sure you have Flask installed: `pip install flask flask-cors`
```bash
cd api
python app.py
```
The server will boot on `http://localhost:5000` and automatically load the ClinVar SNPs and Gene coordinates on initialization.

### 4. Running the Frontend
There is absolutely no internal dependency bundler required for the UI. Simply open the file:
```text
frontend/webpage.html
```
in your preferred web browser or via a VS Code Live Server.

## ⚙️ REST API Endpoints Overview

The following core boundaries are actively hooked into the Frontend UI:

* **`POST /api/pathogen/identify`**
  Uses the *Suffix Tree* database to classify unknown viral sequences down to exact pattern substrings.
* **`POST /api/mutations/scan`** 
  Utilizes the *Trie with Errors* and *Segment Tree* to scan an unknown genome string along a sliding window for exact or fuzzy mutations. 
* **`POST /api/evolution/cluster`** 
  Executes the *Union-Find* logic to identify common similarities between inputted matrices, powering D3 phylogenetic clustering.
* **`POST /api/snp/fuzzy`** 
  Maps disease risk evaluations natively from the *Trie DFS* to predict inherited dangers (like Breast Cancer BRCA1 strains).
* **`POST /api/align`** 
  Scores two comparative genomes over the *Needleman-Wunsch* DP matrix, yielding Matches, Mismatches, % Identity, and raw alignments.

## 🛠 Team Responsibilities
- **Integration Lead**: Non-DS Architecture, `trie_errors.c`, Makefile execution, `fetch_data.py`, Frontend binding hooks.
- **Member 1**: `suffix_tree.c`
- **Member 2**: `segment_tree.c`
- **Member 3**: `interval_tree.c`
- **Member 4**: `union_find.c`

*(Empty functional shell stubs exist for parallel teammate integrations so that the `Makefile` and `c_bridge.py` endpoints safely route payloads without crashing the server!)*
