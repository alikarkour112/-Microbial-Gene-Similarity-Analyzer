
# 🧬 Microbial Gene Similarity Analyzer

An interactive, high-performance web application designed for microbiologists, geneticists, and bioinformaticians. This tool enables users to compare gene sequences across multiple microorganism strains using **Global Sequence Alignment**, instantly computing similarity matrices and rendering publication-ready heatmaps.

---

## 🌟 Key Features

* **Self-Healing Bootstrap:** Upon its very first launch, the script dynamically verifies the active environment and automatically installs any missing dependencies (`biopython`, `seaborn`, `pandas`, `matplotlib`) using pip, removing manual setup bottlenecks.
* **Dual Input Channels:**
* **Reference Strain:** Supports official NCBI GenBank files (`.gb`, `.gbk`) to seamlessly parse the baseline sequence.
* **Query Strains:** Accepts multi-file batch uploads in FASTA format (`.fasta`, `.fa`) for parallel strain processing.


* **Built-in Interactive Demo:** Launches instantly out of the box with embedded synthetic data (WildType vs. Mutants and Divergent strains) so users can test and experience the interface without needing external files.
* **Dynamic Sidebar Calibration:** Provides full slider/input control over the Needleman-Wunsch alignment variables (Match Score, Mismatch Penalty, Gap Open, and Gap Extend).
* **Advanced Visual Dashboard:**
* An automated summary dataframe tracking strain sequences and their length in base pairs (bp).
* A beautifully formatted similarity matrix equipped with live `viridis` color gradient mapping.
* A high-resolution Seaborn/Matplotlib heatmap plot using the `YlGnBu` spectrum, adjusted with tailored axis rotation for direct integration into academic papers.


* **Embedded Biological Interpretation:** Displays an interpretive evaluation framework at the footer of the page to guide users through their structural findings.

---

## 🛠️ Tech Stack & Architecture

* **Programming Language:** Python 3.11+
* **User Interface:** Streamlit (Configured in Wide Layout mode)
* **Bioinformatics Engine:** Biopython (`Bio.SeqIO` & `Bio.Align.PairwiseAligner`)
* **Data Processing:** NumPy & Pandas
* **Data Visualization:** Matplotlib & Seaborn

---

## 🚀 Installation & Local Setup

### 1. Prerequisites

Ensure you have Python 3.8 or a newer version installed on your operating system.

### 2. Download the Code

Save the application script into a file named **`app.py`** inside your local project directory.

### 3. Run the Streamlit Server

Open your terminal or command prompt, navigate into the directory containing `app.py`, and run the following command:

```bash
streamlit run app.py

```

> 💡 **Automated Dependency Tip:** You do not need to manually install dependencies. The script's dynamic startup framework checks for packages and automatically handles downloading them in the background during the initial launch.

---

## 📖 Step-by-Step Usage Guide

```
[ Application Interface Map ]
 ├── 📊 Sidebar Panel: Calibrate rewards & penalties (Match, Mismatches, Gap behaviors)
 ├── 📂 File Ingestion Zones:
 │    ├── 1. Drop a single GenBank File (Reference Sequence)
 │    └── 2. Drop multiple FASTA Files (Query Target Sequences)
 └── 📈 Analytics Workspace: Triggers automatically (Summary -> Color Matrix -> Heatmap Canvas)

```

1. **Fine-Tune Weights (Optional):** Use the sidebar to update open or extension gap penalties according to the specific evolutionary distance or structural deletion expectations of your target organism.
2. **Review Default Mockup:** When opening the application, explore the synthetic data preview matrices to understand how data outputs change based on alignment conditions.
3. **Analyze Custom Genomes:** Drag your raw `.gb`/`.gbk` reference file into the first portal and your query `.fasta` sequences into the second. A live progress bar will trace the multi-sequence computing matrix.

---

## 🧬 Core Algorithmic & Biological Logic

### Alignment Engine

The analytical kernel is governed by a global pairwise iteration via Biopython's alignment suite. It aligns sequences end-to-end to optimize regions of structural agreement while calculating penalties for:

* **Point Mutations (Substitutions):** Weighted via the `Mismatch Penalty` setting.
* **Insertions & Deletions (Indels):** Controlled tightly via the `Gap Open` and `Gap Extend` settings.

### Interpretation Metric System

* **Diagonal Control Values (100.0%):** Represents an isolated alignment of each sequence against itself, establishing a static control benchmark.
* **High Similarity (> 95%):** Points to highly conserved sequences within the exact same species. Variances generally signal simple Single Nucleotide Polymorphisms (SNPs).
* **Moderate Similarity (70% - 95%):** Points to evolutionary divergence. This degree of variance is standard when evaluating differing variants of concern (VOCs), separate serotypes, or lineages adapting to unique selective environments.
* **Low Similarity (< 70%):** Reveals massive genetic drift, large-scale gene dropouts, or indicates that the compared strains belong to completely isolated genus classes altogether.
