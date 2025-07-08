# VaxLab-report

**RNA sequence evaluation and HTML report generation toolkit**

VaxLab-report evaluates mRNA or RNA sequences using diverse biophysical and sequence-based metrics related to translation efficiency, structural stability, and potential immunogenicity. It produces interactive HTML reports summarizing these findings.

---

## ✨ Features

- Supports evaluation of RNA sequences in FASTA format
- Calculates multiple metrics, including:
  - Codon Adaptation Index (CAI)
  - Minimum Free Energy (MFE) and secondary structure (ViennaRNA)
  - GC content (overall and per position)
  - Uridine content
  - Secondary structure features: stems, loops, start codon accessibility
  - Tandem repeats
  - Bicodon usage bias
  - (Optional) DegScore: decay propensity scoring
  - (Optional) IDT gBlock complexity scoring with API integration
- Generates rich HTML reports with:
  - Metric summaries
  - Interactive positional plots
  - RNA secondary structure visualization (Forna)
  - IDT complexity analysis with emoji status indicators
  - Downloadable source files
- Supports preset customization via JSON config
- Outputs raw scores in JSON/TSV format for further analysis

---

## ⚙️ Installation

```bash
# 1. Clone the repository
git clone https://github.com/ChangLabSNU/VaxLab-report.git
cd VaxLab-report

# 2. Set up the environment (Python 3.9)
conda create -y -n vaxpress_report python=3.9.22
conda activate vaxpress_report

# 3. Install required packages
conda install -y -c bioconda viennarna     # ViennaRNA for structure prediction
pip install -e .                           # Installs Python dependencies
```

---

## 🚀 Usage

### Step 1: Sequence evaluation

```bash
python vaxlab_report/evaluate_only.py \
  -i path/to/input.fasta \
  -o output_directory \
  --preset path/to/parameters.json \
  --species "Homo sapiens" \
  --token "YOUR_IDT_API_TOKEN"  # Optional: for IDT complexity scoring
```

**Generates:**
- `evaluation_result.json`
- `checkpoints.tsv`
- `evaluate_only_log.txt`

---

### Step 2: Report generation

```bash
python vaxlab_report/report_only.py \
  -i path/to/input.fasta \
  -o output_directory \
  --forna qbio  # Optional: choose Forna server (qbio/tbi)
```

**Generates:**
- `report.html`

---

## 🔧 Advanced Options

### IDT Complexity Scoring

To enable IDT gBlock complexity analysis, obtain an API token from IDT and use the `--token` parameter:

```bash
python vaxlab_report/evaluate_only.py \
  -i input.fasta \
  -o output/ \
  --token "YOUR_IDT_API_TOKEN"
```

**IDT Score Interpretation:**
- ≤5: Good 😊 (synthesis ready)
- 5-15: Caution ⚠️ (may need optimization)
- ≥15: Poor 😞 (difficult to synthesize)

### Forna Structure Visualization

Choose between different Forna servers for RNA structure visualization:

```bash
# Use qbio server (default when --forna is provided)
python vaxlab_report/report_only.py -i input.fasta -o output/ --forna qbio

# Use TBI server  
python vaxlab_report/report_only.py -i input.fasta -o output/ --forna tbi

# Disable structure visualization
python vaxlab_report/report_only.py -i input.fasta -o output/
```

---

## 🧪 Preset Configuration

Use a JSON file to customize metrics and evaluation parameters:

```json
{
  "global_metrics": ["gc", "cai", "mfe"],
  "local_metrics": ["gc", "degscore", "aup"],
  "fitness": {
    "cai": {"codon_table": "standard"},
    "degscore": {"structure_source": "vienna"}
  }
}
```

---

## 📂 Output Files

| File | Description |
|------|-------------|
| `evaluation_result.json` | Detailed metrics for each sequence |
| `checkpoints.tsv` | Summary of global metrics |
| `report.html` | Interactive report |
| `evaluate_only_log.txt` | Execution log |

---

## 🧾 License

Distributed under the [MIT License](LICENSE).
