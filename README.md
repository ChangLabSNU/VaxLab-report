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

### Install from Source

```bash
# 1. Clone the repository
git clone https://github.com/ChangLabSNU/VaxLab-report.git
cd vaxlab-report

# 2. Set up the environment (Python 3.9+)
conda create -y -n vaxlab_report python=3.9
conda activate vaxlab_report

# 3. Install ViennaRNA
conda install -y -c bioconda viennarna

# 4. Install VaxLab-report and dependencies
pip install -e .
```


**Dependencies:**
- Python 3.9+
- ViennaRNA (for RNA secondary structure prediction)
- BioPython (for sequence handling)
- Jinja2 (for HTML report templating)
- Requests (for IDT API integration)
- NumPy (for numerical computations)

---

## 🚀 Usage

### Step 1: Sequence evaluation

```bash
python vaxlab_report/evaluate_only.py \
  -i path/to/input.fasta \
  -o output_directory \
  --preset path/to/parameters.json \
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

## 🧪 Testing

### Regression Tests

Run automated regression tests to ensure the pipeline works correctly:

```bash
python test_regression.py
```

**Test Coverage:**
- ✅ `evaluate_only.py` - CDS and mRNA sequence evaluation
- ✅ `report_only.py` - HTML report generation  
- ✅ Forna visualization - Both qbio and tbi server options
- ✅ File output validation - Ensures all expected files are created
- ✅ JSON structure validation - Checks evaluation result format

**Test Input:** Uses `test_data.fasta` (HA sequence with 5'UTR and 3'UTR)

**Expected Output:**
```
🧪 Running VaxLab-report regression tests...

Testing evaluate_only.py...
✅ evaluate_only.py test passed

Testing report_only.py...
✅ report_only.py test passed

Testing report_only.py with --forna option...
✅ Forna options test passed

📊 Test Results: 3/3 tests passed
🎉 All regression tests passed!
```

### Manual Testing

For manual testing with your own data:

```bash
# Test evaluation only
python vaxlab_report/evaluate_only.py -i your_file.fasta -o test_output/

# Test with IDT token
python vaxlab_report/evaluate_only.py -i your_file.fasta -o test_output/ --token "YOUR_TOKEN"

# Test report generation
python vaxlab_report/report_only.py -i your_file.fasta -o test_output/ --forna qbio
```

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

**IDT Complexity Score Interpretation:**
- <7: Low 😊 (synthesis ready)
- 7-20: Moderate ⚠️ (may need optimization)
- ≥20: High 😞 (difficult to synthesize)

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

## 🔧 Troubleshooting

### Common Issues

**Import Error: No module named 'RNA'**
```bash
# Install ViennaRNA
conda install -c bioconda viennarna
```

**IDT API Authentication Failed**
- Check your IDT API token is valid
- Obtain a new token from IDT if expired
- Re-run evaluation with `--token "YOUR_NEW_TOKEN"`

**Structure Visualization Not Working**
- Use `--forna qbio` or `--forna tbi` to enable visualization
- Check internet connectivity for Forna server access

**Test Failures**
```bash
# Run tests with verbose output
python test_regression.py 2>&1 | tee test_log.txt
```

---

## 🧾 License

Distributed under the [MIT License](LICENSE).
