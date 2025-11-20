# BioQuery NoLocal 🧬

**Cloud-Powered Natural Language Bioinformatics Tool to do EMBOSS Analyses**

A Streamlit web application that uses natural language AI (Google Gemini) to convert user requests into EMBOSS bioinformatics analyses. Perfect for researchers who want to run bioinformatics tools without learning command-line syntax or managing local AI models.

## Features

- 🤖 **Natural Language Interface**: Ask questions in plain English, get bioinformatics analysis
- 🧬 **EMBOSS Integration**: Access 258+ powerful EMBOSS tools (dynamically discovered)
- 🧪 **Gene-Based Queries**: Query genes by symbol (e.g., "translate ALKBH1") with Ensembl API integration
- 🔬 **BLAST Integration**: Run NCBI BLAST searches on sequences or genes remotely via BioPython
- 🧬 **BLAT Integration**: Search sequences against genomes using UCSC BLAT for near-exact matches
- 🧬 **BEDTools Support**: Find overlaps between genomic regions (SNPs, genes, regulatory elements)
- 📊 **Multi-Step Workflows**: Chain operations together (e.g., "find gene ALKBH1 then BLAST it")
- 🧬 **DNA/RNA Conversion**: Convert between DNA and RNA sequences (T↔U)
- 🎯 **AI-Powered**: Google Gemini API for intelligent tool selection and parameter extraction
- 🌐 **Web Interface**: Beautiful Streamlit UI with 5 integrated tabs
- 💾 **Results Export**: Download analysis results as text files

## Quick Start

### 1. Environment Setup

```bash
# Create conda environment
conda create -n bioquery python=3.9 -y
conda activate bioquery

# Install bioinformatics tools
conda install -c bioconda emboss biopython bedtools -y
conda install -c conda-forge streamlit pandas -y

# Install Python packages
pip install biopython google-generativeai requests
```

### 2. Get Google Gemini API Key

1. Go to https://ai.google.dev/
2. Click "Get API Key" in Google AI Studio
3. Create a new API key (it's free!)
4. Set the environment variable:

**On Linux/macOS:**
```bash
export GOOGLE_API_KEY='your_api_key_here'
```

**On Windows PowerShell:**
```powershell
$env:GOOGLE_API_KEY = "your_api_key_here"
```

**Or create a `.env` file:**
```bash
# Copy the example file
cp .env.example .env

# Edit .env and add your key
GOOGLE_API_KEY=your_api_key_here
```

### 3. Run the App

**Linux/macOS:**
```bash
conda activate bioquery
./run.sh
```

**Windows:**
```bash
conda activate bioquery
run.bat
```

Or manually:
```bash
streamlit run src/app.py
```

The app will open at `http://localhost:8501`

## Usage

### Method 1: Natural Language Query (Recommended)
1. Go to the **"🤖 Natural Language Query"** tab
2. Type what you want to do, e.g.:
   - "Translate this DNA to protein: ATGAAATTTCCC"
   - "What's the reverse complement of GCTA?"
   - "Find all open reading frames in ATGAAATTTCCCGGGAAATTT"
   - "Find gene info for ALKBH1"
   - "Translate ALKBH1 gene"
   - "Convert this DNA to RNA: ATGCCC"
   - "Find ALKBH1 gene then BLAST it"
3. Click **"🚀 Analyze"**
4. Results appear below, you can download them

### Multi-Step Queries
Chain multiple operations together using "then" or "and then":
- "Find gene ALKBH1, then calculate its GC content"
- "Get transcript sequence for TP53 and then BLAST it"
- "Translate BRCA1 then find restriction sites"

### Method 2: Manual Tool Selection
1. Go to the **"🔧 Manual Tool Selection"** tab
2. Choose a tool from the dropdown
3. Enter your sequence(s)
4. Click the tool button
5. Results appear as formatted output

## Available Tools

| Tool | Description | Input | Example Query |
|------|-------------|-------|---------------|
| **translate** | Translate DNA to protein | DNA sequence or gene name, reading frame (1-3) | "Translate ALKBH1" |
| **reverse** | Reverse complement DNA | DNA sequence | "Reverse complement ATGC" |
| **orf** | Find open reading frames | DNA sequence, min size | "Find ORFs in ATGCCC" |
| **align** | Needleman-Wunsch alignment | Two sequences | "Align ATGC and ATGG" |
| **restriction** | Find restriction sites | DNA sequence, enzyme (optional) | "Find restriction sites in ATGC" |
| **shuffle** | Randomize sequence | Sequence | "Shuffle ATGC" |
| **info** | Get sequence statistics | Sequence | "Get info for ATGC" |
| **sixframe** | All 6 translation frames | DNA sequence | "Six frame translate ATGC" |
| **cusp** | Codon usage statistics | DNA/CDS sequence | "Calculate codon usage for ALKBH1" |
| **pepstats** | Protein statistics (MW, pI, composition) | Protein sequence or gene | "Get protein stats for TP53" |
| **wordcount** | Oligonucleotide frequencies | DNA sequence | "Count oligonucleotides in sequence" |
| **gene_query** | Get gene information from Ensembl | Gene symbol (e.g., ALKBH1) | "Find gene info for TP53" |
| **blast** | Run NCBI BLAST search | DNA/protein sequence or gene name | "BLAST this sequence: ATGC" |
| **blat** | UCSC BLAT search (near-exact) | DNA sequence | "BLAT search ATGCCC in hg38" |
| **bedtools** | Find genomic overlaps | Two BED files or regions | "Find SNPs overlapping tRNA genes" |
| **dna_to_rna** | Convert DNA to RNA (T→U) | DNA sequence | "Convert ATGC to RNA" |
| **rna_to_dna** | Convert RNA to DNA (U→T) | RNA sequence | "Convert AUGC to DNA" |

**Plus 250+ additional EMBOSS tools** available through dynamic discovery!

## Project Structure

```
bme110/
├── src/
│   ├── app.py                 # Main Streamlit application
│   ├── emboss_wrapper.py      # EMBOSS tool wrapper & bioinformatics engine
│   └── nlp_handler.py         # Natural language processing handler
├── tests/
│   ├── test_e2e.py            # End-to-end integration tests
│   ├── test_gene_query.py     # Gene query tests
│   ├── test_gene_direct.py    # Direct gene lookup tests
│   ├── test_genome_query.py   # Genome analysis tests
│   └── test_nlp_gene.py       # NLP gene recognition tests
├── docs/
│   ├── ARCHITECTURE.md        # System design & diagrams (9 Mermaid diagrams)
│   ├── GETTING_STARTED.md     # Complete setup & installation guide
│   └── CONFIGURATION.md       # Configuration & troubleshooting
├── run.sh                     # Start script (Linux/macOS)
├── run.bat                    # Start script (Windows)
├── setup.sh                   # Environment setup script
├── setup_windows.ps1          # Environment setup script (Windows)
├── requirements.txt           # Python dependencies
├── LICENSE                    # Project license
├── .gitignore                 # Git ignore rules
└── README.md                  # This file
```

## How It Works

1. **User Input** → Streamlit UI collects natural language query or gene symbol
2. **NLP Processing** → Ollama LLM interprets the query and selects appropriate tool(s)
3. **Parameter Extraction** → LLM extracts sequence data, gene names, and parameters from query
4. **Gene Resolution** (if applicable) → Ensembl API fetches gene/transcript information
5. **EMBOSS Execution** → EMBOSSWrapper runs the selected EMBOSS tool(s) or BLAST
6. **Multi-Step Chaining** (if applicable) → Automatically chains results between steps
7. **Results Display** → Streamlit displays formatted results with download option

### Key Components

- **258+ EMBOSS Tools**: Dynamically discovered bioinformatics tools
- **Ensembl Integration**: Query genes by symbol, get transcript sequences
- **BioPython BLAST**: Remote NCBI BLAST searches (Bio.Blast.NCBIWWW, Bio.Blast.NCBIXML)
- **Multi-Step Execution**: Chain operations with automatic result passing
- **Smart Parameter Handling**: Automatically fetches sequences when gene names are provided

## System Requirements

- **OS**: Windows (WSL2), macOS, or Linux
- **Python**: 3.9+
- **RAM**: 4GB minimum (for model loading)
- **Disk**: 5GB (for models and EMBOSS tools)
- **Internet**: For initial model download

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│          Streamlit Web Interface (src/app.py)                   │
│                     5 Tabs: NLP | Manual | Genome | Batch | Docs│
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────────────┐  ┌─────────────────────────────────┐ │
│  │  NLP Handler         │  │ EMBOSS Wrapper                  │ │
│  │  (nlp_handler.py)    │  │ (emboss_wrapper.py)             │ │
│  │                      │  │                                 │ │
│  │  • Query Parsing     │  │ • 258+ EMBOSS tools             │ │
│  │  • Tool Selection    │  │ • Gene query (Ensembl)          │ │
│  │  • Multi-Step        │  │ • BLAST integration             │ │
│  │  • Param Extract     │  │ • DNA/RNA conversion            │ │
│  └────────┬─────────────┘  │ • Multi-step execution          │ │
│           │                │ • Smart parameter chaining      │ │
│           └────┬──────────►└────┬────────────────────────────┘ │
└────────────────┼───────────────┼──────────────────────────────┘
                 │               │
        ┌────────▼─────┐  ┌──────▼────────────┐
        │   Ollama     │  │  External APIs    │
        │ (gemma3:4b)  │  │  • EMBOSS tools   │
        │              │  │  • Ensembl API    │
        │ • LLM magic  │  │  • NCBI BLAST     │
        └──────────────┘  │  • BioPython      │
                          └───────────────────┘
```

### Data Flow Example: Gene-Based BLAST
```
User: "Find ALKBH1 gene then BLAST it"
  ↓
NLP Handler: Detects multi-step query
  Step 1: gene_query (gene_name: ALKBH1)
  Step 2: blast (sequence: from_previous_step)
  ↓
EMBOSS Wrapper:
  1. Queries Ensembl for ALKBH1 → gets transcript ID
  2. Fetches transcript sequence → 2597 bp
  3. Passes sequence to BLAST → NCBI remote search
  4. Parses BLAST XML results → top alignments
  ↓
Display: Formatted BLAST results with download option
```

## Troubleshooting

### Ollama Connection Issues

**Error: "Failed to connect to Ollama"**
- Make sure Ollama is running: `ollama serve` on Windows PowerShell
- Check the host IP: On WSL2, use `ip route | grep default | awk '{print $3}'`
- Update `nlp_handler.py` line 13 with correct IP if needed

### EMBOSS Not Found

**Error: "EMBOSS not found"**
```bash
conda activate bioquery
conda install -c bioconda emboss -y
```

### Model Too Large

If your model is taking too long to download:
```bash
ollama pull phi3:mini  # Smallest model (2GB)
```

## Performance Notes

- **First run**: May be slow while loading model into memory (~1-2 minutes)
- **Subsequent runs**: Much faster as model stays loaded
- **Sequence size**: Works well with sequences up to 10MB
- **Alignment**: Large sequence alignment may take time

## Contributing

Feel free to:
- Add more EMBOSS tool wrappers in `emboss_wrapper.py`
- Improve NLP parsing in `nlp_handler.py`
- Enhance UI in `app.py`
- Add new features or bug fixes

## License

This is an educational project for BME110 at UCSC.

## References

- [EMBOSS Documentation](http://emboss.open-bio.org/)
- [BioPython Documentation](https://biopython.org/)
- [Ollama Documentation](https://ollama.com)
- [Streamlit Documentation](https://docs.streamlit.io/)

---

**Built with:** EMBOSS | BioPython | Ollama | Streamlit

**Author:** Cagn Steinbrecher  
**Course:** BME110  
**Institution:** UC Santa Cruz
