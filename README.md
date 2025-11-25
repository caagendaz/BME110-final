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
- 🧬 **GTEx Integration**: Query gene expression across 54 human tissues from GTEx Portal
- 📊 **Multi-Step Workflows**: Chain operations together (e.g., "find gene ALKBH1 then BLAST it")
- 🧬 **DNA/RNA Conversion**: Convert between DNA and RNA sequences (T↔U)
- 📋 **Command Logging**: Detailed execution logs for debugging and documentation (NEW!)
- 🛡️ **Enhanced Safety Filters**: Improved handling of scientific terminology in AI queries
- 🎯 **AI-Powered**: Google Gemini API for intelligent tool selection and parameter extraction
- 🌐 **Web Interface**: Beautiful Streamlit UI with 6 integrated tabs
- 💾 **Results Export**: Download analysis results and execution logs as text files

## Quick Start

### 1. Environment Setup

```bash
# Create conda environment with Python 3.12
conda create -n bioquery python=3.12 -y
conda activate bioquery

# Install core packages first (conda-forge)
conda install -c conda-forge streamlit pandas -y

# Then install bioinformatics tools (bioconda)
conda install -c bioconda emboss bedtools -y

# Install Python API packages
pip install biopython google-generativeai requests
```

**Note**: Installing conda-forge packages before bioconda avoids dependency conflicts.

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

### Multiple Question Mode
Paste entire assignments with numbered questions - each will be processed independently:
```
1. What is the GC content of SOCS3?
2. Translate ALKBH1 gene
3. What tissues express TP53 highest?
```

Supports formats:
- `1.`, `2.`, `3.` (numbered list)
- `Question 1:`, `Question 2:` (labeled format)

Each question gets its own result section with expand/collapse capability.

### Method 2: Manual Tool Selection
1. Go to the **"🔧 Manual Tool Selection"** tab
2. Choose a tool from the dropdown
3. Enter your sequence(s)
4. Click the tool button
5. Results appear as formatted output

### Method 3: Command Log (Debugging & Documentation)
1. Go to the **"📋 Command Log"** tab
2. View all executed commands with:
   - Timestamp and tool name
   - All parameters used (with smart truncation for long sequences)
   - Result previews
   - Success/failure status
   - Error messages (if any)
3. Download log as text file for assignment documentation
4. Clear log to start fresh session

**Use this for:**
- Verifying what commands were actually executed
- Debugging unexpected results
- Documenting your workflow for assignments
- Understanding multi-step operations

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
| **gtex** | Get tissue expression data | Gene symbol | "What tissues express SOCS3 highest?" |
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
2. **NLP Processing** → Google Gemini API (gemini-2.5-flash) interprets the query and selects appropriate tool(s)
   - Uses academic/scientific context wrapper to bypass false-positive safety filters
3. **Parameter Extraction** → Gemini extracts sequence data, gene names, and parameters from query
4. **Gene Resolution** (if applicable) → Ensembl API fetches gene/transcript information
5. **Tool Execution** → EMBOSSWrapper runs EMBOSS tools, BLAST, BEDTools, or other APIs
   - **NEW**: All executions are automatically logged with timestamp, parameters, and results
6. **Multi-Step Chaining** (if applicable) → Automatically chains results between steps
7. **Results Display** → Streamlit displays formatted results with download option
8. **Command Logging** → View detailed execution history in Command Log tab

### Key Components

- **258+ EMBOSS Tools**: Dynamically discovered bioinformatics tools
- **Ensembl Integration**: Query genes by symbol, get transcript sequences
- **BioPython BLAST**: Remote NCBI BLAST searches (Bio.Blast.NCBIWWW, Bio.Blast.NCBIXML)
- **Multi-Step Execution**: Chain operations with automatic result passing
- **Smart Parameter Handling**: Automatically fetches sequences when gene names are provided
- **Command Logging System**: Tracks all executions for debugging and documentation

## System Requirements

- **OS**: Windows (WSL2), macOS, or Linux
- **Python**: 3.12+ (tested on 3.12.12)
- **RAM**: 2GB minimum
- **Disk**: 2GB (for EMBOSS tools and packages)
- **Internet**: Required for Google Gemini API, NCBI BLAST, Ensembl, GTEx, UCSC APIs

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│          Streamlit Web Interface (src/app.py)                   │
│         6 Tabs: NLP | Manual | Genome | Batch | Log | Docs     │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────────────┐  ┌─────────────────────────────────┐ │
│  │  NLP Handler         │  │ EMBOSS Wrapper                  │ │
│  │  (nlp_handler.py)    │  │ (emboss_wrapper.py)             │ │
│  │                      │  │                                 │ │
│  │  • Query Parsing     │  │ • 258+ EMBOSS tools             │ │
│  │  • Tool Selection    │  │ • Gene query (Ensembl)          │ │
│  │  • Multi-Question    │  │ • BLAST integration (BioPython) │ │
│  │  • Multi-Step        │  │ • BLAT search (UCSC)            │ │
│  │  • Param Extract     │  │ • BEDTools (genomic overlaps)   │ │
│  │  • Safety Wrapper    │  │ • GTEx expression links         │ │
│  └────────┬─────────────┘  │ • DNA/RNA conversion            │ │
│           │                │ • Protein analysis (pepstats)   │ │
│           │                │ • Command logging (NEW!)        │ │
│           └────┬──────────►└────┬────────────────────────────┘ │
└────────────────┼───────────────┼──────────────────────────────┘
                 │               │
        ┌────────▼─────────┐  ┌──▼──────────────────────────┐
        │  Google Gemini   │  │  External APIs & Tools      │
        │ gemini-2.5-flash │  │  • EMBOSS (local)           │
        │                  │  │  • BEDTools (local)         │
        │ • NLP parsing    │  │  • Ensembl REST API         │
        │ • Tool routing   │  │  • NCBI BLAST (remote)      │
        │ • Multi-question │  │  • UCSC Genome Browser API  │
        │ • Safety bypass  │  │  • GTEx Portal (links)      │
        └──────────────────┘  └─────────────────────────────┘
```

### Data Flow Example: Gene-Based BLAST
```
User: "Find ALKBH1 gene then BLAST it"
  ↓
NLP Handler (Gemini API):
  Detects multi-step query
  Step 1: gene_query (gene_name: ALKBH1)
  Step 2: blast (gene_name: ALKBH1) - auto-resolves sequence
  ↓
EMBOSS Wrapper:
  1. Queries Ensembl REST API for ALKBH1 → gets transcript ID
  2. Fetches transcript sequence → 2597 bp CDS
  3. Submits to NCBI BLAST via BioPython → remote search
  4. Parses BLAST XML results → top alignments with E-values
  ↓
Streamlit Display:
  Formatted BLAST results with:
  • Hit descriptions
  • E-values and bit scores
  • Identity percentages
  • Download button for results
```

## Troubleshooting

### Google API Key Issues

**Error: "Google API key required"**
- Get a free API key from https://ai.google.dev/
- Set it in your environment: `export GOOGLE_API_KEY='your_key'`
- Or create a `.env` file with `GOOGLE_API_KEY=your_key`

**Error: "Rate limit exceeded"**
- Gemini free tier: 10 requests/minute, 1500 requests/day
- Wait ~60 seconds before retrying
- For multi-question queries, use smaller batches (3-5 questions at a time)

### EMBOSS Not Found

**Error: "EMBOSS not found"**
```bash
conda activate bioquery
conda install -c bioconda emboss -y
```

### BEDTools Not Found

**Error: "bedtools: command not found"**
```bash
conda activate bioquery
conda install -c bioconda bedtools -y
```

## Performance Notes

- **NLP queries**: ~1-3 seconds per question (Gemini API)
- **BLAST searches**: 10-60 seconds (NCBI remote server)
- **Local EMBOSS tools**: <1 second for most operations
- **Gene queries**: 1-2 seconds (Ensembl API)
- **Multi-question mode**: Process N questions in ~N seconds + rate limiting
- **Rate limits**: Gemini free tier = 10 requests/minute
- **Sequence size**: Works well with sequences up to 10MB
- **Alignment**: Large sequence alignment may take additional time

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
- [Google Gemini API Documentation](https://ai.google.dev/)
- [Streamlit Documentation](https://docs.streamlit.io/)
- [BEDTools Documentation](https://bedtools.readthedocs.io/)
- [GTEx Portal](https://gtexportal.org/)
- [UCSC Genome Browser API](https://genome.ucsc.edu/goldenPath/help/api.html)

---

**Built with:** EMBOSS | BioPython | Google Gemini | Streamlit | BEDTools

**Author:** Cagn Steinbrecher  
**Course:** BME110  
**Institution:** UC Santa Cruz  
**Version:** 2.0 (Python 3.12, Gemini 2.5 Flash)
