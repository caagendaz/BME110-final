# BioQuery NoLocal 🧬

**Cloud-Powered Natural Language Bioinformatics Tool to do EMBOSS Analyses**

A Streamlit web application that uses natural language AI (Google Gemini) to convert user requests into EMBOSS bioinformatics analyses. Perfect for researchers who want to run bioinformatics tools without learning command-line syntax or managing local AI models.

## Features

- 🤖 **Natural Language Interface**: Ask questions in plain English, get bioinformatics analysis
- 🧬 **258+ EMBOSS Tools**: Access entire EMBOSS suite with AI-powered dynamic discovery ⚡NEW!
- 🎯 **Smart Tool Resolution**: AI automatically finds the right tool even if not explicitly mapped ⚡NEW!
- 🏠 **Local or Cloud AI**: Choose between Gemini (cloud) or Ollama (fully local) for tool resolution ⚡NEW!
- 🚀 **Optimized Performance**: Sequences processed separately from NLP for 10x faster queries ⚡NEW!
- 📊 **Enhanced Graphics**: SVG output for crisp, scalable dotplots and visualizations ⚡NEW!
- 🧪 **Gene-Based Queries**: Query genes by symbol (e.g., "translate ALKBH1") with Ensembl API integration
- 🔬 **BLAST Integration**: Run NCBI BLAST searches on sequences or genes remotely via BioPython
- 🧬 **BLAT Integration**: Search sequences against genomes using UCSC BLAT for near-exact matches
- 🧬 **BEDTools Support**: Find overlaps between genomic regions (SNPs, genes, regulatory elements)
- 🧬 **GTEx Integration**: Query gene expression across 54 human tissues from GTEx Portal
- 📊 **Multi-Step Workflows**: Chain operations together (e.g., "find gene ALKBH1 then BLAST it")
- 🧬 **DNA/RNA Conversion**: Convert between DNA and RNA sequences (T↔U)
- 🧪 **18 Biopython Tools**: Advanced analysis including phylogenetic trees, motif finding, primer design ⚡NEW!
- 📋 **Command Logging**: Detailed execution logs for debugging and documentation
- 🛡️ **Enhanced Safety Filters**: Improved handling of scientific terminology in AI queries
- 🎯 **AI-Powered**: Google Gemini 2.0 Flash OR Ollama (llama3.2, qwen2.5, mistral) ⚡NEW!
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

### 2. Choose Your AI Backend

#### Option A: Cloud Mode (Gemini) - Recommended for Best Accuracy

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

#### Option B: Local Mode (Ollama) - Fully Private, No API Keys ⚡NEW!

For complete privacy with no internet required for AI queries:

1. **Install Ollama**: Download from https://ollama.ai
2. **Pull a model** (one-time):
   ```bash
   ollama pull llama3.2      # Recommended: 2GB, fast, accurate
   # OR
   ollama pull qwen2.5       # Alternative: good for technical tasks
   # OR  
   ollama pull mistral       # Alternative: 4GB, very capable
   ```
3. **Ollama runs automatically** on `localhost:11434`
4. **Select "Local Mode"** in the app sidebar

**Ollama Mode Features:**
- ✅ Complete privacy - all AI processing stays on your machine
- ✅ No API keys required
- ✅ No internet needed for tool resolution
- ✅ Works with llama3.2, llama3.1, qwen2.5, mistral, and more
- ⚠️ Slightly less accurate than Gemini for complex queries

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

### EMBOSS Tools (258+ available via AI-powered discovery)
The system can dynamically discover and use any EMBOSS tool. Here are commonly used ones:

| Tool | Description | Input | Example Query |
|------|-------------|-------|---------------|
| **translate** | Translate DNA to protein | DNA sequence or gene name | "Translate ALKBH1" |
| **reverse** | Reverse complement DNA | DNA sequence | "Reverse complement ATGC" |
| **orf** | Find open reading frames | DNA sequence, min size | "Find ORFs in ATGCCC" |
| **align/needle** | Needleman-Wunsch global alignment | Two sequences | "Align ATGC and ATGG" |
| **water** | Smith-Waterman local alignment | Two sequences | "Use water to align these" |
| **stretcher** | Fast global alignment ⚡NEW | Two sequences | "Use stretcher on sequences" |
| **matcher** | Find all local matches ⚡NEW | Two sequences | "Use matcher to compare" |
| **dotplot** | Dot plot (SVG format) ⚡NEW | Two sequences | "Make a dot plot" |
| **newcpgseek** | Find CpG islands ⚡NEW | DNA sequence | "Find CpG islands" |
| **pepnet** | Helical net plot ⚡NEW | Protein sequence | "Use pepnet to visualize" |
| **antigenic** | Predict antigenic sites ⚡NEW | Protein sequence | "Find antigenic sites" |
| **restriction** | Find restriction sites | DNA sequence | "Find restriction sites" |
| **sixframe** | All 6 reading frames | DNA sequence | "Six frame translate" |
| **cusp** | Codon usage statistics | DNA/CDS | "Calculate codon usage" |
| **pepstats** | Protein statistics | Protein sequence | "Get protein stats" |
| **wossname** | Search EMBOSS docs ⚡NEW | Search term | "Find alignment tools" |
| **showdb** | List databases ⚡NEW | None | "Show available databases" |

### Biopython Tools (18 specialized analyses)
Advanced bioinformatics analysis powered by Biopython:

| Tool | Description | Input | Example Query |
|------|-------------|-------|---------------|
| **phylo** | Phylogenetic tree (UPGMA) | Multiple sequences | "Build phylogenetic tree" |
| **motif** | Find conserved motifs | Multiple sequences | "Find motifs in sequences" |
| **primer** | Primer analysis (Tm, GC%) | DNA primer | "Analyze this primer" |
| **restriction_batch** | Batch restriction analysis | DNA sequence | "Find all restriction sites" |
| **pairwise_compare** | Detailed comparison | Two sequences | "Compare sequences in detail" |
| **secondary_structure** | Predict 2° structure | Protein sequence | "Predict secondary structure" |
| **hydrophobicity** | Hydrophobicity plot | Protein sequence | "Plot hydrophobicity" |
| **entropy** | Shannon entropy | Sequence | "Calculate entropy" |
| **complexity** | Linguistic complexity | Sequence | "Analyze complexity" |

### External API Integration
| Tool | Description | Input | Example Query |
|------|-------------|-------|---------------|
| **gene_query** | Ensembl gene info | Gene symbol | "Find gene info for TP53" |
| **blast** | NCBI BLAST search | DNA/protein or gene | "BLAST this sequence" |
| **blat** | UCSC BLAT (near-exact) | DNA sequence | "BLAT search in hg38" |
| **bedtools** | Genomic region overlaps | BED regions | "Find overlapping SNPs" |
| **gtex** | Tissue expression data | Gene symbol | "Which tissues express SOCS3?" |

⚡ **NEW**: Over 240 additional EMBOSS tools available through AI-powered dynamic discovery!

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
2. **Sequence Extraction** ⚡NEW → System extracts sequences from query BEFORE NLP processing (10x faster)
3. **NLP Processing** → Google Gemini API (gemini-2.0-flash-exp) interprets query and selects tool(s)
   - Uses academic/scientific context wrapper to bypass false-positive safety filters
   - **AI Tool Resolution** ⚡NEW: If tool not in map, asks Gemini "what EMBOSS tool is X?" and caches result
4. **Parameter Extraction** → Gemini extracts parameters from query (without sequences)
5. **Sequence Injection** ⚡NEW → Extracted sequences injected into parameters AFTER tool selection
6. **Gene Resolution** (if applicable) → Ensembl API fetches gene/transcript information
7. **Tool Execution** → EMBOSSWrapper runs tools with support for:
   - Single-sequence tools (`sequence` parameter)
   - Two-sequence tools (`seq1`, `seq2` parameters) ⚡NEW
   - No-sequence tools (`wossname`, `showdb`) ⚡NEW
   - Graphics tools (SVG output for quality) ⚡NEW
8. **Multi-Step Chaining** (if applicable) → Automatically chains results between steps
9. **Results Display** → Streamlit displays formatted results with download option
10. **Command Logging** → View detailed execution history in Command Log tab

### Key Components

- **258+ EMBOSS Tools**: Full suite dynamically discovered via AI
- **18 Biopython Tools**: Phylogenetics, motif finding, primers, restriction analysis, protein analysis
- **AI-Powered Tool Resolution**: Discovers unmapped tools on-the-fly via Gemini ⚡NEW
- **Optimized NLP Flow**: Sequences processed separately for speed ⚡NEW
- **Multi-File Upload**: Process multiple FASTA files, automatic sequence detection
- **Two-Sequence Comparison**: Intelligent routing for alignment/dotplot tools ⚡NEW
- **SVG Graphics**: High-quality scalable vector graphics ⚡NEW
- **Ensembl Integration**: Query genes by symbol, get transcript sequences
- **BioPython BLAST**: Remote NCBI BLAST searches
- **Multi-Step Execution**: Chain operations with automatic result passing
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

MIT License

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
