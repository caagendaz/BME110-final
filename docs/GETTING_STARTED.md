# 🧬 BioQuery Local - Complete Project Guide

## 📋 Table of Contents

1. [Project Overview](#project-overview)
2. [What We Built](#what-we-built)
3. [File Structure](#file-structure)
4. [How to Get Started](#how-to-get-started)
5. [Usage Guide](#usage-guide)
6. [GitHub Deployment](#github-deployment)
7. [Next Steps](#next-steps)

---

## Project Overview

**BioQuery Local** is an AI-driven bioinformatics analysis tool that lets researchers:
- Ask questions in plain English about DNA/protein sequences
- Get automatic tool selection and analysis
- See results instantly in a web browser
- All running **locally** without internet dependency

### Key Innovation
Instead of learning EMBOSS command syntax, users can simply ask:
- "Translate this DNA: ATGAAATTT"
- "Find ORFs in this sequence"
- "What's the reverse complement?"

The AI (Ollama LLM) understands and runs the right tool automatically.

---

## What We Built

### ✅ Core Components

#### 1. **EMBOSS Wrapper** (`src/emboss_wrapper.py` - 1334 lines)
Comprehensive Python interface to bioinformatics tools:
- **258+ EMBOSS tools** dynamically discovered
- **Gene-based operations** with Ensembl API integration
- **BLAST integration** via BioPython (NCBI remote searches)
- **Multi-step execution** with automatic result chaining
- **DNA/RNA conversion** (T↔U)
- Smart parameter handling and sequence fetching

**Key Features:**
- Translate DNA → Protein (with gene name support)
- Reverse & Complement DNA
- Find Open Reading Frames (ORFs)
- Global Sequence Alignment
- Restriction Site Finding
- Gene queries (Ensembl API)
- BLAST searches (NCBI)
- Multi-step workflows
- And 250+ more EMBOSS tools!

#### 2. **NLP Handler** (`src/nlp_handler.py` - 402 lines)
Converts natural language to bioinformatics operations:
- Ollama integration for local AI
- **Multi-step query support** (chain operations with "then")
- JSON-based tool selection
- Parameter extraction from natural language
- **Gene name recognition** (ALKBH1, TP53, BRCA1, etc.)
- Smart disambiguation (e.g., "translate to RNA" vs "translate to protein")
- Sequence validation

#### 3. **Web Interface** (`src/app.py` - 540 lines)
Beautiful Streamlit application with **5 integrated tabs**:
- **Natural language query interface** (AI-powered with multi-step support)
- **Manual tool selection** (user-driven)
- **Genome browser query** (genomic regions)
- **Batch analysis** (FASTA file upload)
- **Documentation & examples**
- Real-time results with formatted output
- Download functionality
- Multi-step workflow visualization

#### 4. **Comprehensive Testing** (`tests/` directory)
- End-to-end integration tests
- Gene query tests
- Direct gene lookup tests
- Genome analysis tests
- NLP gene recognition tests

---

## File Structure

```
bme110/
│
├── CORE APPLICATION FILES (src/)
│   ├── app.py                 # Main Streamlit web app (5 tabs)
│   ├── emboss_wrapper.py      # Bioinformatics engine (258+ tools)
│   └── nlp_handler.py         # NLP + Ollama integration
│
├── TESTING SUITE (tests/)
│   ├── test_e2e.py            # End-to-end integration tests
│   ├── test_gene_query.py     # Gene query tests
│   ├── test_gene_direct.py    # Direct gene lookup tests
│   ├── test_genome_query.py   # Genome analysis tests
│   └── test_nlp_gene.py       # NLP gene recognition tests
│
├── DOCUMENTATION (docs/)
│   ├── GETTING_STARTED.md     # This file - complete setup guide
│   ├── CONFIGURATION.md       # Configuration & troubleshooting
│   └── ARCHITECTURE.md        # System design (9 Mermaid diagrams)
│
├── SETUP & DEPENDENCIES
│   ├── requirements.txt       # Python package versions
│   ├── run.sh                 # Start script (Linux/macOS)
│   ├── run.bat                # Start script (Windows)
│   ├── setup.sh               # Linux/macOS automated setup
│   └── setup_windows.ps1      # Windows PowerShell setup
│
├── GIT & VERSION CONTROL
│   └── .gitignore             # Files to ignore in git
│
└── MISC
    ├── bme110final.py         # (empty, for reference)
    └── BME110-AI-Alternative-Project-Tool Presentation.pdf  # Original project PDF
```

**Total: 1,260+ lines of production code**

---

## How to Get Started

### Option 1: Fresh Setup (Recommended)

```bash
# 1. Clone the repository (or navigate to directory)
cd /mnt/c/Users/cagns/OneDrive/Desktop/bme110

# 2. Run setup script
bash setup.sh                    # Linux/macOS
# OR
powershell -ExecutionPolicy Bypass -File setup_windows.ps1  # Windows

# 3. Start Ollama (in a new terminal)
# Windows PowerShell:
$env:OLLAMA_HOST = "0.0.0.0:11434"
ollama serve

# 4. Pull a model (in another terminal)
ollama pull gemma3:4b            # Already pulled, but can re-pull other models

# 5. Run the app
conda activate bioquery
streamlit run app.py

# 6. Open browser to http://localhost:8501
```

### Option 2: Manual Setup

```bash
# Create environment
conda create -n bioquery python=3.9 -y
conda activate bioquery

# Install bioinformatics tools
conda install -c bioconda emboss biopython -y
conda install -c conda-forge streamlit pandas -y

# Install Python packages
pip install biopython ollama

# Run the app
streamlit run app.py
```

### Option 3: Quick Test (Everything Already Set Up)

```bash
# Assuming bioquery environment exists and Ollama is running
conda activate bioquery
streamlit run app.py
```

---

## Usage Guide

### Interface Tab 1: Natural Language Query (Recommended)

```
User: "Translate this DNA to protein: ATGAAATTTCCCGGG"
                ↓
AI parses the query
                ↓
Selects: tool=translate, sequence=ATGAAATTTCCCGGG
                ↓
Runs EMBOSS transeq
                ↓
Shows results: "MKF PG" (protein sequence)
```

**Examples you can try:**
- "Find ORFs in ATGAAATTTCCCGGGAAATTTAAAGGG"
- "What's the reverse complement of GCTA?"
- "Align ATGAAA with ATGCCC"
- "Give me info about this sequence: ATGAAATTTCCCGGG"
- "Show all 6 reading frames for ATGAAATTTCCCGGGAAATTTAAAGGG"
- **"Find gene info for ALKBH1"** (Gene queries!)
- **"Translate ALKBH1 gene"** (Gene-based operations!)
- **"Convert this DNA to RNA: ATGCCC"** (DNA/RNA conversion!)
- **"Find ALKBH1 gene then BLAST it"** (Multi-step workflows!)

### Interface Tab 2: Manual Tool Selection

1. Choose a tool from dropdown
2. Paste your sequence(s) **or enter a gene name**
3. Adjust parameters if needed
4. Click the tool button
5. Results appear below

**Good for:** When you know exactly which tool you need

### Interface Tab 3: Genome Browser Query

1. Enter genomic coordinates (e.g., chr1:1000-2000)
2. Select genome assembly (hg38, hg19, etc.)
3. Query UCSC Genome Browser API
4. Get region information

**Good for:** Exploring genomic regions

### Interface Tab 4: Batch Analysis

1. Upload a FASTA file
2. Analyze multiple sequences at once
3. Compare results

**Good for:** Processing multiple sequences together

### Interface Tab 5: Documentation

- Complete feature overview
- Example queries
- Technology explanation
- Troubleshooting links

---

## GitHub Deployment

### Repository Already Set Up! ✅

This project is already on GitHub at:
**https://github.com/caagendaz/BME110-final**

The repository includes:
- ✅ All source code in `src/` directory
- ✅ Complete test suite in `tests/` directory
- ✅ Documentation in `docs/` directory
- ✅ Setup scripts (`run.sh`, `run.bat`)
- ✅ Requirements and configuration files

### To Clone and Use:

```bash
# Clone the repository
git clone https://github.com/caagendaz/BME110-final.git
cd BME110-final

# Follow setup instructions in README.md
conda create -n bioquery python=3.9 -y
conda activate bioquery
conda install -c bioconda emboss biopython -y
conda install -c conda-forge streamlit pandas -y
pip install ollama

# Run the application
./run.sh  # Linux/macOS
# or
run.bat   # Windows
```

---

## Next Steps

### For Demonstration

- [x] ✅ App fully functional with 258+ tools
- [x] ✅ Natural language queries working
- [x] ✅ Gene-based operations (ALKBH1, TP53, etc.)
- [x] ✅ BLAST integration (NCBI remote)
- [x] ✅ Multi-step workflows
- [x] ✅ DNA/RNA conversion
- [x] ✅ Repository organized and on GitHub

**Demo Examples to Show:**
1. Simple query: "Translate this DNA: ATGAAATTT"
2. Gene query: "Find gene info for ALKBH1"
3. Gene-based tool: "Translate TP53 gene"
4. Multi-step: "Find ALKBH1 then BLAST it"
5. DNA/RNA conversion: "Convert ATGCCC to RNA"

### Current Features ✅

1. **258+ EMBOSS Tools** - All discovered dynamically
2. **Gene Queries** - Ensembl API integration
3. **BLAST Integration** - Remote NCBI searches via BioPython
4. **Multi-Step Workflows** - Chain operations with automatic result passing
5. **DNA/RNA Conversion** - Simple T↔U conversion
6. **Smart NLP** - Handles gene names, ambiguous queries, multi-step commands
7. **5-Tab Interface** - NLP, Manual, Genome, Batch, Documentation

### Future Enhancements (Optional)

1. **Enhanced Visualization**
   - Sequence alignment plots
   - Gene structure diagrams
   - BLAST hit visualization

2. **Additional Data Sources**
   - UniProt protein data
   - PDB structure information
   - More genome assemblies

3. **Advanced Features**
   - Custom pipelines
   - Result comparison tools
   - Export to different formats
   - GC content graphs

3. **Advanced Features** (Week 4+)
   - Analysis history
   - Saved workflows
   - Collaboration features

---

## Quick Reference

### Starting the App

```bash
conda activate bioquery
streamlit run app.py
# Open: http://localhost:8501
```

### Common Issues & Fixes

| Issue | Fix |
|-------|-----|
| "Connection refused" | Make sure Ollama is running: `ollama serve` |
| "EMBOSS not found" | Run: `conda install -c bioconda emboss -y` |
| "Model too large" | Use smaller model: `ollama pull phi3:mini` |
| "Port 8501 in use" | Kill Streamlit: `pkill streamlit` or use different port |
| "Out of memory" | Close other apps or use smaller model |

### Useful Commands

```bash
# Check Ollama models
ollama list

# Test EMBOSS
embossversion

# Check conda environment
conda info --envs

# Activate environment
conda activate bioquery

# View logs
# Windows: Check TaskManager for Ollama process
# Linux/macOS: tail -f ~/.ollama/logs/server.log

# Stop Streamlit
pkill streamlit

# Test connection (from WSL2)
python -c "from ollama import Client; c = Client(host='http://192.168.128.1:11434'); print(c.list())"
```

---

## Project Statistics

- **Lines of Code**: 1,260+
- **Functions Implemented**: 25+
- **EMBOSS Tools Wrapped**: 13+
- **Documentation Pages**: 5
- **Setup Scripts**: 2 (bash + PowerShell)
- **Development Time**: ~2 hours (setup) + testing

---

## Technology Stack

```
Frontend:
  └─ Streamlit 1.48.1
     ├─ Custom CSS styling
     ├─ File upload handling
     └─ Results export

Backend:
  ├─ Python 3.9
  ├─ EMBOSS 6.5.7 (bioinformatics)
  ├─ BioPython 1.85
  ├─ Pandas 2.3.1
  └─ Ollama 0.6.0 (AI)
     └─ Gemma3-4B model

Integration:
  ├─ Subprocess for EMBOSS
  ├─ Tempfile for I/O
  ├─ JSON for configuration
  └─ HTTP for Ollama connection
```

---

## Success Criteria ✅

- [x] Natural language query processing works
- [x] EMBOSS integration functional
- [x] Ollama LLM successfully selects tools
- [x] Streamlit web interface responsive
- [x] Results export working
- [x] Cross-platform support (Windows, Linux, macOS)
- [x] Comprehensive documentation
- [x] Setup automated
- [x] Ready for GitHub deployment
- [x] Ready for demonstration

---

## Final Checklist Before Submission

- [ ] Test all features in the web app
- [ ] Verify all EMBOSS tools produce correct output
- [ ] Confirm Ollama connection working
- [ ] Check documentation is complete
- [ ] Initialize GitHub repository
- [ ] Push all code to GitHub
- [ ] Create a short demo script (optional)
- [ ] Prepare presentation materials

---

## Contact & Support

For issues or questions:
1. Check **CONFIGURATION.md** for advanced setup
2. See **README.md** for feature overview
3. Review **PROJECT_SUMMARY.md** for technical details
4. Check GitHub issues if using GitHub

---

**🎉 Project Complete and Ready for Deployment!**

Start the app with: `streamlit run app.py`

Good luck with your presentation! 🚀
