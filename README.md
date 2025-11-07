# BioQuery Local 🧬

**Self-Contained Natural Language Bioinformatics Tool to do EMBOSS Analyses**

A Streamlit web application that uses natural language AI (Ollama) to convert user requests into EMBOSS bioinformatics analyses. Perfect for researchers who want to run bioinformatics tools without learning command-line syntax.

## Features

- 🤖 **Natural Language Interface**: Ask questions in plain English, get bioinformatics analysis
- 🧬 **EMBOSS Integration**: Access 13+ powerful EMBOSS tools
- 🎯 **AI-Powered**: Local LLM (Ollama) for intelligent tool selection
- 🌐 **Web Interface**: Beautiful Streamlit UI for easy interaction
- 📊 **Multiple Tools**: Translate, Reverse, ORF finding, Alignment, Restriction sites, and more
- 💾 **Results Export**: Download analysis results as text files

## Quick Start

### 1. Environment Setup

```bash
# Create conda environment
conda create -n bioquery python=3.9 -y
conda activate bioquery

# Install bioinformatics tools
conda install -c bioconda emboss biopython -y
conda install -c conda-forge streamlit pandas -y

# Install Python packages
pip install biopython ollama
```

### 2. Install Ollama

**On Windows:**
- Download from https://ollama.com/download/windows
- Run the installer
- Run in PowerShell: `$env:OLLAMA_HOST = "0.0.0.0:11434"; ollama serve`

**On macOS:**
```bash
brew install ollama
ollama serve
```

**On Linux:**
```bash
curl -fsSL https://ollama.ai/install.sh | sh
ollama serve
```

### 3. Pull a Model

In a new terminal:
```bash
ollama pull gemma3:4b
# Or: ollama pull phi3:mini
# Or: ollama pull llama3.2:3b
```

### 4. Run the App

```bash
conda activate bioquery
streamlit run app.py
```

The app will open at `http://localhost:8501`

## Usage

### Method 1: Natural Language Query (Recommended)
1. Go to the **"🤖 Natural Language Query"** tab
2. Type what you want to do, e.g.:
   - "Translate this DNA to protein: ATGAAATTTCCC"
   - "What's the reverse complement of GCTA?"
   - "Find all open reading frames in ATGAAATTTCCCGGGAAATTT"
3. Click **"🚀 Analyze"**
4. Results appear below, you can download them

### Method 2: Manual Tool Selection
1. Go to the **"🔧 Manual Tool Selection"** tab
2. Choose a tool from the dropdown
3. Enter your sequence(s)
4. Click the tool button
5. Results appear as formatted output

## Available Tools

| Tool | Description | Input |
|------|-------------|-------|
| **translate** | Translate DNA to protein | DNA sequence, reading frame (1-3) |
| **reverse** | Reverse complement DNA | DNA sequence |
| **orf** | Find open reading frames | DNA sequence, min size |
| **align** | Needleman-Wunsch alignment | Two sequences |
| **restriction** | Find restriction sites | DNA sequence, enzyme (optional) |
| **shuffle** | Randomize sequence | Sequence |
| **info** | Get sequence statistics | Sequence |
| **sixframe** | All 6 translation frames | DNA sequence |

## Project Structure

```
bme110/
├── app.py                 # Main Streamlit application
├── emboss_wrapper.py      # EMBOSS tool wrapper
├── nlp_handler.py         # Natural language processing handler
├── README.md              # This file
└── requirements.txt       # Python dependencies (optional)
```

## How It Works

1. **User Input** → Streamlit UI collects natural language query
2. **NLP Processing** → Ollama LLM interprets the query and selects appropriate tool
3. **Parameter Extraction** → LLM extracts sequence data and parameters from query
4. **EMBOSS Execution** → EMBOSSWrapper runs the selected EMBOSS tool
5. **Results Display** → Streamlit displays formatted results with download option

## System Requirements

- **OS**: Windows (WSL2), macOS, or Linux
- **Python**: 3.9+
- **RAM**: 4GB minimum (for model loading)
- **Disk**: 5GB (for models and EMBOSS tools)
- **Internet**: For initial model download

## Architecture

```
┌─────────────────────────────────────────┐
│     Streamlit Web Interface (app.py)    │
├─────────────────────────────────────────┤
│                                         │
│  ┌──────────────────┐  ┌──────────────┐ │
│  │  NLP Handler     │  │ EMBOSS       │ │
│  │  (nlp_handler)   │  │ Wrapper      │ │
│  │                  │  │ (emboss_     │ │
│  │  • Query Parse   │  │  wrapper)    │ │
│  │  • Tool Select   │  │              │ │
│  │  • Param Extract │  │ • translate  │ │
│  └────────┬─────────┘  │ • reverse    │ │
│           │            │ • orf        │ │
│           │            │ • align      │ │
│           └────┬──────►│ • ... etc    │ │
│                │       └────┬─────────┘ │
└────────────────┼────────────┼───────────┘
                 │            │
        ┌────────▼─────┐ ┌────▼──────┐
        │    Ollama    │ │  EMBOSS   │
        │ (AI Model)   │ │  Tools    │
        └──────────────┘ └───────────┘
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

**Author:** BME110 Student  
**Course:** BME110 - Bioinformatics Tool Project  
**Institution:** UC Santa Cruz
