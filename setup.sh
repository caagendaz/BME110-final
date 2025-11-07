#!/bin/bash
# Setup script for BioQuery Local
# Run this script to set up the entire project

set -e  # Exit on error

echo "🧬 BioQuery Local - Setup Script"
echo "=================================="
echo ""

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Miniconda or Anaconda first."
    echo "   Download from: https://docs.conda.io/projects/miniconda/en/latest/"
    exit 1
fi

echo "✓ Conda found"

# Create environment
echo ""
echo "📦 Creating conda environment 'bioquery' with Python 3.9..."
conda create -n bioquery python=3.9 -y

# Activate environment
echo "✓ Environment created"
echo ""
echo "📦 Installing bioinformatics tools..."
conda install -n bioquery -c bioconda emboss biopython -y

echo "✓ EMBOSS and BioPython installed"
echo ""
echo "📦 Installing Streamlit and pandas..."
conda install -n bioquery -c conda-forge streamlit pandas -y

echo "✓ Streamlit and pandas installed"
echo ""
echo "📦 Installing Python packages..."
conda run -n bioquery pip install biopython ollama

echo "✓ Python packages installed"
echo ""
echo "=================================="
echo "✅ Setup Complete!"
echo "=================================="
echo ""
echo "Next steps:"
echo "1. Start Ollama server:"
echo "   Windows: Run 'ollama serve' in PowerShell"
echo "   macOS/Linux: Run 'ollama serve' in terminal"
echo ""
echo "2. Pull a model (in another terminal):"
echo "   ollama pull gemma3:4b"
echo "   (or phi3:mini for a smaller model)"
echo ""
echo "3. Run the app:"
echo "   conda activate bioquery"
echo "   streamlit run app.py"
echo ""
echo "4. Open browser to: http://localhost:8501"
echo ""
