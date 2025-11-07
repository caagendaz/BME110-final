# Setup script for BioQuery Local (Windows PowerShell)
# Run this script with: powershell -ExecutionPolicy Bypass -File setup_windows.ps1

Write-Host "🧬 BioQuery Local - Windows Setup Script" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Check if conda is installed
try {
    conda --version | Out-Null
    Write-Host "✓ Conda found" -ForegroundColor Green
} catch {
    Write-Host "❌ Conda not found. Please install Miniconda or Anaconda first." -ForegroundColor Red
    Write-Host "   Download from: https://docs.conda.io/projects/miniconda/en/latest/" -ForegroundColor Yellow
    exit 1
}

# Create environment
Write-Host ""
Write-Host "📦 Creating conda environment 'bioquery' with Python 3.9..." -ForegroundColor Yellow
conda create -n bioquery python=3.9 -y

Write-Host "✓ Environment created" -ForegroundColor Green
Write-Host ""
Write-Host "📦 Installing bioinformatics tools..." -ForegroundColor Yellow
conda install -n bioquery -c bioconda emboss biopython -y

Write-Host "✓ EMBOSS and BioPython installed" -ForegroundColor Green
Write-Host ""
Write-Host "📦 Installing Streamlit and pandas..." -ForegroundColor Yellow
conda install -n bioquery -c conda-forge streamlit pandas -y

Write-Host "✓ Streamlit and pandas installed" -ForegroundColor Green
Write-Host ""
Write-Host "📦 Installing Python packages..." -ForegroundColor Yellow
conda run -n bioquery pip install biopython ollama

Write-Host "✓ Python packages installed" -ForegroundColor Green
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "✅ Setup Complete!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Cyan
Write-Host "1. Start Ollama server (in a new PowerShell):" -ForegroundColor White
Write-Host '   $env:OLLAMA_HOST = "0.0.0.0:11434"' -ForegroundColor Gray
Write-Host "   ollama serve" -ForegroundColor Gray
Write-Host ""
Write-Host "2. Pull a model (in another PowerShell terminal):" -ForegroundColor White
Write-Host "   ollama pull gemma3:4b" -ForegroundColor Gray
Write-Host "   (or: ollama pull phi3:mini for smaller model)" -ForegroundColor Gray
Write-Host ""
Write-Host "3. Run the app (in WSL or PowerShell):" -ForegroundColor White
Write-Host "   conda activate bioquery" -ForegroundColor Gray
Write-Host "   streamlit run app.py" -ForegroundColor Gray
Write-Host ""
Write-Host "4. Open browser to: http://localhost:8501" -ForegroundColor Green
Write-Host ""
