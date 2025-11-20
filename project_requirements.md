
1

Automatic Zoom
Shell
"BioQuery Local" - Self-Contained Natural Language 
Bioinformatics Tool to do EMBOSS Analyses 
 
This is a suggested AI-driven alternative to Final Group Project.  Feel free to do this as 
described, add, change, or improve.  Email me if you decide to select this alternative to 
the Web Analysis Tool project.  You can work on your own, or in groups of two for this 
project.  
 
There will be all the same “Demo” requirements and timeline as defined in main project, but also requires (1) 
puting your finished work in GitHub,  and (2) anyone else in the class (including me!) should be able to 
download and run what you present. 
Initial Setup (Pre-Week 1) 
# Create conda environment 
conda create -n bioquery python=3.10 
conda activate bioquery 
 
# Install bioinformatics tools 
conda install -c bioconda emboss biopython 
conda install -c conda-forge streamlit pandas 
 
# Install Ollama (macOS) 
brew install ollama 
 
# Pull a lightweight model (Phi-3 is small and fast) 
ollama pull phi3:mini  # 2GB model, good for bioinformatics 
# Alternative: ollama pull llama3.2:3b  # Slightly larger but more capable 
 
# Install Python Ollama client 
pip install ollama 
Week 1: EMBOSS Integration & Basic Framework 
Create a wrapper for EMBOSS tools with natural language mapping: 
Python
# emboss_wrapper.py 
import subprocess 
import tempfile 
import os 
from pathlib import Path 
 
class EMBOSSWrapper: 
    """Wrapper for common EMBOSS tools""" 
     
    def __init__(self): 
        # Map natural language to EMBOSS commands 
        self.tool_map = { 
            'translate': 'transeq', 
            'reverse': 'revseq', 
            'orf': 'getorf', 
            'align': 'needle', 
            'pattern': 'fuzznuc', 
            'restriction': 'restrict', 
            'shuffle': 'shuffleseq', 
            'info': 'infoseq', 
            'sixframe': 'sixpack' 
        } 
         
        # Verify EMBOSS installation 
        self.check_emboss() 
     
    def check_emboss(self): 
        """Verify EMBOSS tools are available""" 
        try: 
            result = subprocess.run(['embossversion'],  
                                  capture_output=True, text=True) 
            print(f"EMBOSS found: {result.stdout.strip()}") 
            return True 
        except FileNotFoundError: 
            print("EMBOSS not found. Please install with: conda install -c 
bioconda emboss") 
            return False 
     
