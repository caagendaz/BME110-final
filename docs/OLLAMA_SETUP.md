# Ollama Setup for Local AI Tool Resolution

## Overview

BioQuery now supports **fully local AI** using Ollama, enabling you to run tool resolution completely on your machine without internet connectivity or API keys. This is perfect for:

- 🔒 **Privacy-sensitive research** - sequences never leave your computer
- 🌐 **Offline work** - no internet required for AI queries
- 💰 **Cost control** - no API usage fees
- 🏠 **Data sovereignty** - complete control over your computational environment

## How It Works

When you ask "use stretcher to align these sequences" or "find CpG islands in this DNA", the system:

1. **Extracts your sequences** from the query
2. **Sends clean query to Ollama** running locally (e.g., "what EMBOSS tool is stretcher?")
3. **Ollama responds** with the tool name: "stretcher"
4. **System executes** the EMBOSS tool with your sequences
5. **Returns results** - all processing happened on your machine!

## Installation

### 1. Install Ollama

**macOS:**
```bash
brew install ollama
```

**Linux:**
```bash
curl -fsSL https://ollama.ai/install.sh | sh
```

**Windows:**
Download installer from https://ollama.ai

### 2. Pull a Model

Choose based on your needs:

```bash
# Recommended: Fast and accurate (2GB)
ollama pull llama3.2

# Alternative: Good for technical tasks (2-4GB)
ollama pull qwen2.5

# Alternative: Larger, more capable (4GB)
ollama pull mistral

# Advanced: Very capable (7GB)
ollama pull llama3.1:8b
```

### 3. Verify Installation

```bash
# Ollama should auto-start on port 11434
curl http://localhost:11434/api/tags

# Should return JSON with installed models
```

### 4. Configure BioQuery

In the app sidebar, select **"🏠 Local Mode (Ollama + EMBOSS only)"**

That's it! No API keys needed.

## Performance Comparison

| Feature | Gemini (Cloud) | Ollama (Local) |
|---------|---------------|----------------|
| **Accuracy** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Speed** | ~500ms | ~1-3s |
| **Privacy** | Cloud processing | 100% local |
| **Cost** | Free tier limits | Free unlimited |
| **Internet Required** | Yes | No (after setup) |
| **API Key** | Required | Not required |

## Model Recommendations

### llama3.2 (Recommended)
- **Size**: 2GB
- **Speed**: Fast (~1-2s per query)
- **Accuracy**: Excellent for tool resolution
- **Best for**: Most users, laptops, quick queries

### qwen2.5
- **Size**: 2-4GB (variants available)
- **Speed**: Fast
- **Accuracy**: Very good for technical/scientific queries
- **Best for**: Technical users, coding tasks

### mistral
- **Size**: 4GB
- **Speed**: Medium (~2-3s)
- **Accuracy**: Very high
- **Best for**: Complex queries, research

### llama3.1:8b
- **Size**: 7GB
- **Speed**: Slower (~3-5s)
- **Accuracy**: Highest
- **Best for**: Powerful workstations, complex reasoning

## Example Queries

All of these work with local Ollama (no internet needed for AI):

```
✅ "use stretcher to align these two sequences"
✅ "find CpG islands in this DNA: ATGCGCGATCG..."
✅ "use matcher to find local alignments"
✅ "calculate hydrophobicity with pepnet for this protein"
✅ "find antigenic sites in MKTAYIAK"
```

The AI resolves:
- `stretcher` → EMBOSS stretcher tool
- `CpG islands` → newcpgseek tool
- `matcher` → matcher tool (local alignment)
- `hydrophobicity with pepnet` → pepnet tool
- `antigenic sites` → antigenic tool

## Testing Tool Resolution

Run the test script to see Ollama in action:

```bash
cd /mnt/c/Users/cagns/OneDrive/Desktop/bme110
python test_ollama_resolution.py
```

This will show:
1. Which EMBOSS tools Ollama can identify
2. Resolution speed
3. Comparison with Gemini (if API key available)

## Troubleshooting

### "Ollama not available"
```bash
# Check if Ollama is running
curl http://localhost:11434/api/tags

# If not, start Ollama
ollama serve
```

### "Model not found"
```bash
# Pull the model
ollama pull llama3.2

# List installed models
ollama list
```

### Slow responses
```bash
# Use a smaller model
ollama pull llama3.2  # 2GB instead of 7GB

# Or optimize Ollama (macOS/Linux)
export OLLAMA_NUM_PARALLEL=4
```

### Tool resolution fails
- Ollama might not recognize very obscure EMBOSS tools
- Fallback: Use the explicit tool_map or switch to Cloud mode
- Report issues: We can add more training examples

## Privacy & Security

### What stays local?
✅ All sequences  
✅ All EMBOSS processing  
✅ All AI queries for tool resolution  
✅ All results and outputs  

### What requires internet?
❌ External databases (NCBI, Ensembl, UCSC) - only if you use them  
❌ Initial Ollama model download (~2-7GB one-time)  

### Data flow in Local Mode
```
Your Query
    ↓
[Extract sequences] ← Local Python
    ↓
[Send "what tool is X?"] → Ollama (localhost:11434) ← Local AI
    ↓
[Receive "tool Y"] ← Ollama response
    ↓
[Execute EMBOSS tool Y with sequences] ← Local EMBOSS
    ↓
Results (never left your machine!)
```

## Advanced Configuration

### Change default model
Edit `src/emboss_wrapper.py`:
```python
def __init__(self, ai_mode: str = "local", ai_model: str = None):
    # Change default from llama3.2 to your preferred model
    self.ai_model = "qwen2.5:latest" if ai_mode == "local" else "gemini-2.0-flash-exp"
```

### Use custom Ollama port
```bash
# Start Ollama on different port
OLLAMA_HOST=0.0.0.0:8080 ollama serve

# Update in emboss_wrapper.py
self.ollama_url = "http://localhost:8080"
```

### Optimize for low memory
```bash
# Use quantized model (smaller, faster, slightly less accurate)
ollama pull llama3.2:1b  # Only 1GB!
```

## Comparison: When to Use Each Mode

### Use Cloud Mode (Gemini) if:
- You need highest accuracy
- You're okay with API keys
- Internet is always available
- Working with complex multi-step queries
- Need external databases (BLAST, Ensembl, GTEx)

### Use Local Mode (Ollama) if:
- Privacy is critical (medical, proprietary data)
- No internet access required
- No API keys/costs desired
- Working with standard EMBOSS tools
- Running on laptop/offline
- Data sovereignty requirements

## Resources

- **Ollama**: https://ollama.ai
- **Ollama Models**: https://ollama.ai/library
- **Ollama GitHub**: https://github.com/ollama/ollama
- **Model Comparison**: https://ollama.ai/library (see benchmarks)

## Credits

Ollama support added December 2024. Implements the same tool resolution logic as Gemini but runs entirely on-device using open-source LLMs (Llama 3.2, Qwen 2.5, Mistral).
