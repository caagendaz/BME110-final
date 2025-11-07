# BioQuery Local Configuration Guide

## Ollama Configuration

### Host Configuration

If Ollama is running on a different machine or port, edit `nlp_handler.py` line 12:

```python
# Default (Windows from WSL2):
def __init__(self, ollama_host: str = "http://192.168.128.1:11434", model: str = "gemma3:4b"):

# Change to your setup:
def __init__(self, ollama_host: str = "http://your.ip:11434", model: str = "your_model"):
```

### Supported Models

#### Lightweight (Recommended for laptops)
- **phi3:mini** (~2GB)
  - Fast, good for bioinformatics
  - `ollama pull phi3:mini`

- **gemma3:4b** (~3.3GB) ✓ Currently configured
  - Balanced speed and quality
  - `ollama pull gemma3:4b`

#### Medium
- **llama3.2:3b** (~2GB)
  - Good reasoning
  - `ollama pull llama3.2:3b`

- **mistral:7b** (~4GB)
  - More capable
  - `ollama pull mistral:7b`

#### Larger (if you have >8GB RAM)
- **llama3:8b** (~4.7GB)
  - Much better understanding
  - `ollama pull llama3:8b`

### Changing Models

1. Pull the new model:
   ```bash
   ollama pull mistral:7b
   ```

2. Update `nlp_handler.py`:
   ```python
   def __init__(self, ollama_host: str = "http://192.168.128.1:11434", model: str = "mistral:7b"):
   ```

3. Restart the app

## Streamlit Configuration

### Custom Theme

Create `.streamlit/config.toml`:

```toml
[theme]
primaryColor = "#2E86AB"
backgroundColor = "#F0F8FF"
secondaryBackgroundColor = "#E8E8FF"
textColor = "#262730"
font = "sans serif"

[client]
showErrorDetails = true

[logger]
level = "info"
```

### Performance Tuning

For large sequences, edit `app.py`:

```python
# Increase timeouts (line ~200 in emboss_wrapper.py)
result = subprocess.run(
    cmd,
    capture_output=True,
    text=True,
    timeout=60  # Increase from 30 for very large sequences
)

# Cache settings in app.py
@st.cache_resource
def initialize_tools():
    # Can adjust cache behavior here
```

## EMBOSS Configuration

### Custom Paths

If EMBOSS is in a non-standard location, update `emboss_wrapper.py`:

```python
def check_emboss(self) -> bool:
    # Add to PATH or specify full path
    import os
    emboss_path = "/path/to/emboss/bin"
    os.environ['PATH'] = emboss_path + os.pathsep + os.environ.get('PATH', '')
    
    try:
        result = subprocess.run(
            ['embossversion'],
            capture_output=True,
            text=True,
            timeout=5
        )
```

### Temporary File Location

Change temp directory in `emboss_wrapper.py`:

```python
# Default behavior (uses system temp)
with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
    
# Custom location:
import tempfile
tempdir = '/tmp/bioquery'  # or your preferred location
with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, dir=tempdir) as f:
```

## Network Configuration

### WSL2 Networking

If `192.168.128.1` doesn't work, find the correct gateway:

```bash
# In WSL2 terminal:
ip route | grep default | awk '{print $3}'
```

This shows your Windows IP. Update `nlp_handler.py` with this IP.

### Firewall Issues

**Windows Firewall:**
```powershell
# Add exception for Ollama (as Administrator)
New-NetFirewallRule -DisplayName "Ollama" -Direction Inbound -Action Allow -Protocol TCP -LocalPort 11434
```

**Linux Firewall (UFW):**
```bash
sudo ufw allow 8501  # Streamlit
sudo ufw allow 11434  # Ollama (if on same machine)
```

## Environment Variables

### Ollama Environment

```bash
# Set before running 'ollama serve'
set OLLAMA_HOST=0.0.0.0:11434
set OLLAMA_NUM_PARALLEL=4  # Increase for faster processing
set OLLAMA_MAX_LOADED_MODELS=2
set OLLAMA_DEBUG=true  # For troubleshooting
```

### Python Environment

```bash
# Set before running streamlit
set PYTHONUNBUFFERED=1  # Real-time output
set PYTHONPATH=/path/to/project
```

## Performance Tips

1. **Model Loading**
   - First run after restart may take 1-2 minutes
   - Subsequent runs are much faster
   - Keep Ollama running in background

2. **Large Sequences**
   - Break into smaller chunks if possible
   - Increase timeout values
   - Use appropriate tool (not all tools handle large sequences well)

3. **Memory**
   - Close other applications
   - Use smaller model (phi3:mini) if running out of memory
   - Monitor with `top` or Task Manager

4. **Caching**
   - Streamlit caches the tools with `@st.cache_resource`
   - Clear cache if you update code: Delete `.streamlit/` folder

## Logging & Debugging

### Enable Debug Mode

1. Edit `nlp_handler.py`:
   ```python
   # Add after imports
   import logging
   logging.basicConfig(level=logging.DEBUG)
   logger = logging.getLogger(__name__)
   ```

2. Run Streamlit with debug:
   ```bash
   streamlit run app.py --logger.level=debug
   ```

### Check Ollama Logs

```bash
# Windows: Check Task Manager for Ollama process
# Linux/macOS: 
tail -f ~/.ollama/logs/server.log

# Or run in foreground to see output:
ollama serve
```

## Database Integration (Advanced)

For future enhancement, add SQLite support:

```python
# In app.py
import sqlite3

def save_analysis(tool, sequence, result):
    conn = sqlite3.connect('bioquery.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS analyses
                 (id INTEGER PRIMARY KEY, tool TEXT, sequence TEXT, result TEXT, timestamp DATETIME DEFAULT CURRENT_TIMESTAMP)''')
    c.execute("INSERT INTO analyses (tool, sequence, result) VALUES (?, ?, ?)", 
              (tool, sequence, result))
    conn.commit()
    conn.close()
```

## Security Hardening

For production deployment:

1. **Add Authentication**
   - Use Streamlit's built-in auth
   - Or add external auth (OAuth, etc.)

2. **Rate Limiting**
   - Add request limits to prevent abuse
   - Limit file sizes

3. **Input Validation**
   - Already implemented in `nlp_handler.py`
   - Can add more strict validation

4. **HTTPS**
   - Use Streamlit Cloud or reverse proxy
   - Enable SSL certificates

## Troubleshooting Configurations

### Problem: "Connection refused"
**Solution**: Update host IP in `nlp_handler.py`
```bash
# Find correct IP:
ip route | grep default
```

### Problem: Model takes too long to load
**Solution**: Use smaller model
```bash
ollama pull phi3:mini
# Update model in nlp_handler.py
```

### Problem: Out of memory errors
**Solution**: Reduce model or sequence size
- Use phi3:mini (2GB vs 3.3GB)
- Break large sequences into chunks

### Problem: Streamlit crashes on large file
**Solution**: Increase upload limit in `.streamlit/config.toml`
```toml
[client]
maxUploadSize = 200
```

---

For more help, check README.md or PROJECT_SUMMARY.md
