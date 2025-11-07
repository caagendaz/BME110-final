# GitHub Setup Guide for BioQuery Local

## Step 1: Initialize Git Repository

From the project directory (`/mnt/c/Users/cagns/OneDrive/Desktop/bme110`):

```bash
git init
git config user.name "Your Name"
git config user.email "your.email@ucsc.edu"
```

## Step 2: Add Files to Git

```bash
git add .
git status  # Review what will be committed
```

## Step 3: Make Initial Commit

```bash
git commit -m "Initial commit: BioQuery Local project

- EMBOSS wrapper for bioinformatics tools
- NLP handler using Ollama for natural language queries
- Streamlit web interface
- Setup scripts for easy installation
- Complete documentation and README"
```

## Step 4: Create GitHub Repository

1. Go to https://github.com/new
2. Create a new repository named `bioquery-local` (or similar)
3. **Do NOT** initialize with README, .gitignore, or license (we already have these)
4. Click **Create repository**

## Step 5: Connect Local to Remote

After creating the GitHub repository, copy the commands it shows and run them:

```bash
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/bioquery-local.git
git push -u origin main
```

Replace `YOUR_USERNAME` with your GitHub username.

## Step 6: Verify on GitHub

- Go to https://github.com/YOUR_USERNAME/bioquery-local
- Check that all files are uploaded
- Verify the README displays correctly

## Step 7: Add Collaborators (Optional)

If working in a group:
1. Go to Settings → Collaborators
2. Add team members by username

## Updating Your Repository

After making changes:

```bash
git add .
git commit -m "Brief description of changes"
git push
```

## File Structure for GitHub

Your repository should have:

```
bioquery-local/
├── README.md           # Project documentation
├── requirements.txt    # Python dependencies
├── app.py              # Main Streamlit app
├── emboss_wrapper.py   # EMBOSS tool wrapper
├── nlp_handler.py      # NLP processing
├── setup.sh            # Linux/macOS setup script
├── setup_windows.ps1   # Windows setup script
├── .gitignore          # Git ignore file
└── project_requirements.md  # Original requirements
```

## Important Notes

1. **Do NOT commit:**
   - `__pycache__/` directories (handled by .gitignore)
   - `.streamlit/` directory
   - Large model files
   - Personal configuration files

2. **Large file handling:**
   - EMBOSS and Ollama are installed separately, not in the repo
   - Only Python scripts go in the repository

3. **Cloning your project:**
   ```bash
   git clone https://github.com/YOUR_USERNAME/bioquery-local.git
   cd bioquery-local
   bash setup.sh  # or setup_windows.ps1 on Windows
   ```

## First Time GitHub Setup

If this is your first time using Git:

```bash
# Set global config
git config --global user.name "Your Name"
git config --global user.email "your.email@ucsc.edu"

# Generate SSH key (recommended)
ssh-keygen -t ed25519 -C "your.email@ucsc.edu"

# Add SSH key to GitHub
# 1. Copy public key: cat ~/.ssh/id_ed25519.pub
# 2. Go to GitHub Settings → SSH and GPG keys
# 3. Click "New SSH key" and paste
```

## Troubleshooting

**"Permission denied" when pushing?**
- Make sure your SSH key is added to GitHub
- Or use HTTPS instead of SSH

**"fatal: not a git repository"?**
- Make sure you're in the project directory
- Run `git init` if needed

**Want to change the remote URL?**
```bash
git remote set-url origin https://github.com/YOUR_USERNAME/bioquery-local.git
```

---

Good luck with your project! 🚀
