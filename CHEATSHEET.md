# BioQuery Quick Reference Cheat Sheet 🧬

## 🚀 Getting Started in 30 Seconds

```bash
conda activate bioquery
streamlit run src/app.py
# Opens at http://localhost:8501
```

**Choose Your Mode:**
- ☁️ **Cloud Mode**: Best accuracy, needs API key (`GOOGLE_API_KEY`)
- 🏠 **Local Mode**: Complete privacy, needs Ollama (`ollama pull llama3.2`)

---

## 💬 Natural Language Queries

### DNA Analysis
```
"Translate ATGAAATTTCCC to protein"
"What's the reverse complement of GCTAGC?"
"Find ORFs in ATGCCCGGGTTTAAA with minimum size 50"
"Calculate GC content of ATGCGCGATCG"
"Find CpG islands in this sequence: ATGCGC..."
```

### Protein Analysis
```
"Calculate molecular weight of MKTAYIAK"
"What's the isoelectric point of MKTAYIAKQRQ?"
"Show protein statistics for MKTAYIAK"
"Find hydrophobic regions in MKTAYIAK"
"Predict antigenic sites in this protein: MKLA..."
```

### Sequence Comparison
```
"Align these two sequences: ATGCCC and ATGAAA"
"Use stretcher to align [seq1] and [seq2]"
"Create a dot plot of these sequences"
"Use water for local alignment of [seq1] and [seq2]"
```

### Gene Queries (Cloud Mode Only)
```
"Translate the ALKBH1 gene"
"What's the GC content of TP53?"
"Find the isoelectric point of BRCA1"
"Get exon information for SOCS3"
```

### Multi-Step Workflows
```
"Find ALKBH1 gene info then BLAST it"
"Translate this DNA then calculate molecular weight"
"Get TP53 sequence then find restriction sites"
```

---

## 🔧 Quick Tool Reference

### DNA Tools
| Task | Natural Language | Tool Name |
|------|-----------------|-----------|
| DNA → Protein | "translate ATGC..." | `transeq` |
| Reverse complement | "reverse ATGC..." | `revseq` |
| Find ORFs | "find orfs in..." | `getorf` |
| GC content | "gc content of..." | `geecee` |
| CpG islands | "find cpg islands..." | `newcpgseek` |
| Restriction sites | "find restriction sites..." | `restrict` |
| Codon usage | "codon usage of..." | `cusp` |

### Protein Tools
| Task | Natural Language | Tool Name |
|------|-----------------|-----------|
| Molecular weight | "molecular weight of..." | `pepstats` |
| Isoelectric point | "isoelectric point of..." | `iep` |
| Protein stats | "protein statistics..." | `pepstats` |
| Hydrophobicity | "hydrophobic regions..." | `pepnet` |
| Antigenic sites | "antigenic sites in..." | `antigenic` |

### Alignment Tools
| Task | Natural Language | Tool Name |
|------|-----------------|-----------|
| Global alignment | "align [seq1] and [seq2]" | `needle` |
| Fast global | "use stretcher to align..." | `stretcher` |
| Local alignment | "use water to align..." | `water` |
| All local matches | "use matcher to find..." | `matcher` |
| Dot plot | "dot plot of [seq1] [seq2]" | `dotmatcher` |

### Search Tools (Cloud Mode)
| Task | Natural Language | Tool Name |
|------|-----------------|-----------|
| DNA BLAST | "blast this sequence..." | `blastn` |
| Protein BLAST | "search for similar proteins..." | `blastp` |
| Gene info | "find gene info for TP53" | `gene_query` |
| Expression data | "where is SOCS3 expressed?" | `gtex` |

---

## 📁 File Upload Tips

### Single File Upload
```
1. Upload FASTA file
2. Select tool from dropdown (e.g., "translate", "gc", "info")
3. Choose "Apply to all sequences"
4. Results for each sequence shown separately
```

### Two-File Comparison
```
1. Upload 2 FASTA files
2. Query: "align the sequences" or "use stretcher"
3. Automatically compares sequences across files
```

### Batch Processing Keywords
```
"translate all sequences"
"calculate gc content for all"
"find orfs in all sequences"
```

---

## 🎯 Pro Tips

### Sequence Formats Accepted
```
✅ Plain text:     ATGAAATTTCCC
✅ FASTA format:   >header\nATGAAATTTCCC
✅ Multi-line:     ATGAAA
                   TTTCCC
✅ With spaces:    ATG AAA TTT CCC (spaces removed)
```

### Speed Optimizations
```
✅ Sequences extracted BEFORE AI processing (10x faster)
✅ Tool resolution cached (2nd use instant)
✅ Use Local Mode for offline work
```

### Graphics Quality
```
✅ Dot plots use SVG (scalable, crisp)
✅ Download plots as .svg files
✅ View in browser or vector graphics editor
```

---

## 🔍 Examples by Use Case

### Molecular Biology Lab
```
"Find restriction sites in this plasmid: ATGC..."
"Translate this cDNA clone: ATGC..."
"Design primers for this gene: TP53"
"What's the melting temperature of ATGCGATCG?"
```

### Sequence Analysis
```
"Compare these two variants: [seq1] vs [seq2]"
"Find conserved regions in these sequences"
"Calculate protein properties for MKTAYIAK"
"Show all 6 reading frames of ATGCCC"
```

### Genomics Research
```
"Get expression data for ALKBH1 across tissues"
"Find genes near chr1:1000000-2000000"
"BLAST this sequence against human genome"
"Get BRCA1 transcript information"
```

### Bioinformatics Education
```
"What does this sequence encode? ATGAAATTT"
"Explain the reading frames of ATGCCC"
"Compare global vs local alignment of [seq1] [seq2]"
"Show codon usage bias in TP53"
```

---

## 🐛 Troubleshooting

### "Tool not found"
```
✅ Try: Use natural language ("find cpg islands") instead of tool name
✅ Try: Switch to Cloud Mode (Gemini more accurate)
✅ Check: Is EMBOSS installed? (`embossversion`)
```

### "Sequence too long"
```
✅ Use file upload instead of pasting
✅ For very long sequences (>10kb), use batch mode
✅ For genome regions, use UCSC query tab
```

### "API key error"
```
✅ Set: export GOOGLE_API_KEY='your_key'
✅ Or: Create .env file with GOOGLE_API_KEY=...
✅ Alternative: Use Local Mode (no API key needed)
```

### "Ollama not responding"
```
✅ Check: curl http://localhost:11434/api/tags
✅ Start: ollama serve
✅ Pull model: ollama pull llama3.2
```

---

## 📊 Feature Comparison

| Feature | Cloud Mode | Local Mode |
|---------|-----------|------------|
| **NLP Accuracy** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| **EMBOSS Tools** | 258+ | 258+ |
| **Biopython Tools** | 18 | 18 |
| **BLAST Search** | ✅ | ❌ |
| **Gene Queries** | ✅ | ❌ |
| **GTEx Expression** | ✅ | ❌ |
| **Internet Required** | Yes | No |
| **API Key** | Required | Not needed |
| **Privacy** | Cloud | 100% local |
| **Cost** | Free tier | Free unlimited |

---

## 🎓 Learning Path

### Beginner (5 minutes)
```
1. "Translate ATGAAATTT to protein"
2. "What's the reverse complement of ATGC?"
3. "Calculate GC content of ATGCGC"
```

### Intermediate (15 minutes)
```
1. Upload a FASTA file with multiple sequences
2. "Calculate molecular weight for all sequences"
3. "Find ORFs in all sequences with minimum size 100"
4. Try two-sequence comparison: "align seq1 and seq2"
```

### Advanced (30 minutes)
```
1. Multi-step: "Find ALKBH1 gene then BLAST it"
2. Custom tools: "use stretcher to align these sequences"
3. Export results: Download command log
4. Try Local Mode: Setup Ollama and compare
```

---

## 🔗 Quick Links

- **Documentation**: `README.md`
- **Ollama Setup**: `docs/OLLAMA_SETUP.md`
- **Architecture**: `docs/ARCHITECTURE.md`
- **Test Script**: `python test_ollama_resolution.py`
- **GitHub**: https://github.com/caagendaz/BME110-final

---

## 💡 Remember

1. **Natural language works best** - Describe what you want, not the tool name
2. **Sequences before AI** - Upload files for speed (10x faster)
3. **Caching helps** - Second query for same tool is instant
4. **Local mode for privacy** - Sensitive data stays on your machine
5. **Multi-step is smart** - Chain operations: "find X then do Y"

---

## 📝 Common Patterns

### Pattern: "From DNA to Protein Analysis"
```
"Translate ATGCCC then calculate molecular weight"
→ Step 1: transeq (DNA → protein)
→ Step 2: pepstats (protein stats)
```

### Pattern: "Gene to Function"
```
"Get TP53 sequence then find domains"
→ Step 1: gene_query (fetch sequence)
→ Step 2: Analysis tool
```

### Pattern: "Compare and Analyze"
```
"Align these sequences then find conserved regions"
→ Step 1: needle (alignment)
→ Step 2: Analysis of alignment
```

---

## 🎯 Quick Decision Tree

```
Need to analyze sequences?
├─ Have raw sequences? → Natural Language Query tab
├─ Have FASTA file? → File Upload tab
├─ Know exact tool? → Manual Tool Selection tab
├─ Need gene info? → Use gene symbols (Cloud Mode)
└─ Want privacy? → Switch to Local Mode
```

---

**Version 2.0** | Updated December 2024 | Supports 258+ EMBOSS tools + 18 Biopython tools
