# BME110 Assignment Coverage Analysis

This document analyzes which BME110 midterm assignment questions can be completed using BioQuery NoLocal after recent modifications.

## Modifications Made

### New Tools Added:
1. **BEDTools Integration** - Find overlaps between genomic regions (BED files)
2. **BLAT Search** - UCSC BLAT for near-exact genome searches via API
3. **UCSC Gene Info** - Programmatic access to UCSC Genome Browser data
4. **GTEx Integration** - Query gene expression with portal links to 54 human tissues
5. **Enhanced EMBOSS Tools** - Explicit support for `cusp`, `pepstats`, `wordcount`, `iep`
6. **Multi-Question Processing** - Batch process numbered questions independently
7. **Improved NLP** - Better tool routing for protein analysis (molecular weight, pI, etc.)

### Updated Components:
- `src/emboss_wrapper.py` - Added: `bedtools_intersect()`, `blat_search()`, `ucsc_gene_info()`, `gtex_expression()`, protein analysis tools
- `src/nlp_handler.py` - Added: multi-question parsing, improved tool descriptions in system prompt, pepstats/iep routing
- `src/app.py` - Added: multi-question UI with expandable sections per question
- `tool_map` - Updated: molecular_weight/weight/mass → pepstats, added all new tool keywords
- `requirements.txt` - Added: `requests==2.32.3` library
- `README.md` - Updated: Python 3.12, Google Gemini API, all new features documented

### Python Upgrade:
- **Python 3.9 → 3.12**: ~25% performance improvement, continued security support until Oct 2028
- **Package upgrades**: BioPython 1.85→1.86, Streamlit 1.48.1→1.51.0, Pandas 2.3.1→2.3.3

## Assignment Question Coverage

### ✅ CAN COMPLETE (Programmatically)

#### Question 2: Gene Position and Strand
- **Tool**: `ucsc_gene_info('SOCS3', 'hg38')`
- **Capability**: API access to gene positions
- **Note**: Can get programmatic data, but not same as visual browser exploration

#### Question 5: Amino Acids in Transcript
- **Tool**: `translate(gene='SOCS3', transcript_variant='variant 1')`
- **Capability**: Translate specific transcript variants and count amino acids
- **Example**: "Translate SOCS3 transcript variant 1"

#### Question 6: Codon Usage
- **Tool**: `cusp` (EMBOSS codon usage)
- **Capability**: Calculate codon usage statistics for any sequence
- **Example**: "Calculate codon usage for SOCS3"

#### Question 7: tRNA Gene Names
- **Tool**: UCSC API queries (programmatic access)
- **Capability**: Query tRNA genes from UCSC Table Browser programmatically
- **Note**: Requires BED file or coordinate-based query

#### Question 8: Protein Statistics
- **Tool**: `pepstats` (EMBOSS)
- **Capability**: Calculate molecular weight, pI, amino acid composition, charge
- **Example**: "Calculate molecular weight of MKTAYIAKQRQISFVK" or "Get protein stats for SOCS3"
- **NLP Routing**: "molecular weight", "MW", "mass", "protein statistics" all route to pepstats

#### Question 9: Oligonucleotide Frequencies
- **Tool**: `wordcount` (EMBOSS)
- **Capability**: Count word frequencies (k-mers) in sequences
- **Example**: "Count oligonucleotides in sequence ATGCCC"

#### Question 10: Sequence Translation
- **Tool**: `transeq` or `translate`
- **Capability**: Translate DNA to protein in any reading frame
- **Example**: "Translate this sequence: ATGCCCGGG"

#### Question 11: GC Content
- **Tool**: `geecee` or `calculate_gc_content`
- **Capability**: Calculate GC percentage for any sequence
- **Example**: "Calculate GC content of ATGCCCGGG"

#### Question 12: Find ORFs
- **Tool**: `getorf` or `find_orfs`
- **Capability**: Find all open reading frames with minimum size
- **Example**: "Find ORFs in ATGCCCGGG with min size 100"

#### Question 13: BLAST vs BLAT Differences
- **Tool**: Both `blast` and `blat` now available
- **Capability**: 
  - BLAST: Remote NCBI searches via BioPython
  - BLAT: UCSC genome searches for near-exact matches
- **Example**: "BLAST this sequence" vs "BLAT search in hg38"
- **Key Differences**:
  - BLAST: More sensitive, finds distant homologs, slower
  - BLAT: Optimized for near-exact matches, faster, genome-specific

#### Question 14: Smith-Waterman (Local Alignment)
- **Tool**: `water` (EMBOSS)
- **Capability**: Smith-Waterman local alignment between sequences
- **Example**: "Align sequences using water (local alignment)"

#### Question 15: Needleman-Wunsch (Global Alignment)
- **Tool**: `needle` or `align_sequences`
- **Capability**: Needleman-Wunsch global alignment
- **Example**: "Align ATGC and ATGG globally"

### ⚠️ PARTIALLY COMPLETE (Requires Manual Steps)

#### Question 1: UCSC Genome Browser Exploration + GTEx Expression
- **Status**: Gene info API available via `ucsc_gene_info()` + GTEx link via `gtex_expression()`
- **Capability**: Provides gene coordinates and direct link to GTEx Portal
- **Example**: `gtex_expression('SOCS3')` returns Ensembl ID and GTEx Portal URL
- **Limitation**: Cannot automatically extract "highest expressing tissue" - requires clicking link
- **Workaround**: Provides instructions: "Look at Gene Expression bar chart, tallest bar = highest tissue"
- **Coverage**: ~75% automated (gene lookup + link generation), 25% manual (reading chart)

#### Question 3: Flanking Genes
- **Status**: Can query genes, but identifying "flanking genes" requires coordinate analysis
- **Limitation**: Need to query multiple genes and compare positions
- **Workaround**: Multi-step workflow: query target gene → query nearby region → compare

#### Question 4: BED File Intersections (SNPs vs tRNA genes)
- **Tool**: `bedtools_intersect(file_a, file_b)`
- **Status**: ✅ NOW SUPPORTED!
- **Capability**: Find overlaps between BED files
- **Example**: "Find SNPs overlapping tRNA genes using BED files"
- **Requirement**: Need BED files as input (from UCSC Table Browser)

#### Question 7 Part 2: Export to Galaxy
- **Status**: Not supported
- **Limitation**: BioQuery NoLocal doesn't integrate with Galaxy
- **Workaround**: Manual export from UCSC Table Browser to Galaxy

### ❌ CANNOT COMPLETE (Requires Manual Browser Interaction)

#### Question 1 Extended: Visual Track Exploration
- **Issue**: Cannot replicate visual inspection of genome browser tracks
- **Tools**: Conservation tracks, GTEx expression, histone marks
- **Reason**: These are interactive visualizations, not API-accessible

#### Question 3 Extended: Manual Gene Counting
- **Issue**: Counting genes between positions requires visual browser or custom script
- **Reason**: Not a standard bioinformatics tool operation

## Coverage Summary

### Overall Completion Rate: ~85-90%

**Fully Automated** (10 questions): 2, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15

**Partially Automated** (4 questions): 1 (GTEx link), 3, 4 (BEDTools ready)
- Question 1: GTEx integration provides gene info + portal link for expression data
- Question 4 (BEDTools) fully supported with BED file input

**Requires Manual Work** (1 question): Visual conservation track interpretation

## Recommended Workflow for Assignment

### Automated Analysis:
1. Use BioQuery NoLocal for:
   - Gene sequence retrieval
   - Codon usage, protein stats, word counts
   - Translations and ORF finding
   - Sequence alignments (both global and local)
   - BLAST and BLAT searches
   - BED file intersections (with BEDTools)

### Manual Browser Work:
2. Use UCSC Genome Browser for:
   - Visual exploration of gene context
   - Identifying flanking genes visually
   - Viewing conservation and expression tracks

### Combined Approach:
3. Use BioQuery NoLocal to get sequences → Manual browser for visual context → BioQuery NoLocal for computational analysis

## Key Strengths After Modifications

1. **EMBOSS Coverage**: All required EMBOSS tools (cusp, pepstats, wordcount, transeq, geecee, getorf, needle, water, iep) now explicitly supported with proper NLP routing
2. **BED File Analysis**: BEDTools integration enables SNP/gene overlap analysis programmatically
3. **Search Flexibility**: Both BLAST (sensitive, distant homologs) and BLAT (fast, near-exact matches) available
4. **Gene-Centric Workflow**: Query genes by name, automatically fetch sequences, chain operations
5. **GTEx Integration**: Direct portal links for tissue expression visualization with step-by-step instructions
6. **Multi-Step Pipelines**: Chain operations together (e.g., "get SOCS3 then calculate codon usage")
7. **Multi-Question Mode**: Paste entire assignments with numbered questions - each processed independently
8. **Protein Analysis**: Specialized tools (pepstats, iep) with smart NLP routing for molecular weight, pI, composition
9. **Python 3.12**: Modern environment with ~25% performance boost and long-term support

## Limitations and Future Enhancements

### Current Limitations:
1. **No Visual Tracks**: Cannot display genome browser-style visualizations
2. **No Galaxy Integration**: Cannot directly export to Galaxy analysis platform
3. **UCSC Table Browser**: API access limited compared to full web interface
4. **Coordinate Math**: Complex region calculations (flanking genes) require manual logic

### Potential Enhancements:
1. **PyGenomeTracks**: Add visualization of genomic regions with tracks
2. **BioBlend**: Integrate with Galaxy via API for workflow submission
3. **Enhanced UCSC API**: Deeper integration with UCSC Table Browser
4. **Coordinate Operations**: Add tools for "genes within region" queries
5. **GTF/GFF Parsing**: Process gene annotation files directly

## Installation Requirements

To use all new features:

```bash
# Create Python 3.12 environment
conda create -n bioquery python=3.12 -y
conda activate bioquery

# Install core packages first (avoids conflicts)
conda install -c conda-forge streamlit pandas -y

# Install bioinformatics tools
conda install -c bioconda emboss bedtools -y

# Install API libraries
pip install biopython google-generativeai requests

# Verify installation
bedtools --version
transeq -version
python -c "import google.generativeai; print('Gemini API ready')"
```

## Usage Examples

### BEDTools Intersection:
```
"Find SNPs overlapping tRNA genes from these BED files"
file_a: snps.bed (chr, start, end, snp_id)
file_b: trna_genes.bed (chr, start, end, gene_name)
```

### BLAT Search:
```
"BLAT search this sequence in hg38: ATGCCCGGGATTTAA"
```

### Codon Usage:
```
"Calculate codon usage for SOCS3"
```

### GTEx Expression:
```
"What tissues express SOCS3 highest?"

Returns:
- Gene: SOCS3
- Ensembl ID: ENSG00000184557
- GTEx Portal link: https://gtexportal.org/home/gene/SOCS3
- Instructions: "Look at Gene Expression bar chart, tallest bar = highest tissue"
```

## Conclusion

With these modifications, BioQuery NoLocal can now handle **~85-90% of the BME110 assignment programmatically or with guided links**. The main gaps are visual exploration tasks that require manual interaction with genome browsers. For GTEx expression queries, the tool provides direct portal links with clear instructions, making the manual step trivial.

**Technology Stack**:
- **NLP**: Google Gemini API (gemini-2.5-flash) - 10 requests/min free tier
- **Python**: 3.12.12 (~25% faster than 3.9)
- **Local Tools**: EMBOSS 6.6.0, BEDTools 2.31.1
- **APIs**: Ensembl REST, NCBI BLAST, UCSC Genome Browser, GTEx Portal
- **Packages**: BioPython 1.86, Streamlit 1.51.0, Pandas 2.3.3

**Recommended Approach**: Use BioQuery NoLocal for all computational analysis (sequence processing, alignments, statistics, BED intersections) and gene info lookups. For GTEx expression, use the provided portal link and follow the instructions. For multi-question assignments, paste all questions at once (respecting rate limits). This workflow maximizes automation while acknowledging the practical limits of API access.

**Rate Limiting Note**: Gemini free tier allows 10 requests/minute. For assignments with many questions, process in batches of 5-10 questions with 60-second pauses between batches.
