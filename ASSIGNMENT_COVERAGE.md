# BME110 Assignment Coverage Analysis

This document analyzes which BME110 midterm assignment questions can be completed using BioQuery NoLocal after recent modifications.

## Modifications Made

### New Tools Added:
1. **BEDTools Integration** - Find overlaps between genomic regions (BED files)
2. **BLAT Search** - UCSC BLAT for near-exact genome searches
3. **UCSC Gene Info** - Programmatic access to UCSC Genome Browser data
4. **Enhanced EMBOSS Tools** - Explicit support for `cusp`, `pepstats`, `wordcount`

### Updated Components:
- `src/emboss_wrapper.py` - Added 3 new methods: `bedtools_intersect()`, `blat_search()`, `ucsc_gene_info()`
- Tool mappings updated to include new keywords
- `requirements.txt` - Added `requests` library
- `README.md` - Documented new capabilities

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
- **Capability**: Calculate molecular weight, pI, amino acid composition
- **Example**: "Get protein stats for SOCS3"

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

#### Question 1: UCSC Genome Browser Exploration
- **Status**: API available via `ucsc_gene_info()` but not full browser UI
- **Limitation**: Cannot replicate visual exploration of tracks, GTEx expression, etc.
- **Workaround**: Provides browser link for manual follow-up

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

### Overall Completion Rate: ~80-85%

**Fully Automated** (10 questions): 2, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15

**Partially Automated** (3 questions): 1, 3, 4
- Question 4 (BEDTools) now fully supported with BED file input

**Requires Manual Work** (2 questions): Visual browser exploration, Galaxy integration

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

1. **EMBOSS Coverage**: All required EMBOSS tools (cusp, pepstats, wordcount, transeq, geecee, getorf, needle, water) now explicitly supported
2. **BED File Analysis**: BEDTools integration enables SNP/gene overlap analysis
3. **Search Flexibility**: Both BLAST (sensitive) and BLAT (fast, exact) available
4. **Gene-Centric Workflow**: Can query genes by name and automatically fetch sequences
5. **Multi-Step Pipelines**: Chain operations together (e.g., "get SOCS3 then calculate codon usage")

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
# Install BEDTools
conda install -c bioconda bedtools

# Install requests library
pip install requests

# Verify installation
bedtools --version
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

### Protein Statistics:
```
"Get protein stats for SOCS3 transcript variant 1"
```

## Conclusion

With these modifications, BioQuery NoLocal can now handle **~80-85% of the BME110 assignment programmatically**. The main gaps are visual exploration tasks that inherently require manual browser interaction. For a typical bioinformatics workflow combining automated analysis with manual exploration, this tool provides excellent coverage of the computational components.

**Recommended Approach**: Use BioQuery NoLocal for all computational analysis (sequence processing, alignments, statistics, BED intersections), then supplement with manual UCSC Browser exploration for visual inspection and context.
