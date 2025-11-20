# BioQuery NoLocal Test Queries

**IMPORTANT**: Gemini free tier allows 10 requests/minute. Wait ~60 seconds between test sets!

## Test Set 1: Basic Sequence Analysis (3 questions)
```
1. What is the GC content of the sequence ATGCGATCGATCGATCG?
2. Translate the DNA sequence ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG to protein.
3. Find the reverse complement of ATGCGATCG.
```

**Expected Results:**
- Q1: GC content ~52.94%
- Q2: Protein sequence MAIVMGR*KGAR*
- Q3: CGATCGCAT

---

## Test Set 2: Protein Analysis (3 questions)
```
1. Calculate the molecular weight of the protein sequence MKTAYIAKQRQISFVK.
2. What is the isoelectric point of MKTAYIAKQRQISFVK?
3. Show protein statistics for MKTAYIAKQRQISFVK.
```

**Expected Results:**
- Q1: Molecular weight from pepstats
- Q2: pI value
- Q3: Detailed amino acid composition

---

## Test Set 3: Database Searches (2 questions - SLOW!)
```
1. Find the codon usage for gene TP53.
2. Search for similar proteins to MYPFIRTARMTV using BLAST.
```

**Expected Results:**
- Q1: Codon frequency table from cusp
- Q2: BLAST hits (auto-detects protein, uses blastp/nr)

**NOTE**: These queries take 20-60 seconds each due to external APIs!

---

## Test Set 4: Genomic Queries (3 questions)
```
1. What gene is located at chromosome 17 position 7676154?
2. Which tissue has the highest expression of the SOCS3 gene?
3. Find genes overlapping with chromosome 1 positions 100000-200000.
```

**Expected Results:**
- Q1: 100bp sequence around position (UCSC API)
- Q2: GTEx portal link with instructions
- Q3: UCSC gene info for chr1 region

---

## Test Set 5: Multi-Step Workflow (2 questions)
```
1. Translate ATGGAGCAGAAACTCATCTCTGAAGAGGATCTG and calculate the molecular weight of the resulting protein.
2. Get the sequence for gene BRCA1, then calculate its GC content.
```

**Expected Results:**
- Q1: Translation + molecular weight in one result
- Q2: Gene sequence retrieval + GC calculation

**NOTE**: Multi-step queries use more Gemini tokens!

---

## Test Set 6: Rapid Simple Queries (5 questions)
```
1. GC content of ATGC?
2. Reverse complement of AAAA?
3. Translate ATG?
4. What is the length of ATGCGATCG?
5. Reverse complement of GCGC?
```

**Expected Results:**
- Q1: 50%
- Q2: TTTT
- Q3: M (single amino acid)
- Q4: 9 bp
- Q5: GCGC (palindrome)

---

## Tips

### Rate Limiting
- **Gemini free tier**: 10 requests/minute
- Each numbered question = 1 request
- Wait 60+ seconds between large test sets
- If you hit the limit, wait the suggested time (usually ~60s)

### Multi-Question Format
The app detects these formats:
- `1.` `2.` `3.` (number + period + space)
- `Question 1:` `Question 2:` (word + number + colon)
- `Q1:` `Q2:` (short form)

### Performance
- Simple queries (GC, reverse, translate): <1 second
- Database queries (BLAST, gene info): 10-60 seconds
- Complex multi-step: Variable (depends on steps)

### Known Issues Fixed
1. ✅ BLAST now auto-detects protein sequences
2. ✅ Molecular weight queries now use `pepstats` (not `info`)
3. ✅ Genome queries handle single positions (creates 100bp window)
4. ❌ Rate limiting is a hard limit (wait between sets)

### Best Practice
Start with **Test Set 1** (3 simple questions) to verify everything works, then try more complex sets!
