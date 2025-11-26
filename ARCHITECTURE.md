# BioQuery Local - Architecture & Flow Diagrams

## System Architecture Overview

```mermaid
graph TD
    UI["🌐 Streamlit Web App<br/>http://localhost:8501"]
    
    subgraph Frontend["Frontend Layer - Dual Mode"]
        NLPRouter{NLP Mode<br/>Selection}
        CloudNLP["Cloud Mode<br/>Google Gemini 2.5 Flash"]
        LocalNLP["Local Mode<br/>Ollama gemma3:4b"]
        FallbackNLP["Ollama Fallback<br/>(on Gemini safety filter)"]
        ManualTool["Manual Tool Selection"]
    end
    
    subgraph Backend["Backend Processing - Enhanced"]
        Wrapper["EMBOSS Wrapper<br/>(emboss_wrapper.py)<br/>258+ tools"]
        GeneResolver["Gene Symbol Resolver<br/>Ensembl API"]
        TranscriptFetcher["Transcript Fetcher<br/>Ensembl API"]
        ToolExecutor["EMBOSS Tool Executor<br/>Subprocess runner"]
        NeighborGenes["Neighboring Genes Finder<br/>UCSC + Ensembl"]
        BLATMapper["BLAT Sequence Mapper<br/>Genome Location"]
    end
    
    subgraph ExternalAPIs["External Databases"]
        Ensembl["Ensembl REST API<br/>Gene/Transcript Info"]
        UCSC["UCSC Genome Browser<br/>Genomic Regions & Tracks"]
        NCBI["NCBI Entrez & BLAST<br/>Sequence Search"]
        GTEx["GTEx Portal<br/>Gene Expression"]
        PubMed["PubMed E-utilities<br/>Literature Search"]
    end
    
    UI -->|Natural Language Query| NLPRouter
    UI -->|Manual Input| ManualTool
    
    NLPRouter -->|Cloud Mode Selected| CloudNLP
    NLPRouter -->|Local Mode Selected| LocalNLP
    CloudNLP -->|Safety Filter Triggered| FallbackNLP
    
    CloudNLP -->|Parse Query| Wrapper
    LocalNLP -->|Parse Query| Wrapper
    FallbackNLP -->|Parse Query| Wrapper
    ManualTool -->|Tool + Params| Wrapper
    
    Wrapper -->|Gene Symbol?| GeneResolver
    Wrapper -->|Find Neighbors?| NeighborGenes
    Wrapper -->|Map Sequence?| BLATMapper
    
    GeneResolver -->|Query Gene Info| Ensembl
    Ensembl -->|Transcripts & IDs| TranscriptFetcher
    TranscriptFetcher -->|Fetch Sequence| Ensembl
    
    NeighborGenes -->|Get Gene Location| Ensembl
    NeighborGenes -->|Query Nearby Genes| UCSC
    
    BLATMapper -->|BLAT Search| UCSC
    BLATMapper -->|Get Annotations| UCSC
    
    Wrapper -->|Run EMBOSS Tool| ToolExecutor
    Wrapper -->|BLAST Search| NCBI
    Wrapper -->|Expression Data| GTEx
    Wrapper -->|Literature Search| PubMed
    
    ToolExecutor -->|Execute| Wrapper
    
    Wrapper -->|Results| UI
    
    style UI fill:#1e88e5,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style NLPRouter fill:#ff6f00,stroke:#e65100,color:#ffffff,font-weight:bold
    style CloudNLP fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style LocalNLP fill:#9c27b0,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style FallbackNLP fill:#ff9800,stroke:#e65100,color:#ffffff,font-weight:bold
    style ManualTool fill:#7b1fa2,stroke:#4a148c,color:#ffffff,font-weight:bold
    style Wrapper fill:#388e3c,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style GeneResolver fill:#0097a7,stroke:#006064,color:#ffffff,font-weight:bold
    style TranscriptFetcher fill:#00838f,stroke:#004d40,color:#ffffff,font-weight:bold
    style ToolExecutor fill:#00796b,stroke:#004d40,color:#ffffff,font-weight:bold
    style NeighborGenes fill:#1976d2,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style BLATMapper fill:#e91e63,stroke:#c2185b,color:#ffffff,font-weight:bold
    style Ensembl fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style UCSC fill:#f57c00,stroke:#e65100,color:#ffffff,font-weight:bold
    style NCBI fill:#fbc02d,stroke:#f57f17,color:#000000,font-weight:bold
    style GTEx fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style PubMed fill:#8bc34a,stroke:#558b2f,color:#ffffff,font-weight:bold
```

## Natural Language Query Flow - Cloud & Local Modes

```mermaid
graph TB
    A["User Query<br/>e.g., 'Characterize mystery sequence'"] -->|Send to NLP| B{Mode<br/>Selection}
    
    B -->|Cloud Mode| C["Google Gemini<br/>2.5 Flash"]
    B -->|Local Mode| D["Ollama<br/>gemma3:4b"]
    
    C -->|Parse with System Prompt| E{Safety<br/>Filter?}
    E -->|Blocked| F["⚠️ Fallback to Ollama<br/>(Cloud features remain active)"]
    E -->|Approved| G["Gemini Processing"]
    
    F -->|Process Query| H{Query Type?}
    G -->|Parse Query| H
    D -->|Parse Query| H
    
    H -->|Characterize/Mystery| I["Multi-Step Workflow:<br/>1. BLAT - map to genome<br/>2. ucsc_table - gene overlaps<br/>3. BLAST - RNA database search"]
    H -->|Gene Query| J["→ gene_query tool<br/>gene_name: ALKBH1"]
    H -->|EMBOSS Tool| K["→ translate tool<br/>gene_name: ALKBH1"]
    H -->|Neighboring Genes| L["→ find_neighboring_genes<br/>gene_name: CARS1"]
    H -->|Genome Region| M["→ genome_query<br/>chrom, start, end"]
    H -->|Raw Sequence| N["→ tool<br/>sequence: ATGC..."]
    H -->|BLAST Search| O["→ blast/blastn/blastp<br/>exclude_taxa, organism, etc."]
    
    I -->|Return Steps JSON| P["Tool + Parameters"]
    J -->|Return JSON| P
    K -->|Return JSON| P
    L -->|Return JSON| P
    M -->|Return JSON| P
    N -->|Return JSON| P
    O -->|Return JSON| P
    
    P -->|Pass to EMBOSS Wrapper| Q["Execute Operation(s)"]
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#ff6f00,stroke:#e65100,color:#ffffff,font-weight:bold
    style C fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style D fill:#9c27b0,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style E fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style F fill:#ff9800,stroke:#e65100,color:#ffffff,font-weight:bold
    style G fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style H fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style I fill:#e91e63,stroke:#c2185b,color:#ffffff,font-weight:bold
    style J fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style K fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style L fill:#1976d2,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style M fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style N fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style O fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style P fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style Q fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
```

## Mystery Sequence Characterization Workflow (NEW)

```mermaid
graph TD
    A["User: Paste assignment with<br/>FASTA sequences + questions"] 
    
    A -->|Smart FASTA Parser| B["Extract Context:<br/>- 'characterize'<br/>- 'mystery'<br/>- 'hg19'<br/>First sequence extracted"]
    
    B -->|Build Query| C["Cleaned Query:<br/>'Characterize this sequence using hg19'<br/>+ extracted sequence"]
    
    C -->|Send to NLP| D["NLP Creates Multi-Step:<br/>1. BLAT<br/>2. BLAST<br/>3. UCSC tracks (manual)"]
    
    D -->|Step 1: BLAT| E["Map to Genome<br/>Returns top hits with coordinates"]
    
    E -->|Extract Best Hit| F["chr6:135,418,563-135,418,717<br/>99.4% identity<br/>155 bp match"]
    
    F -->|Generate UCSC Links| G["UCSC Browser Links<br/>for manual inspection:<br/>- Gene annotations<br/>- Conservation tracks<br/>- ChromHMM states<br/>- POL2/POL3 ChIP-seq"]
    
    D -->|Step 2: BLAST| H["Search RefSeq RNA<br/>database=refseq_rna<br/>organism=Homo sapiens"]
    
    H -->|Top Matches| I["RNA Database Hits:<br/>- Accession IDs<br/>- Identity %<br/>- E-values<br/>- Gene names"]
    
    G -->|Manual Review| J["User clicks links to:<br/>a. Check gene overlaps<br/>b. View conservation (Multiz)<br/>c. Check chromHMM state<br/>d. Check POL2/3 signals"]
    
    I -->|Cross-reference| K["Match with:<br/>RNAcentral, Rfam, GtRNAdb"]
    
    J -->|Answer Questions| L["Complete Analysis:<br/>a. Gene name & function<br/>b. Neighboring genes<br/>c. Conservation across species<br/>d. Chromatin state<br/>e. Transcription evidence<br/>f. Database matches"]
    
    K -->|Combine Results| L
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style C fill:#9c27b0,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style D fill:#e91e63,stroke:#c2185b,color:#ffffff,font-weight:bold
    style E fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style F fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style G fill:#ff9800,stroke:#e65100,color:#ffffff,font-weight:bold
    style H fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style I fill:#8bc34a,stroke:#558b2f,color:#ffffff,font-weight:bold
    style J fill:#ff9800,stroke:#e65100,color:#ffffff,font-weight:bold
    style K fill:#8bc34a,stroke:#558b2f,color:#ffffff,font-weight:bold
    style L fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
```

## Neighboring Genes Discovery Flow (NEW)

```mermaid
sequenceDiagram
    participant User
    participant NLP as NLP Handler
    participant Wrapper as EMBOSS Wrapper
    participant Ensembl as Ensembl API
    participant UCSC as UCSC REST API

    User->>NLP: "What genes neighbor CARS1?"
    NLP->>NLP: Parse → find_neighboring_genes
    NLP->>Wrapper: tool='find_neighboring_genes'<br/>gene_name='CARS1', genome='hg38'
    
    Wrapper->>Ensembl: GET /lookup/symbol/homo_sapiens/CARS1
    Ensembl->>Wrapper: chr11:3,011,986-3,051,922<br/>26 transcripts
    
    Note over Wrapper: Expand region ±500kb
    Wrapper->>UCSC: GET /getData/track?genome=hg38<br/>&track=knownGene<br/>&chrom=chr11<br/>&start=2511986&end=3551922
    
    UCSC->>Wrapper: 47 gene annotations
    
    Note over Wrapper: Filter & process
    Wrapper->>Wrapper: 1. Filter: geneType=='protein_coding'<br/>2. Calculate distances from CARS1<br/>3. Find closest left & right<br/>4. Deduplicate by gene symbol
    
    Wrapper->>User: LEFT: NAP1L4 (~8,536 bp)<br/>RIGHT: OSBPL5 (~29,493 bp)<br/>Plus top 5 each side with distances
```
graph TD
    A["User: 'Translate CARS1<br/>transcript variant 5'"] 
    
    A -->|NLP Parsing| B["Identified:<br/>tool: translate<br/>gene_name: CARS1<br/>transcript_variant: 5"]
    
    B -->|Pass to run_tool| C["EMBOSS Wrapper<br/>run_tool()"]
    
    C -->|Query Gene Info| D["Ensembl API<br/>query_gene_info"]
    D -->|Returns 26 Transcripts| E["Find Transcript 5<br/>ENST00000380525"]
    
    E -->|Fetch CDS| F["get_transcript_sequence<br/>seq_type='cds'<br/>Result: 2496 bp"]
    
    F -->|Translate to Protein| G["translate_sequence<br/>frame=1"]
    G -->|EMBOSS transeq| H["Output: 832 aa<br/>Last 3: FQ*"]
    
    H -->|Format Result| I["Gene CARS1 from transcript variant 5:<br/>Protein Length: 832 aa<br/>Last 3 aa: FQ*"]
    
    I -->|Return to User| J["Display in Streamlit"]
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style C fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
    style D fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style E fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style F fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style G fill:#5e35b1,stroke:#3f51b5,color:#ffffff,font-weight:bold
    style H fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style I fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style J fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
```

## Neighboring Genes Discovery Flow (NEW)

```mermaid
sequenceDiagram
    participant User
    participant NLP as NLP Handler
    participant Wrapper as EMBOSS Wrapper
    participant Ensembl as Ensembl API
    participant UCSC as UCSC REST API

    User->>NLP: "What genes neighbor CARS1?"
    NLP->>NLP: Parse → find_neighboring_genes
    NLP->>Wrapper: tool='find_neighboring_genes'<br/>gene_name='CARS1', genome='hg38'
    
    Wrapper->>Ensembl: GET /lookup/symbol/homo_sapiens/CARS1
    Ensembl->>Wrapper: chr11:3,011,986-3,051,922<br/>26 transcripts
    
    Note over Wrapper: Expand region ±500kb
    Wrapper->>UCSC: GET /getData/track?genome=hg38<br/>&track=knownGene<br/>&chrom=chr11<br/>&start=2511986&end=3551922
    
    UCSC->>Wrapper: 47 gene annotations
    
    Note over Wrapper: Filter & process
    Wrapper->>Wrapper: 1. Filter: geneType=='protein_coding'<br/>2. Calculate distances from CARS1<br/>3. Find closest left & right<br/>4. Deduplicate by gene symbol
    
    Wrapper->>User: LEFT: NAP1L4 (~8,536 bp)<br/>RIGHT: OSBPL5 (~29,493 bp)<br/>Plus top 5 each side with distances
```

## Mystery Sequence Characterization Workflow (NEW)

```mermaid
graph TD
    A["User: Paste assignment with<br/>FASTA sequences + questions"] 
    
    A -->|Smart FASTA Parser| B["Extract Context:<br/>- 'characterize'<br/>- 'mystery'<br/>- 'hg19'<br/>First sequence extracted"]
    
    B -->|Build Query| C["Cleaned Query:<br/>'Characterize this sequence using hg19'<br/>+ extracted sequence"]
    
    C -->|Send to NLP| D["NLP Creates Multi-Step:<br/>1. BLAT<br/>2. BLAST<br/>3. UCSC tracks (manual)"]
    
    D -->|Step 1: BLAT| E["Map to Genome<br/>Returns top hits with coordinates"]
    
    E -->|Extract Best Hit| F["chr6:135,418,563-135,418,717<br/>99.4% identity<br/>155 bp match"]
    
    F -->|Generate UCSC Links| G["UCSC Browser Links<br/>for manual inspection:<br/>- Gene annotations<br/>- Conservation tracks<br/>- ChromHMM states<br/>- POL2/POL3 ChIP-seq"]
    
    D -->|Step 2: BLAST| H["Search RefSeq RNA<br/>database=refseq_rna<br/>organism=Homo sapiens"]
    
    H -->|Top Matches| I["RNA Database Hits:<br/>- Accession IDs<br/>- Identity %<br/>- E-values<br/>- Gene names"]
    
    G -->|Manual Review| J["User clicks links to:<br/>a. Check gene overlaps<br/>b. View conservation (Multiz)<br/>c. Check chromHMM state<br/>d. Check POL2/3 signals"]
    
    I -->|Cross-reference| K["Match with:<br/>RNAcentral, Rfam, GtRNAdb"]
    
    J -->|Answer Questions| L["Complete Analysis:<br/>a. Gene name & function<br/>b. Neighboring genes<br/>c. Conservation across species<br/>d. Chromatin state<br/>e. Transcription evidence<br/>f. Database matches"]
    
    K -->|Combine Results| L
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style C fill:#9c27b0,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style D fill:#e91e63,stroke:#c2185b,color:#ffffff,font-weight:bold
    style E fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style F fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style G fill:#ff9800,stroke:#e65100,color:#ffffff,font-weight:bold
    style H fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style I fill:#8bc34a,stroke:#558b2f,color:#ffffff,font-weight:bold
    style J fill:#ff9800,stroke:#e65100,color:#ffffff,font-weight:bold
    style K fill:#8bc34a,stroke:#558b2f,color:#ffffff,font-weight:bold
    style L fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
```

## EMBOSS Tool Resolution

```mermaid
graph LR
    A["User Tool Request<br/>e.g., 'gc' or 'iep'"]
    
    A -->|Check tool_map| B{Is in<br/>tool_map?}
    
    B -->|Yes| C["Route to Specific<br/>Implementation<br/>e.g., calculate_gc_content"]
    B -->|No| D["Generic EMBOSS<br/>Tool Handler<br/>_run_generic_emboss_tool"]
    
    C -->|Execute| E["Run EMBOSS Command"]
    D -->|Discover Tool| F["Check if tool exists<br/>in system PATH"]
    F -->|Execute| E
    
    E -->|Output Files| G["Parse Results"]
    G -->|Format| H["Return to User"]
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style C fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style D fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style E fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style G fill:#3f51b5,stroke:#283593,color:#ffffff,font-weight:bold
    style H fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
```

## Streamlit App Tab Flow

```mermaid
graph TB
    UI["Streamlit App<br/>http://localhost:8501"]
    
    UI --> TAB1["Tab 1:<br/>Natural Language Query"]
    UI --> TAB2["Tab 2:<br/>Manual Tool Selection"]
    UI --> TAB3["Tab 3:<br/>Genome Browser Query"]
    UI --> TAB4["Tab 4:<br/>Sequence Analysis"]
    UI --> TAB5["Tab 5:<br/>Command Log NEW"]
    UI --> TAB6["Tab 6:<br/>Documentation"]
    
    TAB1 --> NLP["Enter natural language<br/>e.g., 'Find GC of ALKBH1'"]
    NLP --> NLPRUN["NLP Handler processes<br/>& calls EMBOSS Wrapper"]
    NLPRUN --> NLP_OUT["Display formatted results"]
    
    TAB2 --> MANUAL["Select tool from dropdown<br/>Enter gene name or sequence"]
    MANUAL --> MANUAL_RUN["Execute selected tool"]
    MANUAL_RUN --> MANUAL_OUT["Show output"]
    
    TAB3 --> GENOME["Enter genome coordinates<br/>hg38, chr1, 1000-2000"]
    GENOME --> GENOME_RUN["Query UCSC API"]
    GENOME_RUN --> GENOME_OUT["Display region info"]
    
    TAB4 --> BATCH["Upload FASTA or enter<br/>multiple sequences"]
    BATCH --> BATCH_RUN["Process batch"]
    BATCH_RUN --> BATCH_OUT["Downloadable results"]
    
    TAB5 --> LOG["View command execution history"]
    LOG --> LOG_DETAIL["Timestamps, parameters,<br/>results, errors"]
    LOG_DETAIL --> LOG_DOWNLOAD["Download log or clear"]
    
    TAB6 --> DOCS["View documentation<br/>& API reference"]
    
    style UI fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style TAB1 fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style TAB2 fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
    style TAB3 fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style TAB4 fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style TAB5 fill:#00838f,stroke:#006064,color:#ffffff,font-weight:bold
    style TAB6 fill:#5e35b1,stroke:#3f51b5,color:#ffffff,font-weight:bold
```

## Command Logging System Flow (NEW)

```mermaid
sequenceDiagram
    participant User
    participant UI as Streamlit UI
    participant Wrapper as EMBOSS Wrapper
    participant Logger as Command Logger
    participant LogTab as Command Log Tab

    User->>UI: Execute any query/tool
    UI->>Wrapper: run_tool(tool_name, **params)
    
    Note over Wrapper: Wrapper function
    Wrapper->>Wrapper: _run_tool_internal(tool_name, **params)
    Wrapper->>Wrapper: Execute actual tool logic
    
    alt Success
        Wrapper->>Logger: _log_command(tool, params, result)
        Logger->>Logger: Store: timestamp, tool, params,<br/>result preview, success=True
    else Error
        Wrapper->>Logger: _log_command(tool, params, "", error=msg)
        Logger->>Logger: Store: timestamp, tool, params,<br/>error message, success=False
    end
    
    Wrapper->>UI: Return result
    UI->>User: Display result
    
    Note over User,LogTab: Later...
    User->>LogTab: Open Command Log tab
    LogTab->>Logger: get_command_log()
    Logger->>LogTab: Return all log entries
    LogTab->>User: Display formatted log with<br/>metrics, expandable entries
    
    User->>LogTab: Click Download Log
    LogTab->>Logger: get_formatted_log()
    Logger->>LogTab: Return text format
    LogTab->>User: Download .txt file
```

## Data Flow: Gene Query to Result

```mermaid
sequenceDiagram
    participant User
    participant Streamlit as Streamlit UI
    participant NLP as NLP Handler
    participant Wrapper as EMBOSS Wrapper
    participant Ensembl as Ensembl API
    participant EMBOSS as EMBOSS Tools

    User->>Streamlit: "Give protein length of CARS1 transcript 5"
    Streamlit->>NLP: parse_user_query()
    NLP->>NLP: Format system prompt + query
    NLP->>NLP: Call Ollama gemma3:4b
    NLP->>Streamlit: return {tool: 'translate', parameters: {gene_name: 'CARS1', transcript_variant: 'transcript variant 5'}}
    
    Streamlit->>Wrapper: run_tool('translate', gene_name='CARS1', transcript_variant='transcript variant 5')
    Wrapper->>Wrapper: Detect transcript_variant parameter
    Wrapper->>Ensembl: query_gene_info('CARS1')
    Ensembl->>Wrapper: [26 transcripts with IDs]
    Wrapper->>Wrapper: Parse "Transcript 5: ENST00000380525"
    
    Wrapper->>Ensembl: get_transcript_sequence('ENST00000380525', seq_type='cds')
    Ensembl->>Wrapper: return 2496 bp CDS sequence
    
    Wrapper->>EMBOSS: transeq -sequence input.fasta -outseq output.fasta
    EMBOSS->>Wrapper: 832 aa protein sequence
    
    Wrapper->>Streamlit: "Gene CARS1 from transcript variant 5:\nProtein Length: 832 aa\nLast 3: FQ*"
    Streamlit->>User: Display formatted result
```

## Key Features at a Glance

```mermaid
mindmap
  root((BioQuery Local))
    Natural Language Interface
      Google Gemini 2.5 Flash
      Automatic tool selection
      Gene symbol recognition
      Transcript variant support
      Academic context wrapper NEW
    EMBOSS Integration
      258+ bioinformatics tools
      Dynamic tool discovery
      Generic fallback for any tool
      Real-time execution
    Gene-Based Operations
      Ensembl API integration
      Transcript resolution
      CDS/cDNA sequence fetching
      Variant-specific analysis
    Genomic Data Access
      UCSC Genome Browser API
      Region queries
      No data storage needed
      Live API streaming
    Streamlit Interface
      6 integrated tabs NEW
      Manual & NLP modes
      Batch processing
      Download results
      Command logging NEW
    Advanced Capabilities
      Protein translation
      Reverse complement
      Restriction sites
      GC content analysis
      Isoelectric point calculation
      Execution tracking NEW
  
  %%{init: { 'theme': 'base', 'primaryColor':'#42a5f5', 'primaryTextColor':'#fff', 'primaryBorderColor':'#1e88e5', 'secondBkgColor':'#66bb6a', 'tertiaryColor':'#ef5350', 'tertiaryTextColor':'#fff', 'textPlacement': 'center', 'mindmapBkg':'transparent', 'nodeBkg':'transparent'} }%%
```

## Safety Filter Handling & Ollama Fallback (NEW)

```mermaid
graph TB
    A["User Query:<br/>'BLAST FBL mRNA excluding primates'"]
    
    A -->|Send to Cloud Mode| B["NLP Handler<br/>(Cloud Mode Active)"]
    
    B -->|Try Gemini First| C["Google Gemini 2.5 Flash"]
    
    C -->|Process Query| D{Safety<br/>Filter<br/>Triggered?}
    
    D -->|BLOCKED| E["⚠️ Safety filter detected:<br/>- Genetic sequences<br/>- Scientific terminology<br/>- Domain/taxonomy terms"]
    
    E -->|Auto-Fallback| F["Switch to Ollama<br/>(gemma3:4b)<br/>Cloud features remain active"]
    
    F -->|Reprocess Query| G["Ollama Processing<br/>No content restrictions"]
    
    D -->|APPROVED| H["Gemini Processing"]
    
    G -->|Parse Query| I["Extract Parameters:<br/>tool: blastn<br/>exclude_taxa: primates<br/>database: nt"]
    
    H -->|Parse Query| I
    
    I -->|Execute| J["BLAST Search<br/>with filters applied"]
    
    J -->|Results sorted by bit score| K["Display Results"]
    
    Note over F: Cloud features still work:<br/>- gene_query<br/>- genome_query<br/>- BLAST<br/>- GTEx<br/>- PubMed
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style C fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style D fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style E fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style F fill:#9c27b0,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style G fill:#9c27b0,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style H fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style I fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style J fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style K fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
```

```mermaid
graph LR
    A["User Query:<br/>'BLAST FBL in Archaea domain'"]
    
    A -->|Send to NLP| B["NLP Handler"]
    B -->|Wrap in context| C["Academic/Scientific<br/>Context Wrapper"]
    
    C -->|Generate prompt| D["ACADEMIC BIOINFORMATICS EXERCISE<br/>Context: university-level molecular biology<br/>Scientific Objective: taxonomic classification<br/>Technical Request: [user query]"]
    
    D -->|Send to Gemini| E{Safety<br/>Filter}
    
    E -->|Blocked| F["Error Message:<br/>Switch to Local (Ollama) mode<br/>- No content restrictions<br/>- Better for scientific queries"]
    E -->|Approved| G["Gemini Processing"]
    
    G -->|Parse query| H["Return:<br/>tool: blastn<br/>organism: Archaea"]
    
    H -->|Execute| I["BLAST Search"]
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style C fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style D fill:#00838f,stroke:#006064,color:#ffffff,font-weight:bold
    style E fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style F fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style G fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
    style H fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style I fill:#00838f,stroke:#006064,color:#ffffff,font-weight:bold
```

## File Structure

```mermaid
graph TB
    subgraph Code["Source Code"]
        APP["app.py<br/>Streamlit web interface"]
        NLP["nlp_handler.py<br/>Ollama NLP integration"]
        EMBOSS["emboss_wrapper.py<br/>EMBOSS tool wrapper"]
        APP -->|imports| NLP
        APP -->|imports| EMBOSS
    end
    
    subgraph Config["Configuration"]
        REQ["requirements.txt<br/>Python dependencies"]
        README["README.md<br/>Project documentation"]
        GETTING["GETTING_STARTED.md<br/>Setup instructions"]
    end
    
    subgraph Setup["Setup Scripts"]
        SETUP_SH["setup.sh<br/>Linux/macOS setup"]
        SETUP_PS["setup_windows.ps1<br/>Windows setup"]
    end
    
    NLP -->|HTTP requests| Ollama[("Ollama<br/>gemma3:4b")]
    EMBOSS -->|subprocess| Tools[("EMBOSS<br/>Tools")]
    EMBOSS -->|HTTP requests| Ensembl[("Ensembl<br/>API")]
    
    style APP fill:#42a5f5,stroke:#1976d2,color:#ffffff,font-weight:bold
    style NLP fill:#66bb6a,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style EMBOSS fill:#ef5350,stroke:#c62828,color:#ffffff,font-weight:bold
    style REQ fill:#ab47bc,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style README fill:#ab47bc,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style GETTING fill:#ab47bc,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style SETUP_SH fill:#ffa726,stroke:#f57c00,color:#ffffff,font-weight:bold
    style SETUP_PS fill:#ffa726,stroke:#f57c00,color:#ffffff,font-weight:bold
    style Ollama fill:#42a5f5,stroke:#1976d2,color:#ffffff,font-weight:bold
    style Tools fill:#66bb6a,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style Ensembl fill:#ef5350,stroke:#c62828,color:#ffffff,font-weight:bold
```

## Technology Stack - Enhanced v3.0

```mermaid
graph LR
    subgraph Frontend["Frontend"]
        Streamlit["Streamlit 1.48.1<br/>6 Tabs + File Upload"]
    end
    
    subgraph Backend["Backend"]
        Python["Python 3.9-3.12"]
        Bio["BioPython 1.85"]
    end
    
    subgraph Tools["Bioinformatics Tools"]
        EMBOSS["EMBOSS 6.6.0<br/>258+ tools"]
    end
    
    subgraph LLM["Dual LLM Engine"]
        Gemini["Google Gemini<br/>2.5 Flash Model<br/>(Cloud Mode)"]
        Ollama["Ollama gemma3:4b<br/>(Local Mode +<br/>Auto Fallback)"]
    end
    
    subgraph APIs["External APIs"]
        Ensembl["Ensembl REST API<br/>Gene/Transcript Data"]
        UCSC["UCSC Genome Browser<br/>Genomic Regions & Tracks<br/>BLAT Mapping"]
        NCBI["NCBI BLAST<br/>Sequence Homology<br/>megablast support"]
        GTEx["GTEx Portal<br/>Gene Expression"]
        PubMed["PubMed E-utilities<br/>Literature Search"]
    end
    
    subgraph Env["Environment"]
        Conda["conda<br/>bioquery environment"]
    end
    
    subgraph NewFeatures["New Features v3.0"]
        Logger["Command Logger<br/>Execution Tracking"]
        SafetyFallback["Ollama Auto-Fallback<br/>on Gemini safety filter"]
        NeighborFunc["Neighboring Genes<br/>Discovery Function"]
        MysteryWorkflow["Mystery Sequence<br/>Multi-Step Analysis"]
        BLATEnhanced["Enhanced BLAT Parser<br/>Direct UCSC Links"]
        FASTAParser["Smart FASTA Parser<br/>Multi-sequence support"]
        BlastFixes["BLAST Enhancements<br/>megablast + bit score sort"]
    end
    
    Streamlit -.->|runs on| Python
    Python -->|imports| Bio
    Python -->|subprocess calls| EMBOSS
    Python -->|HTTP requests| Gemini
    Python -->|HTTP requests| Ollama
    Gemini -.->|on block| SafetyFallback
    SafetyFallback -.->|routes to| Ollama
    Python -->|HTTP requests| Ensembl
    Python -->|HTTP requests| UCSC
    Python -->|HTTP requests| NCBI
    Python -->|HTTP requests| GTEx
    Python -->|HTTP requests| PubMed
    Conda -.->|manages| Python
    Python -->|integrates| Logger
    Python -->|includes| NeighborFunc
    Python -->|implements| MysteryWorkflow
    Python -->|uses| BLATEnhanced
    Python -->|uses| FASTAParser
    Python -->|applies| BlastFixes
    
    style Streamlit fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style Python fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style Bio fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style EMBOSS fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style Gemini fill:#4caf50,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style Ollama fill:#9c27b0,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style Ensembl fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style UCSC fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style NCBI fill:#e57100,stroke:#d84315,color:#ffffff,font-weight:bold
    style GTEx fill:#00838f,stroke:#006064,color:#ffffff,font-weight:bold
    style PubMed fill:#8bc34a,stroke:#558b2f,color:#ffffff,font-weight:bold
    style Conda fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
    style Logger fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
    style SafetyFallback fill:#ff9800,stroke:#e65100,color:#ffffff,font-weight:bold
    style NeighborFunc fill:#1976d2,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style MysteryWorkflow fill:#e91e63,stroke:#c2185b,color:#ffffff,font-weight:bold
    style BLATEnhanced fill:#e91e63,stroke:#c2185b,color:#ffffff,font-weight:bold
    style FASTAParser fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style BlastFixes fill:#00acc1,stroke:#00838f,color:#ffffff,font-weight:bold
```
