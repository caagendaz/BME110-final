# BioQuery NoLocal - Architecture & Flow Diagrams

## System Architecture Overview

```mermaid
graph TD
    UI["🌐 Streamlit Web App<br/>http://localhost:8501"]
    
    subgraph Frontend["Frontend Layer"]
        NLP["NLP Handler<br/>(nlp_handler.py)<br/>Google Gemini API<br/>gemini-2.5-flash"]
        ManualTool["Manual Tool Selection"]
        GenomeQuery["Genome Browser Query"]
    end
    
    subgraph Backend["Backend Processing"]
        Wrapper["EMBOSS Wrapper<br/>(emboss_wrapper.py)<br/>258+ tools + APIs"]
        GeneResolver["Gene Symbol Resolver<br/>Ensembl REST API"]
        TranscriptFetcher["Transcript Fetcher<br/>Ensembl API"]
        ToolExecutor["EMBOSS Tool Executor<br/>Subprocess runner"]
        BEDTools["BEDTools Executor<br/>Genomic overlaps"]
        ProteinAnalyzer["Protein Analyzer<br/>pepstats, iep"]
    end
    
    subgraph ExternalAPIs["External Databases & Tools"]
        Ensembl["Ensembl REST API<br/>Gene/Transcript/CDS"]
        UCSC["UCSC Genome Browser<br/>DAS API & BLAT"]
        NCBI["NCBI BLAST<br/>BioPython Remote"]
        GTEx["GTEx Portal<br/>Expression Links"]
    end
    
    UI -->|Natural Language Query| NLP
    UI -->|Manual Input| ManualTool
    UI -->|Genome Coordinates| GenomeQuery
    
    NLP -->|Parse with Gemini| Wrapper
    ManualTool -->|Tool + Params| Wrapper
    GenomeQuery -->|Coordinates| Wrapper
    
    Wrapper -->|Gene Symbol?| GeneResolver
    GeneResolver -->|Query Gene Info| Ensembl
    Ensembl -->|Transcripts & IDs| TranscriptFetcher
    TranscriptFetcher -->|Fetch CDS/cDNA| Ensembl
    
    Wrapper -->|DNA/RNA Tool| ToolExecutor
    Wrapper -->|Protein Analysis| ProteinAnalyzer
    Wrapper -->|BED Intersect| BEDTools
    ToolExecutor -->|Execute EMBOSS| Wrapper
    ProteinAnalyzer -->|pepstats/iep| Wrapper
    
    Wrapper -->|BLAST Query| NCBI
    Wrapper -->|Genomic Region| UCSC
    Wrapper -->|Expression Data| GTEx
    
    NCBI -->|Alignment Results| Wrapper
    UCSC -->|Sequence/BLAT| Wrapper
    GTEx -->|Portal Link| Wrapper
    
    Wrapper -->|Results| UI
    
    style UI fill:#1e88e5,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style NLP fill:#ff6f00,stroke:#e65100,color:#ffffff,font-weight:bold
    style ManualTool fill:#7b1fa2,stroke:#4a148c,color:#ffffff,font-weight:bold
    style GenomeQuery fill:#0288d1,stroke:#01579b,color:#ffffff,font-weight:bold
    style Wrapper fill:#388e3c,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style GeneResolver fill:#0097a7,stroke:#006064,color:#ffffff,font-weight:bold
    style TranscriptFetcher fill:#00838f,stroke:#004d40,color:#ffffff,font-weight:bold
    style ToolExecutor fill:#00796b,stroke:#004d40,color:#ffffff,font-weight:bold
    style BEDTools fill:#5e35b1,stroke:#4527a0,color:#ffffff,font-weight:bold
    style ProteinAnalyzer fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style Ensembl fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style UCSC fill:#f57c00,stroke:#e65100,color:#ffffff,font-weight:bold
    style NCBI fill:#fbc02d,stroke:#f57f17,color:#000000,font-weight:bold
    style GTEx fill:#7cb342,stroke:#558b2f,color:#ffffff,font-weight:bold
```

## Natural Language Query Flow

```mermaid
graph LR
    A["User Query<br/>e.g., 'Calculate molecular weight of MKTAYIAK'"] -->|Send to Gemini API| B["NLP Handler<br/>gemini-2.5-flash"]
    
    B -->|Parse with System Prompt| C{Query Type?}
    
    C -->|Gene Query| D["→ gene_query tool<br/>gene_name: ALKBH1"]
    C -->|EMBOSS Tool| E["→ translate tool<br/>gene_name: ALKBH1"]
    C -->|Protein Analysis| F["→ pepstats/iep<br/>sequence: MKTAYIAK"]
    C -->|Genome Region| G["→ genome_query<br/>chrom, start, end"]
    C -->|BLAST/BLAT| H["→ blast/blat<br/>sequence: ATGC..."]
    C -->|Multi-Step| I["→ steps array<br/>[gene_query, translate]"]
    
    D -->|Return JSON| J["Tool + Parameters"]
    E -->|Return JSON| J
    F -->|Return JSON| J
    G -->|Return JSON| J
    H -->|Return JSON| J
    I -->|Return JSON| J
    
    J -->|Pass to EMBOSS Wrapper| K["Execute Operation(s)"]
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#ff6f00,stroke:#e65100,color:#ffffff,font-weight:bold
    style C fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style D fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style E fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style F fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style G fill:#0288d1,stroke:#01579b,color:#ffffff,font-weight:bold
    style H fill:#5e35b1,stroke:#4527a0,color:#ffffff,font-weight:bold
    style I fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style J fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style K fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
```

## Gene-Based Tool Execution Pipeline

```mermaid
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
    UI --> TAB4["Tab 4:<br/>Documentation"]
    
    TAB1 --> NLP["Enter natural language<br/>Supports multi-question mode<br/>e.g., '1. GC of ALKBH1<br/>2. Translate TP53'"]
    NLP --> NLPRUN["NLP Handler (Gemini)<br/>processes & routes to tools"]
    NLPRUN --> NLP_OUT["Display formatted results<br/>Expandable sections per question"]
    
    TAB2 --> MANUAL["Select tool from dropdown<br/>258+ EMBOSS tools<br/>+ pepstats, iep, cusp, blast"]
    MANUAL --> MANUAL_RUN["Execute selected tool<br/>with gene name or sequence"]
    MANUAL_RUN --> MANUAL_OUT["Show output with download"]
    
    TAB3 --> GENOME["Enter genome coordinates<br/>hg38, chr1, 1000-2000<br/>Or single position"]
    GENOME --> GENOME_RUN["Query UCSC DAS API<br/>100bp window if single pos"]
    GENOME_RUN --> GENOME_OUT["Display sequence in FASTA"]
    
    TAB4 --> DOCS["View documentation<br/>Tool descriptions<br/>Usage examples"]
    
    style UI fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style TAB1 fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style TAB2 fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
    style TAB3 fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style TAB4 fill:#5e35b1,stroke:#3f51b5,color:#ffffff,font-weight:bold
```

## Data Flow: Gene Query to Result

```mermaid
sequenceDiagram
    participant User
    participant Streamlit as Streamlit UI
    participant NLP as NLP Handler<br/>(Gemini API)
    participant Wrapper as EMBOSS Wrapper
    participant Ensembl as Ensembl REST API
    participant EMBOSS as EMBOSS Tools

    User->>Streamlit: "Give protein length of CARS1 transcript 5"
    Streamlit->>NLP: parse_user_query()
    NLP->>NLP: Format system prompt + query
    NLP->>NLP: Call Gemini gemini-2.5-flash
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
    Streamlit->>User: Display formatted result with download button
```

## Key Features at a Glance

```mermaid
mindmap
  root((BioQuery NoLocal))
    Natural Language Interface
      Google Gemini API gemini-2.5-flash
      Automatic tool selection
      Gene symbol recognition
      Transcript variant support
      Multi-question batch processing
    EMBOSS Integration
      258+ bioinformatics tools
      Dynamic tool discovery
      Generic fallback for any tool
      Real-time execution
      Protein analysis pepstats iep cusp
    Gene-Based Operations
      Ensembl REST API integration
      Transcript resolution
      CDS/cDNA sequence fetching
      Variant-specific analysis
      Multi-step workflows
    Genomic Data Access
      UCSC Genome Browser DAS API
      BLAT near-exact search
      Region queries with auto-windowing
      No data storage needed
      Live API streaming
    Database Searches
      NCBI BLAST via BioPython
      Auto-detect protein vs DNA
      BEDTools genomic overlaps
      GTEx tissue expression links
    Streamlit Interface
      4 integrated tabs
      Manual & NLP modes
      Multi-question support
      Download results
    Advanced Capabilities
      Protein translation
      Molecular weight and pI
      Reverse complement
      Restriction sites
      GC content analysis
      Codon usage statistics
  
  %%{init: { 'theme': 'base', 'primaryColor':'#42a5f5', 'primaryTextColor':'#fff', 'primaryBorderColor':'#1e88e5', 'secondBkgColor':'#66bb6a', 'tertiaryColor':'#ef5350', 'tertiaryTextColor':'#fff', 'textPlacement': 'center', 'mindmapBkg':'transparent', 'nodeBkg':'transparent'} }%%
```

## File Structure

```mermaid
graph TB
    subgraph Code["Source Code"]
        APP["app.py<br/>Streamlit web interface<br/>Multi-question UI"]
        NLP["nlp_handler.py<br/>Google Gemini integration<br/>Multi-question parser"]
        EMBOSS["emboss_wrapper.py<br/>EMBOSS + BEDTools + APIs<br/>1665 lines"]
        APP -->|imports| NLP
        APP -->|imports| EMBOSS
    end
    
    subgraph Config["Configuration & Docs"]
        REQ["requirements.txt<br/>Python 3.12 deps"]
        README["README.md<br/>Complete guide"]
        COVERAGE["ASSIGNMENT_COVERAGE.md<br/>BME110 coverage"]
        TESTS["TEST_QUERIES.md<br/>Test examples"]
    end
    
    subgraph Setup["Setup Scripts"]
        SETUP_SH["setup.sh<br/>Linux/macOS setup"]
        RUN_SH["run.sh<br/>Start script"]
    end
    
    NLP -->|HTTPS| Gemini[("Google Gemini<br/>gemini-2.5-flash")]
    EMBOSS -->|subprocess| Tools[("EMBOSS 6.6.0<br/>BEDTools 2.31.1")]
    EMBOSS -->|HTTPS| Ensembl[("Ensembl<br/>REST API")]
    EMBOSS -->|HTTPS| UCSC[("UCSC<br/>DAS & BLAT")]
    EMBOSS -->|HTTPS| NCBI[("NCBI<br/>BLAST")]
    EMBOSS -->|Links| GTEx[("GTEx<br/>Portal")]
    
    style APP fill:#42a5f5,stroke:#1976d2,color:#ffffff,font-weight:bold
    style NLP fill:#66bb6a,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style EMBOSS fill:#ef5350,stroke:#c62828,color:#ffffff,font-weight:bold
    style REQ fill:#ab47bc,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style README fill:#ab47bc,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style COVERAGE fill:#ab47bc,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style TESTS fill:#ab47bc,stroke:#6a1b9a,color:#ffffff,font-weight:bold
    style SETUP_SH fill:#ffa726,stroke:#f57c00,color:#ffffff,font-weight:bold
    style RUN_SH fill:#ffa726,stroke:#f57c00,color:#ffffff,font-weight:bold
    style Gemini fill:#ffd54f,stroke:#ffa000,color:#000000,font-weight:bold
    style Tools fill:#66bb6a,stroke:#2e7d32,color:#ffffff,font-weight:bold
    style Ensembl fill:#ef5350,stroke:#c62828,color:#ffffff,font-weight:bold
    style UCSC fill:#ff7043,stroke:#d84315,color:#ffffff,font-weight:bold
    style NCBI fill:#ffa726,stroke:#f57c00,color:#ffffff,font-weight:bold
    style GTEx fill:#9ccc65,stroke:#689f38,color:#ffffff,font-weight:bold
```

## Technology Stack

```mermaid
graph LR
    subgraph Frontend["Frontend"]
        Streamlit["Streamlit 1.51.0<br/>4 tabs, multi-question UI"]
    end
    
    subgraph Backend["Backend"]
        Python["Python 3.12.12<br/>~25% faster than 3.9"]
        Bio["BioPython 1.86<br/>BLAST integration"]
        Pandas["Pandas 2.3.3<br/>Data handling"]
        Requests["Requests 2.32.3<br/>API calls"]
    end
    
    subgraph Tools["Bioinformatics Tools"]
        EMBOSS["EMBOSS 6.6.0<br/>258+ tools"]
        BEDTools["BEDTools 2.31.1<br/>Genomic overlaps"]
    end
    
    subgraph LLM["Cloud NLP"]
        Gemini["Google Gemini API<br/>gemini-2.5-flash<br/>10 req/min free"]
    end
    
    subgraph APIs["External APIs"]
        Ensembl["Ensembl REST API<br/>Gene/Transcript/CDS"]
        UCSC["UCSC DAS & BLAT<br/>Genomic regions"]
        NCBI["NCBI BLAST<br/>Sequence similarity"]
        GTEx["GTEx Portal<br/>Expression links"]
    end
    
    subgraph Env["Environment"]
        Conda["conda bioquery312<br/>Python 3.12 env"]
    end
    
    Streamlit -.->|runs on| Python
    Python -->|imports| Bio
    Python -->|imports| Pandas
    Python -->|imports| Requests
    Python -->|subprocess| EMBOSS
    Python -->|subprocess| BEDTools
    Python -->|HTTPS API calls| Gemini
    Python -->|HTTPS| Ensembl
    Python -->|HTTPS| UCSC
    Python -->|HTTPS via BioPython| NCBI
    Python -->|generates links| GTEx
    Conda -.->|manages| Python
    
    style Streamlit fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style Python fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style Bio fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style Pandas fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style Requests fill:#5e35b1,stroke:#4527a0,color:#ffffff,font-weight:bold
    style EMBOSS fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style BEDTools fill:#7b1fa2,stroke:#4a148c,color:#ffffff,font-weight:bold
    style Gemini fill:#ffd54f,stroke:#ffa000,color:#000000,font-weight:bold
    style Ensembl fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style UCSC fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style NCBI fill:#e57100,stroke:#d84315,color:#ffffff,font-weight:bold
    style GTEx fill:#7cb342,stroke:#558b2f,color:#ffffff,font-weight:bold
    style Conda fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
```
