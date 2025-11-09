# BioQuery Local - Architecture & Flow Diagrams

## System Architecture Overview

```mermaid
graph TD
    UI["🌐 Streamlit Web App<br/>http://localhost:8501"]
    
    subgraph Frontend["Frontend Layer"]
        NLP["NLP Handler<br/>(nlp_handler.py)<br/>Ollama gemma3:4b"]
        ManualTool["Manual Tool Selection"]
    end
    
    subgraph Backend["Backend Processing"]
        Wrapper["EMBOSS Wrapper<br/>(emboss_wrapper.py)<br/>258+ tools"]
        GeneResolver["Gene Symbol Resolver<br/>Ensembl API"]
        TranscriptFetcher["Transcript Fetcher<br/>Ensembl API"]
        ToolExecutor["EMBOSS Tool Executor<br/>Subprocess runner"]
    end
    
    subgraph ExternalAPIs["External Databases"]
        Ensembl["Ensembl REST API<br/>Gene/Transcript Info"]
        UCSC["UCSC Genome Browser<br/>Genomic Regions"]
        NCBI["NCBI Entrez<br/>Sequence Download"]
    end
    
    UI -->|Natural Language Query| NLP
    UI -->|Manual Input| ManualTool
    
    NLP -->|Parse Query| Wrapper
    ManualTool -->|Tool + Params| Wrapper
    
    Wrapper -->|Gene Symbol?| GeneResolver
    GeneResolver -->|Query Gene Info| Ensembl
    Ensembl -->|Transcripts & IDs| TranscriptFetcher
    TranscriptFetcher -->|Fetch Sequence| Ensembl
    
    Wrapper -->|Run EMBOSS Tool| ToolExecutor
    ToolExecutor -->|Execute| Wrapper
    
    Wrapper -->|Genomic Region?| UCSC
    Wrapper -->|Download Sequence?| NCBI
    
    Wrapper -->|Results| UI
    
    style UI fill:#1e88e5,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style NLP fill:#ff6f00,stroke:#e65100,color:#ffffff,font-weight:bold
    style ManualTool fill:#7b1fa2,stroke:#4a148c,color:#ffffff,font-weight:bold
    style Wrapper fill:#388e3c,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style GeneResolver fill:#0097a7,stroke:#006064,color:#ffffff,font-weight:bold
    style TranscriptFetcher fill:#00838f,stroke:#004d40,color:#ffffff,font-weight:bold
    style ToolExecutor fill:#00796b,stroke:#004d40,color:#ffffff,font-weight:bold
    style Ensembl fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style UCSC fill:#f57c00,stroke:#e65100,color:#ffffff,font-weight:bold
    style NCBI fill:#fbc02d,stroke:#f57f17,color:#000000,font-weight:bold
```

## Natural Language Query Flow

```mermaid
graph LR
    A["User Query<br/>e.g., 'Translate ALKBH1'"] -->|Send to Ollama| B["NLP Handler<br/>gemma3:4b Model"]
    
    B -->|Parse with System Prompt| C{Query Type?}
    
    C -->|Gene Query| D["→ gene_query tool<br/>gene_name: ALKBH1"]
    C -->|EMBOSS Tool| E["→ translate tool<br/>gene_name: ALKBH1"]
    C -->|Genome Region| F["→ genome_query<br/>chrom, start, end"]
    C -->|Raw Sequence| G["→ tool<br/>sequence: ATGC..."]
    
    D -->|Return JSON| H["Tool + Parameters"]
    E -->|Return JSON| H
    F -->|Return JSON| H
    G -->|Return JSON| H
    
    H -->|Pass to EMBOSS Wrapper| I["Execute Operation"]
    
    style A fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style B fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style C fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style D fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style E fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style F fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style G fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style H fill:#00695c,stroke:#004d40,color:#ffffff,font-weight:bold
    style I fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
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
    UI --> TAB4["Tab 4:<br/>Sequence Analysis"]
    UI --> TAB5["Tab 5:<br/>Documentation"]
    
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
    
    TAB5 --> DOCS["View documentation<br/>& API reference"]
    
    style UI fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style TAB1 fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style TAB2 fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
    style TAB3 fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style TAB4 fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style TAB5 fill:#5e35b1,stroke:#3f51b5,color:#ffffff,font-weight:bold
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
      Ollama gemma3:4b LLM
      Automatic tool selection
      Gene symbol recognition
      Transcript variant support
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
      5 integrated tabs
      Manual & NLP modes
      Batch processing
      Download results
    Advanced Capabilities
      Protein translation
      Reverse complement
      Restriction sites
      GC content analysis
      Isoelectric point calculation
  
  %%{init: { 'theme': 'base', 'primaryColor':'#42a5f5', 'primaryTextColor':'#fff', 'primaryBorderColor':'#1e88e5', 'secondBkgColor':'#66bb6a', 'tertiaryColor':'#ef5350', 'tertiaryTextColor':'#fff', 'textPlacement': 'center', 'mindmapBkg':'transparent', 'nodeBkg':'transparent'} }%%
```

## File Structure

```mermaid
graph TB
    subgraph Code["Source Code"]
        direction TB
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
    
    APP -->|HTTP requests| Ollama[("Ollama<br/>gemma3:4b")]
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

## Technology Stack

```mermaid
graph LR
    subgraph Frontend["Frontend"]
        Streamlit["Streamlit 1.48.1"]
    end
    
    subgraph Backend["Backend"]
        Python["Python 3.9"]
        Bio["BioPython 1.85"]
    end
    
    subgraph Tools["Bioinformatics Tools"]
        EMBOSS["EMBOSS 6.6.0<br/>258+ tools"]
    end
    
    subgraph LLM["LLM Engine"]
        Ollama["Ollama 0.6.0<br/>gemma3:4b Model"]
    end
    
    subgraph APIs["External APIs"]
        Ensembl["Ensembl REST API<br/>Gene/Transcript Data"]
        UCSC["UCSC DAS API<br/>Genomic Regions"]
        NCBI["NCBI Entrez<br/>Sequence Database"]
    end
    
    subgraph Env["Environment"]
        Conda["conda<br/>bioquery environment"]
    end
    
    Streamlit -.->|runs on| Python
    Python -->|imports| Bio
    Python -->|subprocess calls| EMBOSS
    Uses["uses"]
    Python -->|HTTP requests| Uses
    Uses -->|HTTP requests| Ollama
    Python -->|HTTP requests| Ensembl
    Python -->|HTTP requests| UCSC
    Python -->|HTTP requests| NCBI
    LLMLabel["interfaces<br/>with"]
    Ollama -.-> LLMLabel
    LLMLabel -.-> LLM
    Conda -.->|manages| Python
    
    style Streamlit fill:#f57f17,stroke:#e65100,color:#ffffff,font-weight:bold
    style Python fill:#1565c0,stroke:#0d47a1,color:#ffffff,font-weight:bold
    style Bio fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style EMBOSS fill:#2e7d32,stroke:#1b5e20,color:#ffffff,font-weight:bold
    style Ollama fill:#ffc400,stroke:#f57f17,color:#000000,font-weight:bold
    style Ensembl fill:#c62828,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style UCSC fill:#d32f2f,stroke:#b71c1c,color:#ffffff,font-weight:bold
    style NCBI fill:#e57100,stroke:#d84315,color:#ffffff,font-weight:bold
    style Conda fill:#6a1b9a,stroke:#4a148c,color:#ffffff,font-weight:bold
    style Uses fill:#ffffff,stroke:#cccccc,color:#000000,font-weight:bold
    style LLMLabel fill:#ffffff,stroke:#cccccc,color:#000000,font-weight:bold
```
