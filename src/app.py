"""
BioQuery Local - Streamlit App
Interactive web interface for natural language bioinformatics analysis
"""

import streamlit as st
from emboss_wrapper import EMBOSSWrapper
from nlp_handler import NLPHandler
import traceback
import os
from datetime import datetime
from PIL import Image


def display_result(result):
    """Display result, handling both text and image outputs"""
    try:
        if isinstance(result, str) and result.startswith("IMAGE_FILE:"):
            # Extract image path
            image_path = result.replace("IMAGE_FILE:", "")
            
            try:
                # Read the image file
                with open(image_path, "rb") as file:
                    image_data = file.read()
                
                st.success("✓ Dot plot image generated successfully!")
                
                # Try to display the image
                try:
                    from PIL import Image
                    image = Image.open(image_path)
                    st.image(image, caption="Dot Plot", use_container_width=True)
                except:
                    # If PIL fails, just show the download button
                    st.info("Image generated. Use the download button below to view it.")
                
                # Provide download button
                st.download_button(
                    label="📥 Download Dot Plot Image",
                    data=image_data,
                    file_name="dotplot.png",
                    mime="image/png"
                )
                
                # Clean up temp file after reading
                try:
                    os.remove(image_path)
                except:
                    pass
                    
            except Exception as e:
                st.error(f"Error handling image: {str(e)}")
                st.text(f"Image path: {image_path}")
        else:
            # Regular text output - ensure it's a string
            if not isinstance(result, str):
                result = str(result)
            st.code(result, language="text")
    except RecursionError:
        st.error("Display error: maximum recursion depth exceeded")
        st.text(str(result)[:500])  # Show first 500 chars as fallback
    except Exception as e:
        st.error(f"Error displaying result: {str(e)}")
        st.text(str(result)[:500])


# Page configuration
st.set_page_config(
    page_title="BioQuery Local",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state for mode selection
if 'mode' not in st.session_state:
    st.session_state.mode = 'cloud'  # Default to cloud mode

# Mode selector in sidebar
with st.sidebar:
    st.title("⚙️ Configuration")
    
    mode = st.radio(
        "Select Mode:",
        options=['cloud', 'local'],
        format_func=lambda x: "☁️ Cloud Mode (Gemini + APIs)" if x == 'cloud' else "🏠 Local Mode (Ollama + EMBOSS only)",
        index=0 if st.session_state.mode == 'cloud' else 1,
        help="Cloud mode uses Google Gemini API and external databases. Local mode uses Ollama and only local EMBOSS tools."
    )
    
    # Update session state if changed
    if mode != st.session_state.mode:
        st.session_state.mode = mode
        st.rerun()
    
    st.markdown("---")
    
    if mode == 'cloud':
        st.info("**Cloud Mode Features:**\n- Google Gemini AI\n- NCBI BLAST\n- Ensembl API\n- UCSC Genome Browser\n- GTEx Expression\n- All EMBOSS tools")
    else:
        st.info("**Local Mode Features:**\n- Ollama (local LLM)\n- EMBOSS tools only\n- No internet required\n- Complete privacy")
    
    st.markdown("---")

# Custom CSS for better styling
st.markdown("""
<style>
    .main {
        padding: 2rem;
    }
    .stTabs [data-baseweb="tab-list"] button {
        font-size: 1.1rem;
    }
    .tool-box {
        border: 2px solid #2E86AB;
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
        background-color: #f0f8ff;
    }
    .error-box {
        border: 2px solid #d32f2f;
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
        background-color: #ffebee;
    }
    .success-box {
        border: 2px solid #388e3c;
        border-radius: 8px;
        padding: 1rem;
        margin: 1rem 0;
        background-color: #e8f5e9;
    }
</style>
""", unsafe_allow_html=True)


@st.cache_resource
def initialize_tools(mode='cloud'):
    """Initialize EMBOSS and NLP tools (cached for performance)
    
    Args:
        mode: 'cloud' for Gemini API or 'local' for Ollama
    """
    emboss = EMBOSSWrapper()
    nlp = NLPHandler(mode=mode)
    return emboss, nlp


def parse_fasta(fasta_content):
    """Parse FASTA content and return list of (header, sequence) tuples
    
    Args:
        fasta_content: String containing FASTA formatted sequences
    
    Returns:
        List of tuples: [(header1, sequence1), (header2, sequence2), ...]
    """
    sequences = []
    current_header = None
    current_seq = []
    
    # Handle different line endings and strip BOM if present
    fasta_content = fasta_content.replace('\r\n', '\n').replace('\r', '\n')
    if fasta_content.startswith('\ufeff'):
        fasta_content = fasta_content[1:]
    
    lines = fasta_content.split('\n')
    
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue
        
        # Check for FASTA header (supports both > and ; as some formats use semicolon)
        if line.startswith('>') or (line.startswith(';') and current_header is None):
            # Save previous sequence if exists
            if current_header is not None:
                seq_string = ''.join(current_seq)
                if seq_string:  # Only add if sequence is not empty
                    sequences.append((current_header, seq_string))
                else:
                    st.warning(f"Empty sequence for header: {current_header}")
            
            # Start new sequence
            current_header = line[1:].strip()  # Remove '>' or ';'
            current_seq = []
        elif line.startswith(';'):
            # Comment line in FASTA (some formats use semicolons for comments)
            continue
        else:
            # Sequence line - validate it contains valid characters
            if current_header is None:
                # Found sequence data before any header - might not be FASTA format
                st.warning(f"Line {line_num}: Found sequence data before FASTA header (line starts with: '{line[:20]}...')")
                continue
            
            # Remove spaces and numbers (some FASTA files have line numbers)
            clean_line = ''.join(c for c in line if c.isalpha() or c == '-' or c == '*')
            if clean_line:
                current_seq.append(clean_line)
    
    # Don't forget the last sequence
    if current_header is not None:
        seq_string = ''.join(current_seq)
        if seq_string:
            sequences.append((current_header, seq_string))
        else:
            st.warning(f"Empty sequence for header: {current_header}")
    
    # If no sequences found, provide helpful error
    if not sequences and fasta_content.strip():
        st.error("❌ No valid FASTA sequences found. FASTA files must:\n"
                 "- Start headers with '>' character\n"
                 "- Have at least one header line\n"
                 "- Contain sequence data after headers")
        st.text("First 200 characters of file:")
        st.code(fasta_content[:200])
    
    return sequences


def main():
    """Main Streamlit application"""
    
    # Initialize session state for sequence storage
    if 'sequences' not in st.session_state:
        st.session_state.sequences = {}
    
    # Get current mode from session state
    current_mode = st.session_state.mode
    
    # Header
    col1, col2 = st.columns([1, 4])
    with col1:
        st.title("🧬")
    with col2:
        if current_mode == 'cloud':
            st.title("BioQuery NoLocal")
            st.caption("☁️ Cloud-powered natural language bioinformatics analysis tool")
        else:
            st.title("BioQuery Local")
            st.caption("🏠 Privacy-first local bioinformatics analysis tool")
    
    st.markdown("---")
    
    # Initialize tools with selected mode
    with st.spinner("Initializing bioinformatics tools..."):
        try:
            emboss, nlp = initialize_tools(mode=current_mode)
            emboss_ready = emboss.check_emboss()
            nlp_ready = nlp.test_connection()
        except Exception as e:
            st.error(f"Failed to initialize tools: {str(e)}")
            return
    
    # Sidebar with status and info
    with st.sidebar:
        st.header("📊 System Status")
        
        # Status indicators
        col1, col2 = st.columns(2)
        with col1:
            if emboss_ready:
                st.success("✓ EMBOSS")
            else:
                st.error("✗ EMBOSS")
        with col2:
            if nlp_ready:
                if current_mode == 'cloud':
                    st.success("✓ Gemini")
                else:
                    st.success("✓ Ollama")
            else:
                if current_mode == 'cloud':
                    st.error("✗ Gemini")
                else:
                    st.error("✗ Ollama")
        
        st.markdown("---")
        
        # Available tools
        st.subheader("📊 Available Tools")
        if current_mode == 'local':
            st.caption("Local mode: EMBOSS tools only")
        else:
            st.caption("Cloud mode: EMBOSS + external APIs")
        
        tools = emboss.get_available_tools()
        for tool_name, description in tools.items():
            with st.expander(f"**{tool_name}**"):
                st.write(description)
        
        st.markdown("---")
        
        # About
        st.subheader("ℹ️ About")
        if current_mode == 'cloud':
            st.info(
                "BioQuery NoLocal uses:\n"
                "- **EMBOSS**: Bioinformatics analysis\n"
                "- **Google Gemini**: AI for natural language\n"
                "- **External APIs**: BLAST, Ensembl, UCSC, GTEx\n"
                "- **Streamlit**: Interactive interface"
            )
        else:
            st.info(
                "BioQuery Local uses:\n"
                "- **EMBOSS**: Bioinformatics analysis\n"
                "- **Ollama**: Local AI (no internet)\n"
                "- **Streamlit**: Interactive interface\n\n"
                "✓ Complete privacy - all processing is local"
            )
    
    # Main interface tabs
    if current_mode == 'cloud':
        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
            "🤖 Natural Language Query",
            "🔧 Manual Tool Selection",
            "🌐 Genome Browser Query",
            "📁 FASTA Upload & Batch",
            "📋 Command Log",
            "📖 Documentation"
        ])
    else:
        # Local mode: hide genome browser tab (requires internet)
        tab1, tab2, tab4, tab5, tab6 = st.tabs([
            "🤖 Natural Language Query",
            "🔧 Manual Tool Selection",
            "📁 FASTA Upload & Batch",
            "📋 Command Log",
            "📖 Documentation"
        ])
        tab3 = None  # No genome browser in local mode
    
    # TAB 1: Natural Language Query
    with tab1:
        st.subheader("Ask a bioinformatics question")
        st.write("Describe what you want to do with your sequences in natural language.")
        
        # Optional FASTA file upload (supports multiple files)
        with st.expander("📁 Optional: Upload FASTA file(s) to reference in your query"):
            uploaded_files = st.file_uploader(
                "Upload FASTA file(s) (optional):",
                type=["fasta", "fa", "faa", "fna", "txt"],
                help="Upload one or more FASTA files. Use for batch analysis or comparing sequences from different files.",
                key="nlp_fasta_upload",
                accept_multiple_files=True
            )
            
            if uploaded_files:
                all_sequences = []
                file_info = []
                
                for uploaded_file in uploaded_files:
                    content = uploaded_file.read().decode('utf-8')
                    sequences = parse_fasta(content)
                    
                    if sequences:
                        file_info.append((uploaded_file.name, len(sequences)))
                        # Tag sequences with their source file
                        for header, seq in sequences:
                            all_sequences.append((f"{uploaded_file.name}: {header}", seq))
                
                if all_sequences:
                    total_files = len(uploaded_files)
                    total_seqs = len(all_sequences)
                    st.success(f"✓ Loaded {total_seqs} sequence(s) from {total_files} file(s)")
                    
                    # Store in session state for use in queries
                    st.session_state['uploaded_sequences'] = all_sequences
                    st.session_state['uploaded_files_info'] = file_info
                    
                    # Show preview grouped by file
                    with st.expander("Preview sequences"):
                        for filename, seq_count in file_info:
                            st.markdown(f"**{filename}** ({seq_count} sequences)")
                        
                        st.markdown("---")
                        st.caption("First 3 sequences:")
                        for i, (header, seq) in enumerate(all_sequences[:3], 1):
                            st.text(f">{header}")
                            st.text(f"{seq[:80]}{'...' if len(seq) > 80 else ''}")
                        if len(all_sequences) > 3:
                            st.caption(f"... and {len(all_sequences) - 3} more")
                    
                    st.info("💡 Now you can reference these sequences in your query!\n\n"
                           "Examples:\n"
                           "- 'Find the GC content of the sequences in the file'\n"
                           "- 'Translate all sequences in the uploaded files'\n"
                           "- 'Make a dot plot comparing sequences from file 1 and file 2'\n"
                           "- 'Calculate protein stats for the sequences'")
                else:
                    st.warning("No sequences found. Make sure files are in FASTA format.")
        
        col1, col2 = st.columns([3, 1])
        with col1:
            user_query = st.text_area(
                "Your query:",
                placeholder="E.g., 'Translate this DNA sequence to protein: ATGAAATTTCCC' or 'Find GC content of the uploaded file'",
                height=100,
                key="nlp_query"
            )
        
        with col2:
            submit_button = st.button("🚀 Analyze", key="analyze_btn", use_container_width=True)
        
        if submit_button and user_query.strip():
            with st.spinner("Processing your query with AI..."):
                try:
                    # First, check if there's a long sequence in the query (FASTA or plain sequence)
                    # This prevents sending sequences to Gemini which triggers content filters
                    extracted_sequence = None
                    cleaned_query = user_query
                    
                    # Check for PDB structure queries first (bypass NLP entirely)
                    import re
                    pdb_match = re.search(r'\b(pdb|structure|protein structure)\s+([0-9][A-Za-z0-9]{3})\b', user_query, re.IGNORECASE)
                    if pdb_match:
                        pdb_id = pdb_match.group(2).upper()
                        st.markdown('<div class="tool-box">', unsafe_allow_html=True)
                        st.success(f"✓ Detected PDB structure query: {pdb_id}")
                        st.markdown('</div>', unsafe_allow_html=True)
                        
                        st.info("**Action:** Get protein structure information from PDB")
                        st.write(f"**PDB ID:** {pdb_id}")
                        
                        try:
                            result = emboss.run_tool('protein_structure', pdb_id=pdb_id)
                            st.markdown("#### Results:")
                            display_result(result)
                            
                            # Mark as handled
                            success = True
                            result = {'handled': True}
                        except Exception as e:
                            st.error(f"Error: {str(e)}")
                            success = False
                    
                    # Look for FASTA format (>header followed by sequence)
                    elif '>' in user_query:
                        import re
                        # Extract all FASTA sequences
                        fasta_pattern = r'>([^\n]+)\n([ATCGNatcgn\s]+?)(?=>|\Z)'
                        matches = re.findall(fasta_pattern, user_query, re.IGNORECASE | re.DOTALL)
                        
                        if matches:
                            # Get text before first > as the instruction
                            first_gt = user_query.index('>')
                            text_before = user_query[:first_gt].strip()
                            
                            # PRIORITY: Check for two-sequence tool keywords first (dotplot, align, compare)
                            two_seq_keywords = ['dot plot', 'dotplot', 'align', 'compare', 'similarity', 'match']
                            has_two_seq_request = any(kw.lower() in text_before.lower() for kw in two_seq_keywords)
                            
                            # Check if we have multiple sequences
                            has_multiple_seqs = len(matches) >= 2
                            
                            # Build cleaned query from context
                            if has_two_seq_request and has_multiple_seqs:
                                # User wants to compare/align/dotplot two sequences
                                cleaned_query = text_before if text_before else "make a dot plot of these two sequences"
                            else:
                                # Extract key context words for characterization
                                context_keywords = ['characterize', 'analyze', 'identify', 'what is', 'mystery', 
                                                  'genome browser', 'hg19', 'hg38', 'functional genomics']
                                found_keywords = [kw for kw in context_keywords if kw.lower() in text_before.lower()]
                                
                                if found_keywords:
                                    cleaned_query = f"Characterize this sequence using hg19"
                                else:
                                    cleaned_query = text_before if text_before else "analyze this sequence"
                            
                            # Extract sequences based on context
                            if has_two_seq_request and has_multiple_seqs:
                                # Extract both sequences for two-sequence tools
                                header1, seq1 = matches[0]
                                header2, seq2 = matches[1]
                                extracted_sequence = {
                                    'seq1': seq1.replace('\n', '').replace(' ', '').strip().upper(),
                                    'seq2': seq2.replace('\n', '').replace(' ', '').strip().upper()
                                }
                            else:
                                # Take first sequence for single-sequence tools
                                if matches:
                                    header, seq = matches[0]
                                    extracted_sequence = seq.replace('\n', '').replace(' ', '').strip().upper()
                    
                    # Look for long continuous sequences without FASTA format
                    if not extracted_sequence:
                        import re
                        
                        # Check if user wants dotplot/align with two sequences
                        two_seq_keywords = ['dot plot', 'dotplot', 'align', 'compare', 'similarity', 'match']
                        has_two_seq_request = any(kw.lower() in user_query.lower() for kw in two_seq_keywords)
                        
                        # Look for multiple sequence blocks (lines with gene IDs followed by sequences)
                        # Pattern: gene identifier line followed by nucleotide sequences
                        sequence_blocks = re.findall(r'(?:^|\n)([A-Z0-9]+(?:\.[0-9]+)?\s*\([^)]+\)[^\n]*)\n([ATCGNatcgn\s]+?)(?=\n[A-Z0-9]+(?:\.[0-9]+)?\s*\(|\Z)', 
                                                     user_query, re.MULTILINE | re.IGNORECASE)
                        
                        if has_two_seq_request and len(sequence_blocks) >= 2:
                            # Extract two sequences for dotplot/align
                            header1, seq1 = sequence_blocks[0]
                            header2, seq2 = sequence_blocks[1]
                            extracted_sequence = {
                                'seq1': ''.join(seq1.split()).upper(),
                                'seq2': ''.join(seq2.split()).upper()
                            }
                            # Extract instruction text before sequences
                            first_seq_pos = user_query.index(sequence_blocks[0][0])
                            cleaned_query = user_query[:first_seq_pos].strip()
                            if not cleaned_query:
                                cleaned_query = "make a dot plot of these two sequences"
                        else:
                            # Single sequence - find first long sequence
                            seq_match = re.search(r'\b([ATCGNatcgn]{30,})\b', user_query)
                            if seq_match:
                                extracted_sequence = seq_match.group(1).upper()
                                # Remove the sequence from the query
                                cleaned_query = user_query.replace(seq_match.group(0), '').strip()
                                if not cleaned_query:
                                    cleaned_query = "characterize this sequence using hg19"
                    
                    # Check if query references uploaded file
                    file_keywords = ['file', 'uploaded', 'fasta', 'sequences in', 'all sequences']
                    references_file = any(keyword in user_query.lower() for keyword in file_keywords)
                    has_uploaded = 'uploaded_sequences' in st.session_state
                    
                    # If we extracted a sequence from the query, send to NLP WITHOUT sequence data for speed
                    if extracted_sequence and not references_file:
                        st.markdown('<div class="tool-box">', unsafe_allow_html=True)
                        
                        # Handle both single sequence (string) and two sequences (dict)
                        if isinstance(extracted_sequence, dict):
                            st.success(f"✓ Detected two sequences: seq1 ({len(extracted_sequence['seq1'])} bp), seq2 ({len(extracted_sequence['seq2'])} bp)")
                        else:
                            st.success(f"✓ Detected sequence in query ({len(extracted_sequence)} bp)")
                        
                        st.markdown('</div>', unsafe_allow_html=True)
                        
                        # Parse with NLP using ONLY the cleaned query (no sequence data)
                        # This is MUCH faster - sequences get injected AFTER tool selection
                        success, result = nlp.parse_user_query(cleaned_query)
                        
                        if success:
                            # Inject the extracted sequence(s) into parameters
                            if 'steps' in result:
                                # Multi-step workflow
                                for step in result['steps']:
                                    if isinstance(extracted_sequence, dict):
                                        # Two sequences - always inject (overwrite placeholders)
                                        step['parameters']['seq1'] = extracted_sequence['seq1']
                                        step['parameters']['seq2'] = extracted_sequence['seq2']
                                    else:
                                        # Single sequence - always inject (overwrite placeholder)
                                        step['parameters']['sequence'] = extracted_sequence
                            elif 'parameters' in result:
                                # Single step
                                if isinstance(extracted_sequence, dict):
                                    # Two sequences - always inject (overwrite placeholders)
                                    result['parameters']['seq1'] = extracted_sequence['seq1']
                                    result['parameters']['seq2'] = extracted_sequence['seq2']
                                else:
                                    # Single sequence - always inject (overwrite placeholder)
                                    result['parameters']['sequence'] = extracted_sequence
                        else:
                            # NLP failed, fall back to simple keyword detection
                            query_lower = cleaned_query.lower()
                            detected_operations = []
                            detected_tools = set()
                            
                            # Define keyword mappings to tools
                            tool_keywords = [
                                # Two-sequence tools (check first)
                                ('dot plot', 'dotplot', {}, 'Create dot plot comparison'),
                                ('dotplot', 'dotplot', {}, 'Create dot plot comparison'),
                                ('align', 'align', {}, 'Align two sequences'),
                                # Single-sequence tools
                                ('gc content', 'gc_content', {}, 'Calculate GC content'),
                                ('reverse complement', 'reverse', {}, 'Reverse complement'),
                                ('isoelectric point', 'iep', {}, 'Calculate isoelectric point'),
                                ('molecular weight', 'pepstats', {}, 'Get protein statistics including molecular weight'),
                                ('protein stats', 'pepstats', {}, 'Get protein statistics'),
                                ('open reading', 'orf', {}, 'Find open reading frames'),
                                ('gc', 'gc_content', {}, 'Calculate GC content'),
                                ('translate', 'translate', {}, 'Translate to protein'),
                                ('reverse', 'reverse', {}, 'Reverse complement'),
                                ('isoelectric', 'iep', {}, 'Calculate isoelectric point'),
                                ('pepstats', 'pepstats', {}, 'Get protein statistics'),
                                ('orf', 'orf', {}, 'Find open reading frames'),
                            ]
                            
                            # Check which operations are mentioned
                            for keyword, tool, params, explanation in tool_keywords:
                                if keyword in query_lower and tool not in detected_tools:
                                    detected_operations.append({
                                        'tool': tool,
                                        'parameters': params.copy(),
                                        'explanation': explanation
                                    })
                                    detected_tools.add(tool)
                            
                            if detected_operations:
                                st.info(f"**Detected {len(detected_operations)} operation(s)**")
                                
                                # Process each detected operation
                                for op_idx, op_info in enumerate(detected_operations, 1):
                                    if len(detected_operations) > 1:
                                        st.markdown(f"### Operation {op_idx}: {op_info['explanation']}")
                                    else:
                                        st.info(f"**Action:** {op_info['explanation']}")
                                    
                                    st.write(f"**Tool:** {op_info['tool']}")
                                    
                                    try:
                                        # Execute tool with the extracted sequence(s)
                                        params = op_info['parameters'].copy()
                                        
                                        # Check if this is a two-sequence tool
                                        if op_info['tool'] in ['dotplot', 'align'] and isinstance(extracted_sequence, dict):
                                            params['seq1'] = extracted_sequence['seq1']
                                            params['seq2'] = extracted_sequence['seq2']
                                        else:
                                            params['sequence'] = extracted_sequence
                                        
                                        result = emboss.run_tool(op_info['tool'], **params)
                                        
                                        st.markdown("#### Results:")
                                        display_result(result)
                                        
                                    except Exception as e:
                                        st.error(f"Error: {str(e)}")
                                    
                                    if len(detected_operations) > 1 and op_idx < len(detected_operations):
                                        st.markdown("---")
                            else:
                                st.warning("⚠️ Couldn't detect specific operations. Try: 'characterize this sequence', 'get gc content', 'translate', 'reverse complement'")
                            
                            # Mark as handled so we don't show error message below
                            success = True
                            result = {'handled': True}
                    
                    elif references_file and has_uploaded:
                        # Process query for each sequence in uploaded file
                        sequences = st.session_state['uploaded_sequences']
                        filename = st.session_state.get('uploaded_filename', 'uploaded file')
                        
                        st.markdown('<div class="tool-box">', unsafe_allow_html=True)
                        st.success(f"✓ Processing query for {len(sequences)} sequences from {filename}")
                        st.markdown('</div>', unsafe_allow_html=True)
                        
                        # Check if this is a two-sequence comparison tool
                        query_lower = user_query.lower()
                        two_seq_tools = ['align', 'water', 'stretcher', 'matcher', 'dotplot', 'dot plot', 
                                        'compare', 'alignment', 'similarity', 'match']
                        is_comparison = any(tool in query_lower for tool in two_seq_tools)
                        
                        # If we have exactly 2 sequences and query mentions comparison, handle as two-sequence comparison
                        if len(sequences) == 2 and is_comparison:
                            seq1_header, seq1 = sequences[0]
                            seq2_header, seq2 = sequences[1]
                            
                            st.info(f"**Detected two-sequence comparison:**")
                            st.write(f"• Sequence 1: {seq1_header} ({len(seq1)} bp)")
                            st.write(f"• Sequence 2: {seq2_header} ({len(seq2)} bp)")
                            
                            # Parse query WITHOUT sequences to get tool name
                            success, result = nlp.parse_user_query(user_query)
                            
                            if success:
                                tool_name = result.get('tool')
                                explanation = result.get('explanation', '')
                                
                                st.info(f"**Action:** {explanation}")
                                st.write(f"**Tool:** {tool_name}")
                                
                                # Run comparison tool
                                try:
                                    result_data = emboss.run_tool(tool_name, seq1=seq1, seq2=seq2)
                                    
                                    st.markdown("---")
                                    st.subheader("📊 Results")
                                    display_result(result_data)
                                    
                                except Exception as e:
                                    st.error(f"Error running {tool_name}: {str(e)}")
                            else:
                                st.error("Could not determine which tool to use")
                            
                            # Mark as handled
                            success = True
                        else:
                            # Regular batch processing for single-sequence tools or more than 2 sequences
                            # For file-based queries, try to detect operations directly to avoid Gemini filters
                            detected_operations = []
                            detected_tools = set()  # Track which tools we've already added
                            
                            # Define keyword mappings to tools (order matters - longer phrases first)
                            tool_keywords = [
                                ('gc content', 'gc_content', {}, 'Calculate GC content'),
                                ('reverse complement', 'reverse', {}, 'Reverse complement'),
                                ('isoelectric point', 'iep', {}, 'Calculate isoelectric point'),
                                ('molecular weight', 'pepstats', {}, 'Get protein statistics including molecular weight'),
                                ('protein stats', 'pepstats', {}, 'Get protein statistics'),
                                ('open reading', 'orf', {}, 'Find open reading frames'),
                                ('gc', 'gc_content', {}, 'Calculate GC content'),
                                ('translate', 'translate', {}, 'Translate to protein'),
                                ('reverse', 'reverse', {}, 'Reverse complement'),
                                ('isoelectric', 'iep', {}, 'Calculate isoelectric point'),
                                ('pepstats', 'pepstats', {}, 'Get protein statistics'),
                                ('orf', 'orf', {}, 'Find open reading frames'),
                            ]
                            
                            # Check which operations are mentioned (avoid duplicates)
                            for keyword, tool, params, explanation in tool_keywords:
                                if keyword in query_lower and tool not in detected_tools:
                                    detected_operations.append({
                                        'tool': tool,
                                        'parameters': params.copy(),
                                        'explanation': explanation
                                    })
                                    detected_tools.add(tool)
                            
                            # If we detected operations directly, use them instead of calling Gemini
                            if detected_operations:
                                st.info(f"**Detected {len(detected_operations)} operation(s)**")
                                
                                # Process each detected operation
                                for op_idx, op_info in enumerate(detected_operations, 1):
                                    if len(detected_operations) > 1:
                                        st.markdown(f"### Operation {op_idx}: {op_info['explanation']}")
                                    else:
                                        st.info(f"**Action:** {op_info['explanation']}")
                                    
                                    st.write(f"**Tool:** {op_info['tool']}")
                                    
                                    # Progress tracking
                                    progress_bar = st.progress(0)
                                    status_text = st.empty()
                                    
                                    all_results = []
                                    
                                    # Process each sequence
                                    for i, (header, seq) in enumerate(sequences):
                                        progress = (i + 1) / len(sequences)
                                        progress_bar.progress(progress)
                                        status_text.text(f"Processing {i+1}/{len(sequences)}: {header[:50]}...")
                                        
                                        try:
                                            # Override sequence parameter with current sequence
                                            params = op_info['parameters'].copy()
                                            params['sequence'] = seq
                                            
                                            # Execute tool
                                            result = emboss.run_tool(op_info['tool'], **params)
                                            all_results.append({
                                                'header': header,
                                                'result': result
                                            })
                                        except Exception as e:
                                            all_results.append({
                                                'header': header,
                                                'result': f"Error: {str(e)}"
                                            })
                                    
                                    progress_bar.progress(1.0)
                                    status_text.text("✓ Complete!")
                                    
                                    # Display results
                                    st.markdown("#### Results:")
                                    for item in all_results:
                                        with st.expander(f"📄 {item['header']}", expanded=(len(all_results) == 1)):
                                            display_result(item['result'])
                                    
                                    if len(detected_operations) > 1:
                                        st.markdown("---")
                            
                            # Only proceed with Gemini parsing if we didn't detect operations directly
                            if not detected_operations:
                                # Fall back to Gemini if we couldn't detect operations directly
                                # Simplify query by removing file references and adding a placeholder sequence
                                # This prevents Gemini from seeing actual sequence data which can trigger filters
                                simplified_query = user_query.lower()
                                
                                # Remove file references
                                for keyword in ['in the file', 'in the uploaded file', 'from the file', 'in file', 
                                               'from file', 'in uploaded file', 'of the file', 'the sequences',
                                               'of the uploaded', 'the uploaded', 'uploaded', 'of the sequence',
                                               'the sequence in']:
                                    simplified_query = simplified_query.replace(keyword, '')
                                
                                # Clean up extra spaces and "and" at the beginning
                                simplified_query = ' '.join(simplified_query.split())
                                if simplified_query.startswith('and '):
                                    simplified_query = simplified_query[4:]
                                
                                # For multi-operation queries, split them into numbered questions
                                # This uses the existing multi-question parsing which is more reliable
                                if ' and the ' in simplified_query or ' and ' in simplified_query:
                                    # Split on "and the" or "and"
                                    if ' and the ' in simplified_query:
                                        parts = simplified_query.split(' and the ')
                                    else:
                                        parts = simplified_query.split(' and ')
                                    
                                    # Format as numbered questions if we have multiple parts
                                    if len(parts) >= 2:
                                        numbered_questions = []
                                        for i, part in enumerate(parts, 1):
                                            part = part.strip()
                                            # Ensure each part has a verb
                                            if not any(verb in part for verb in ['find', 'get', 'calculate', 'show', 'translate', 'determine']):
                                                part = f"find {part}"
                                            numbered_questions.append(f"{i}. {part}")
                                        simplified_query = " ".join(numbered_questions)
                                
                                # Add a generic sequence placeholder so Gemini knows it's a sequence operation
                                if "sequence" not in simplified_query:
                                    simplified_query = simplified_query.strip() + " of ATGC"
                                
                                # Parse the simplified query to understand what tool to use
                                success, result = nlp.parse_user_query(simplified_query)
                                
                                if success and 'tool' in result:
                                    tool_name = result.get('tool')
                                    parameters = result.get('parameters', {})
                                    explanation = result.get('explanation', '')
                                    
                                    st.info(f"**Action:** {explanation}")
                                    st.write(f"**Tool:** {tool_name}")
                                    
                                    # Progress tracking
                                    progress_bar = st.progress(0)
                                    status_text = st.empty()
                                    
                                    all_results = []
                                    
                                    # Process each sequence
                                    for i, (header, seq) in enumerate(sequences):
                                        progress = (i + 1) / len(sequences)
                                        progress_bar.progress(progress)
                                        status_text.text(f"Processing {i+1}/{len(sequences)}: {header[:50]}...")
                                        
                                        try:
                                            # Override sequence parameter with current sequence
                                            params = parameters.copy()
                                            params['sequence'] = seq
                                            
                                            # Execute tool
                                            if tool_name == 'genome_query':
                                                result_text = emboss.query_ucsc_genome(
                                                    params.get('genome', ''),
                                                    params.get('chrom', ''),
                                                    int(params.get('start', 0)),
                                                    int(params.get('end', 0))
                                                )
                                            elif tool_name == 'gene_query':
                                                result_text = emboss.query_gene_info(
                                                    params.get('gene_name', ''),
                                                    params.get('genome', 'hg38'),
                                                    params.get('track', 'gencode')
                                                )
                                            else:
                                                result_text = emboss.run_tool(tool_name, **params)
                                            
                                            all_results.append((header, result_text))
                                            
                                        except Exception as e:
                                            all_results.append((header, f"ERROR: {str(e)}"))
                                    
                                    status_text.text("✓ Analysis complete!")
                                    progress_bar.progress(1.0)
                                    
                                    # Display results
                                    st.markdown("---")
                                    st.subheader("📊 Results")
                                    
                                    for header, result_text in all_results:
                                        with st.expander(f"📄 {header}"):
                                            st.code(result_text, language="text")
                                    
                                    # Download button
                                    combined_results = ""
                                    for header, result_text in all_results:
                                        combined_results += f"{'='*60}\n"
                                        combined_results += f"Sequence: {header}\n"
                                        combined_results += f"Tool: {tool_name}\n"
                                        combined_results += f"Query: {user_query}\n"
                                        combined_results += f"{'='*60}\n"
                                        combined_results += f"{result_text}\n\n"
                                    
                                    st.download_button(
                                        label="📥 Download All Results",
                                        data=combined_results,
                                        file_name=f"{tool_name}_results.txt",
                                        mime="text/plain",
                                        type="primary"
                                    )
                                else:
                                    st.error("Could not parse query. Please be more specific about what analysis to perform.")
                                    if 'error' in result:
                                        st.error(f"Could not parse query: {result['error']}")
                    
                    elif references_file and not has_uploaded:
                        success = False  # Initialize for error handling below
                        st.warning("⚠️ Your query references a file, but no FASTA file has been uploaded. Please upload a file first using the section above.")
                    
                    else:
                        # Normal query processing (no file reference)
                        # Parse the query
                        success, result = nlp.parse_user_query(user_query)
                    
                    # Only process regular results if 'success' was defined (not using keyword detection)
                    if 'success' in locals() and success:
                        # Skip if already handled by keyword detection
                        if isinstance(result, dict) and result.get('handled'):
                            pass  # Already processed above
                        # Check if this is multiple unrelated questions
                        elif result.get('type') == 'multiple_questions':
                            # Handle multiple questions
                            st.markdown('<div class="tool-box">', unsafe_allow_html=True)
                            st.write(f"**Detected {result['total']} questions**")
                            st.markdown('</div>', unsafe_allow_html=True)
                            
                            # Process each question
                            for q_result in result['questions']:
                                q_num = q_result['question_number']
                                q_text = q_result['question_text']
                                
                                st.markdown(f"### Question {q_num}")
                                st.info(f"**Q{q_num}:** {q_text}")
                                
                                if not q_result.get('success'):
                                    st.error(f"Failed to parse: {q_result.get('error')}")
                                    continue
                                
                                parsed = q_result['parsed']
                                
                                # Handle multi-step for this question
                                if 'steps' in parsed:
                                    with st.expander(f"Q{q_num} - Multi-step workflow", expanded=True):
                                        st.write(f"**Workflow:** {parsed.get('explanation', 'N/A')}")
                                        
                                        use_prev = parsed.get('use_previous_result', False)
                                        success_multi, results = emboss.execute_multi_step(
                                            parsed['steps'],
                                            use_previous_result=use_prev
                                        )
                                        
                                        if success_multi:
                                            formatted_output = emboss.format_multi_step_results(
                                                results,
                                                explanation=parsed.get('explanation', '')
                                            )
                                            st.code(formatted_output, language="text")
                                        else:
                                            st.error("Some steps failed")
                                            for r in results:
                                                if not r.get('success'):
                                                    st.write(f"Step {r.get('step')}: {r.get('error')}")
                                
                                # Handle single-step for this question
                                else:
                                    tool_name = parsed.get('tool')
                                    parameters = parsed.get('parameters', {})
                                    
                                    with st.expander(f"Q{q_num} - {tool_name}", expanded=True):
                                        st.write(f"**Tool:** {tool_name}")
                                        st.write(f"**Explanation:** {parsed.get('explanation', '')}")
                                        
                                        # Execute tool
                                        try:
                                            if tool_name == 'genome_query':
                                                analysis_result = emboss.query_ucsc_genome(
                                                    parameters.get('genome', ''),
                                                    parameters.get('chrom', ''),
                                                    int(parameters.get('start', 0)),
                                                    int(parameters.get('end', 0))
                                                )
                                            elif tool_name == 'gene_query':
                                                analysis_result = emboss.query_gene_info(
                                                    parameters.get('gene_name', ''),
                                                    parameters.get('genome', 'hg38'),
                                                    parameters.get('track', 'gencode')
                                                )
                                            else:
                                                analysis_result = emboss.run_tool(tool_name, **parameters)
                                            
                                            st.code(analysis_result, language="text")
                                        except Exception as e:
                                            st.error(f"Error: {str(e)}")
                                
                                st.markdown("---")
                            
                            # Overall download button
                            all_results = "\n\n".join([
                                f"=== Question {q['question_number']} ===\n{q['question_text']}\n\n[Results would go here]"
                                for q in result['questions']
                            ])
                            st.download_button(
                                label="📥 Download All Results",
                                data=all_results,
                                file_name="multiple_questions_results.txt",
                                mime="text/plain"
                            )
                        
                        # Check if this is a multi-step query
                        elif 'steps' in result:
                            # Multi-step query
                            st.markdown('<div class="tool-box">', unsafe_allow_html=True)
                            st.write(f"**Multi-Step Workflow:** {result.get('explanation', 'N/A')}")
                            st.write(f"**Number of Steps:** {len(result['steps'])}")
                            
                            for idx, step in enumerate(result['steps'], 1):
                                st.write(f"  Step {idx}: {step.get('tool')} with {step.get('parameters', {})}")
                            st.markdown('</div>', unsafe_allow_html=True)
                            
                            # Execute multi-step workflow
                            with st.spinner("Executing multi-step workflow..."):
                                try:
                                    use_prev = result.get('use_previous_result', False)
                                    success_multi, results = emboss.execute_multi_step(
                                        result['steps'],
                                        use_previous_result=use_prev
                                    )
                                    
                                    if success_multi:
                                        st.markdown('<div class="success-box">', unsafe_allow_html=True)
                                        st.success("Multi-step analysis completed!")
                                        
                                        # Format and display results
                                        formatted_output = emboss.format_multi_step_results(
                                            results,
                                            explanation=result.get('explanation', '')
                                        )
                                        
                                        with st.expander("📊 Detailed Results"):
                                            st.code(formatted_output, language="text")
                                        
                                        st.markdown('</div>', unsafe_allow_html=True)
                                        
                                        # Download button
                                        st.download_button(
                                            label="📥 Download Results",
                                            data=formatted_output,
                                            file_name="multi_step_results.txt",
                                            mime="text/plain"
                                        )
                                    else:
                                        st.markdown('<div class="error-box">', unsafe_allow_html=True)
                                        st.error("One or more steps failed")
                                        for result_item in results:
                                            if not result_item.get('success'):
                                                st.write(f"Step {result_item.get('step')}: {result_item.get('error')}")
                                        st.markdown('</div>', unsafe_allow_html=True)
                                
                                except Exception as e:
                                    st.markdown('<div class="error-box">', unsafe_allow_html=True)
                                    st.error(f"Error executing multi-step workflow: {str(e)}")
                                    st.markdown('</div>', unsafe_allow_html=True)
                        
                        else:
                            # Single-step query (original code)
                            tool_name = result.get('tool')
                            parameters = result.get('parameters', {})
                            explanation = result.get('explanation', '')
                            
                            # Display parsed information
                            st.markdown('<div class="tool-box">', unsafe_allow_html=True)
                            st.write(f"**Tool Selected:** {tool_name}")
                            st.write(f"**Explanation:** {explanation}")
                            st.write(f"**Parameters:** {parameters}")
                            st.markdown('</div>', unsafe_allow_html=True)
                            
                            # Run the tool
                            with st.spinner(f"Running {tool_name}..."):
                                try:
                                    summary = None
                                    
                                    # Handle genome queries specially
                                    if tool_name == 'genome_query':
                                        sequence = emboss.query_ucsc_genome(
                                            parameters.get('genome', ''),
                                            parameters.get('chrom', ''),
                                            int(parameters.get('start', 0)),
                                            int(parameters.get('end', 0))
                                        )
                                        analysis_result = sequence
                                    # Handle gene queries specially
                                    elif tool_name == 'gene_query':
                                        gene_info = emboss.query_gene_info(
                                            parameters.get('gene_name', ''),
                                            parameters.get('genome', 'hg38'),
                                            parameters.get('track', 'gencode')
                                        )
                                        
                                        # Get natural language summary
                                        with st.spinner("Generating natural language summary..."):
                                            summary = emboss.summarize_gene_info(gene_info)
                                        
                                        analysis_result = gene_info
                                    else:
                                        analysis_result = emboss.run_tool(tool_name, **parameters)
                                    
                                    st.markdown('<div class="success-box">', unsafe_allow_html=True)
                                    st.success("Analysis completed!")
                                    
                                    # Show natural language summary for gene queries
                                    if summary:
                                        st.subheader("📝 Answer:")
                                        st.info(summary)
                                    
                                    # Expandable section for detailed results
                                    with st.expander("📊 Detailed Results"):
                                        st.subheader("Detailed Information:")
                                        st.code(analysis_result, language="text")
                                    
                                    st.markdown('</div>', unsafe_allow_html=True)
                                    
                                    # Download button
                                    st.download_button(
                                        label="📥 Download Results",
                                        data=analysis_result,
                                        file_name=f"{tool_name}_results.txt",
                                        mime="text/plain"
                                    )
                                except Exception as e:
                                    st.markdown('<div class="error-box">', unsafe_allow_html=True)
                                    st.error(f"Error running tool: {str(e)}")
                                    st.markdown('</div>', unsafe_allow_html=True)
                    else:
                        # Only show error if result variable exists and parsing failed
                        if 'result' in locals():
                            st.markdown('<div class="error-box">', unsafe_allow_html=True)
                            if isinstance(result, dict):
                                st.error(f"Could not parse query: {result.get('error', 'Unknown error')}")
                            else:
                                st.error(f"Could not parse query: {result}")
                            st.markdown('</div>', unsafe_allow_html=True)
                
                except Exception as e:
                    st.markdown('<div class="error-box">', unsafe_allow_html=True)
                    st.error(f"Error processing query: {str(e)}")
                    st.markdown('</div>', unsafe_allow_html=True)
                    st.write(traceback.format_exc())
    
    # TAB 2: Manual Tool Selection
    with tab2:
        st.subheader("Select a specific tool to run")
        
        col1, col2 = st.columns(2)
        
        with col1:
            selected_tool = st.selectbox(
                "Choose a tool:",
                list(emboss.get_available_tools().keys()),
                key="manual_tool"
            )
            st.info(emboss.get_available_tools()[selected_tool])
        
        with col2:
            st.write("Enter sequence(s) below:")
        
        # Tool-specific inputs
        if selected_tool == "translate":
            sequence = st.text_area("DNA Sequence:", key="trans_seq", height=80)
            frame = st.slider("Reading Frame:", 1, 3, 1)
            if st.button("Translate", key="trans_btn"):
                with st.spinner("Translating..."):
                    result = emboss.translate_sequence(sequence, frame)
                    st.code(result, language="text")
        
        elif selected_tool == "reverse":
            sequence = st.text_area("DNA Sequence:", key="rev_seq", height=80)
            if st.button("Reverse Complement", key="rev_btn"):
                with st.spinner("Reversing..."):
                    result = emboss.reverse_complement(sequence)
                    st.code(result, language="text")
        
        elif selected_tool == "orf":
            sequence = st.text_area("DNA Sequence:", key="orf_seq", height=80)
            min_size = st.slider("Minimum ORF Size (bp):", 10, 500, 100)
            if st.button("Find ORFs", key="orf_btn"):
                with st.spinner("Finding ORFs..."):
                    result = emboss.find_orfs(sequence, min_size)
                    st.code(result, language="text")
        
        elif selected_tool == "info":
            sequence = st.text_area("Sequence:", key="info_seq", height=80)
            if st.button("Get Info", key="info_btn"):
                with st.spinner("Analyzing..."):
                    result = emboss.get_sequence_info(sequence)
                    st.code(result, language="text")
        
        elif selected_tool == "sixframe":
            sequence = st.text_area("DNA Sequence:", key="six_seq", height=80)
            if st.button("Show Six Frames", key="six_btn"):
                with st.spinner("Translating..."):
                    result = emboss.get_six_frame_translation(sequence)
                    st.code(result, language="text")
        
        elif selected_tool == "align":
            col1, col2 = st.columns(2)
            with col1:
                seq1 = st.text_area("Sequence 1:", key="align_seq1", height=60)
            with col2:
                seq2 = st.text_area("Sequence 2:", key="align_seq2", height=60)
            if st.button("Align", key="align_btn"):
                with st.spinner("Aligning..."):
                    result = emboss.align_sequences(seq1, seq2)
                    st.code(result, language="text")
        
        elif selected_tool == "restriction":
            sequence = st.text_area("DNA Sequence:", key="rest_seq", height=80)
            enzyme = st.text_input("Enzyme (optional):", key="rest_enz")
            if st.button("Find Sites", key="rest_btn"):
                with st.spinner("Searching..."):
                    result = emboss.find_restriction_sites(sequence, enzyme if enzyme else None)
                    st.code(result, language="text")
    
    # TAB 3: Genome Browser Query (only in cloud mode)
    if tab3 is not None:
        with tab3:
            st.subheader("Query genome sequences via UCSC (no download needed)")
            st.info("Stream genomic data directly from UCSC Genome Browser - no storage required!")
            
            col1, col2 = st.columns(2)
            
            with col1:
                genome = st.selectbox(
                    "Select genome assembly:",
                    ["hg38", "hg37", "mm10", "mm9", "dm6", "dm3", "ce10", "sacCer3"],
                    help="Human (hg), Mouse (mm), Fly (dm), Worm (ce), Yeast (sac)"
                )
                
                chrom = st.text_input(
                    "Chromosome:",
                    value="chr1",
                    placeholder="e.g., chr1, chrX, chr2",
                    help="Include 'chr' prefix"
                )
            
            with col2:
                start_pos = st.number_input(
                    "Start position (0-based):",
                    value=1000000,
                    min_value=0,
                    step=1000,
                    help="Start of genomic region"
                )
                
                end_pos = st.number_input(
                    "End position (0-based):",
                    value=1001000,
                    min_value=0,
                    step=1000,
                    help="End of genomic region"
                )
            
            if st.button("🔍 Query UCSC", key="ucsc_btn"):
                if start_pos >= end_pos:
                    st.error("Start position must be less than end position")
                else:
                    region_size = end_pos - start_pos
                    if region_size > 1000000:
                        st.warning(f"Large region ({region_size:,} bp). This may take a moment...")
                    
                    with st.spinner(f"Querying UCSC for {genome}:{chrom}:{start_pos}-{end_pos}..."):
                        try:
                            sequence = emboss.query_ucsc_genome(genome, chrom, int(start_pos), int(end_pos))
                            
                            if "Error" not in sequence:
                                st.markdown('<div class="success-box">', unsafe_allow_html=True)
                                st.success(f"Successfully retrieved {region_size:,} bp")
                                st.code(sequence, language="text")
                                
                                # Store in session for later use
                                st.session_state.last_sequence = sequence
                                
                                # Option to analyze
                                st.subheader("Quick Analysis")
                                col1, col2, col3 = st.columns(3)
                                
                                with col1:
                                    if st.button("📊 Get Info", key="ucsc_info"):
                                        info = emboss.get_sequence_info(sequence.split('\\n', 1)[1])
                                        st.code(info, language="text")
                                
                                with col2:
                                    if st.button("🧬 GC Content", key="ucsc_gc"):
                                        gc = emboss.calculate_gc_content(sequence.split('\\n', 1)[1])
                                        st.code(gc, language="text")
                                
                                with col3:
                                    if st.button("📥 Download", key="ucsc_down"):
                                        st.download_button(
                                            label="Download FASTA",
                                            data=sequence,
                                            file_name=f"{genome}_{chrom}_{start_pos}_{end_pos}.fasta",
                                            mime="text/plain"
                                        )
                                
                                st.markdown('</div>', unsafe_allow_html=True)
                            else:
                                st.markdown('<div class="error-box">', unsafe_allow_html=True)
                                st.error(sequence)
                                st.markdown('</div>', unsafe_allow_html=True)
                        
                        except Exception as e:
                            st.markdown('<div class="error-box">', unsafe_allow_html=True)
                            st.error(f"Error querying UCSC: {str(e)}")
                            st.markdown('</div>', unsafe_allow_html=True)
            
            st.markdown("---")
            st.subheader("Available Genomes")
            st.write("""
            - **hg38** / **hg37**: Human genome (GRCh38 / GRCh37)
            - **mm10** / **mm9**: Mouse genome (GRCm38 / GRCm37)
            - **dm6** / **dm3**: Drosophila (fruit fly)
            - **ce10**: Caenorhabditis elegans (nematode)
            - **sacCer3**: Saccharomyces cerevisiae (yeast)
            """)
    
    # TAB 4: FASTA File Upload & Batch Analysis
    with tab4:
        st.subheader("Upload and analyze FASTA files")
        
        st.info("Upload FASTA file(s) to analyze multiple sequences. For single file: batch processing. For multiple files: compare sequences across files.")
        
        # File upload (supports multiple)
        uploaded_files = st.file_uploader(
            "Choose FASTA file(s):",
            type=["fasta", "fa", "faa", "fna", "txt"],
            help="Upload one or more files containing sequences in FASTA format",
            accept_multiple_files=True
        )
        
        if uploaded_files:
            num_files = len(uploaded_files)
            st.success(f"✓ {num_files} file(s) uploaded")
            
            # Read and parse all files
            try:
                all_sequences = []
                file_sequence_map = {}
                
                for uploaded_file in uploaded_files:
                    content = uploaded_file.read().decode('utf-8')
                    sequences = parse_fasta(content)
                    
                    if sequences:
                        file_sequence_map[uploaded_file.name] = sequences
                        # Tag sequences with file name for clarity
                        for header, seq in sequences:
                            all_sequences.append((f"[{uploaded_file.name}] {header}", seq))
                
                if not all_sequences:
                    st.warning("No sequences found in files. Make sure they're in FASTA format (headers start with '>').")
                else:
                    total_seqs = len(all_sequences)
                    st.info(f"Found {total_seqs} sequence(s) across {num_files} file(s)")
                    
                    # Show preview grouped by file
                    with st.expander("📄 Preview sequences"):
                        for filename, sequences in file_sequence_map.items():
                            st.markdown(f"**{filename}** ({len(sequences)} sequences)")
                            for i, (header, seq) in enumerate(sequences[:2], 1):  # Show first 2 per file
                                st.text(f"  >{header}")
                                st.text(f"  {seq[:80]}{'...' if len(seq) > 80 else ''}")
                            if len(sequences) > 2:
                                st.caption(f"  ... and {len(sequences) - 2} more from this file")
                            st.markdown("---")
                    
                    # Special handling for 2 files (comparison mode)
                    if num_files == 2 and len(file_sequence_map) == 2:
                        st.markdown("### 🔀 Two-File Comparison Mode")
                        st.info("With 2 files, you can compare sequences between them!")
                        
                        file_names = list(file_sequence_map.keys())
                        if len(file_names) >= 2:
                            file1_seqs = file_sequence_map[file_names[0]]
                            file2_seqs = file_sequence_map[file_names[1]]
                        else:
                            st.warning("One or more files contain no valid sequences.")
                            file1_seqs = []
                            file2_seqs = []
                        
                        comparison_tool = st.selectbox(
                            "Comparison tool:",
                            ["dotplot", "align", "water", "None - batch process each file separately"],
                            help="Choose a tool to compare sequences from the two files"
                        )
                        
                        if comparison_tool != "None - batch process each file separately":
                            col1, col2 = st.columns(2)
                            with col1:
                                seq1_idx = st.selectbox(
                                    f"Sequence from {file_names[0]}:",
                                    range(len(file1_seqs)),
                                    format_func=lambda i: f"{file1_seqs[i][0][:50]}..."
                                )
                            with col2:
                                seq2_idx = st.selectbox(
                                    f"Sequence from {file_names[1]}:",
                                    range(len(file2_seqs)),
                                    format_func=lambda i: f"{file2_seqs[i][0][:50]}..."
                                )
                            
                            if st.button("🔬 Compare Selected Sequences", type="primary"):
                                seq1_header, seq1 = file1_seqs[seq1_idx]
                                seq2_header, seq2 = file2_seqs[seq2_idx]
                                
                                st.markdown("#### Comparison Results:")
                                st.write(f"**Sequence 1:** {seq1_header}")
                                st.write(f"**Sequence 2:** {seq2_header}")
                                st.write(f"**Tool:** {comparison_tool}")
                                
                                try:
                                    with st.spinner(f"Running {comparison_tool}..."):
                                        result = emboss.run_tool(comparison_tool, seq1=seq1, seq2=seq2)
                                    
                                    # Debug output
                                    if result.startswith("IMAGE_FILE:"):
                                        st.success("✓ Image generated!")
                                    else:
                                        st.info(f"Result type: {type(result)}, starts with: {result[:50] if result else 'empty'}")
                                    
                                    display_result(result)
                                except Exception as e:
                                    st.error(f"Error: {str(e)}")
                                    import traceback
                                    st.code(traceback.format_exc())
                            
                            st.markdown("---")
                    
                    st.markdown("### 📊 Batch Processing Mode")
                    
                    # Tool selection for batch processing
                    st.subheader("Select tool to apply")
                    
                    col1, col2 = st.columns([2, 1])
                    
                    with col1:
                        # Get all available tools dynamically
                        available_tools = list(emboss.get_available_tools().keys())
                        batch_tool = st.selectbox(
                            "Choose analysis tool:",
                            available_tools,
                            key="batch_tool"
                        )
                        
                        # Show tool description
                        st.caption(emboss.get_available_tools()[batch_tool])
                    
                    with col2:
                        # Tool-specific parameters
                        params = {}
                        if batch_tool == "translate":
                            params['frame'] = st.slider("Reading frame:", 1, 3, 1, key="batch_frame")
                        elif batch_tool == "orf":
                            params['min_size'] = st.slider("Min ORF size:", 10, 500, 100, key="batch_orf")
                        elif batch_tool == "align":
                            st.info("Note: Alignment requires 2 sequences. Will align consecutive pairs.")
                    
                    # Run analysis button
                    if st.button("🚀 Analyze All Sequences", type="primary", key="batch_run"):
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        results_container = st.container()
                        
                        batch_results = []
                        
                        with results_container:
                            for i, (header, seq) in enumerate(all_sequences):
                                progress = (i + 1) / len(all_sequences)
                                progress_bar.progress(progress)
                                status_text.text(f"Processing {i+1}/{len(all_sequences)}: {header[:50]}...")
                                
                                try:
                                    # Use the generic run_tool method with sequence parameter
                                    result = emboss.run_tool(batch_tool, sequence=seq, **params)
                                    batch_results.append((header, result))
                                    
                                except Exception as e:
                                    batch_results.append((header, f"ERROR: {str(e)}"))
                            
                            status_text.text("✓ Analysis complete!")
                            progress_bar.progress(1.0)
                        
                        # Display results
                        st.markdown("---")
                        st.subheader("📊 Results")
                        
                        for header, result in batch_results:
                            with st.expander(f"📄 {header}"):
                                display_result(result)
                        
                        # Download all results
                        combined_results = ""
                        for header, result in batch_results:
                            combined_results += f"{'='*60}\n"
                            combined_results += f"Sequence: {header}\n"
                            combined_results += f"Tool: {batch_tool}\n"
                            combined_results += f"{'='*60}\n"
                            combined_results += f"{result}\n\n"
                        
                        st.download_button(
                            label="📥 Download All Results",
                            data=combined_results,
                            file_name=f"{batch_tool}_batch_results.txt",
                            mime="text/plain",
                            type="primary"
                        )
            
            except Exception as e:
                st.error(f"Error processing file: {str(e)}")
                st.write(traceback.format_exc())
        
        else:
            # Show example FASTA format
            st.markdown("### FASTA Format Example")
            st.code("""
>sequence1 description here
ATGAAATTTCCCGGGAAATTTAAAGGG
AAATTTCCCGGG
>sequence2 another description
ATGCCCAAAGGGTTTTAA
>sequence3
GCTAGCTAGCTAGCTA
            """.strip(), language="text")
            
            st.markdown("### Supported Tools")
            st.write("""
            **All EMBOSS tools are supported in batch mode!** Including:
            
            **Common DNA Tools:**
            - **translate**: Convert DNA to protein
            - **reverse**: Reverse complement of DNA
            - **gc**: Calculate GC content percentage  
            - **orf**: Find open reading frames
            - **sixframe**: Show all 6 translation frames
            - **restriction**: Find restriction enzyme sites
            - **cusp**: Codon usage statistics
            
            **Protein Analysis:**
            - **pepstats**: Protein statistics (molecular weight, amino acid composition)
            - **iep**: Calculate isoelectric point
            - **info**: General sequence information
            
            **And 250+ more EMBOSS tools** dynamically discovered from your system!
            
            Simply upload your FASTA file and select any tool from the dropdown.
            """)
    
    # TAB 5: Command Log (NEW - Shows detailed execution history)
    with tab5:
        st.subheader("📋 Command Execution Log")
        st.write("View detailed logs of all commands executed in this session for debugging and verification.")
        
        # Get the command log from emboss wrapper
        log_entries = emboss.get_command_log()
        
        if not log_entries:
            st.info("No commands executed yet. Run some queries to see the execution log here.")
        else:
            # Summary stats
            total_commands = len(log_entries)
            successful_commands = sum(1 for entry in log_entries if entry['success'])
            failed_commands = total_commands - successful_commands
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Commands", total_commands)
            with col2:
                st.metric("Successful", successful_commands)
            with col3:
                st.metric("Failed", failed_commands)
            
            st.markdown("---")
            
            # Button to clear log
            if st.button("🗑️ Clear Log"):
                emboss.clear_log()
                st.rerun()
            
            # Download button for log
            log_text = emboss.get_formatted_log()
            st.download_button(
                label="📥 Download Log as Text",
                data=log_text,
                file_name=f"bioquery_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                mime="text/plain"
            )
            
            st.markdown("---")
            
            # Display each log entry
            for i, entry in enumerate(reversed(log_entries), 1):  # Show newest first
                entry_num = len(log_entries) - i + 1
                
                # Color-code by success/failure
                if entry['success']:
                    status_emoji = "✅"
                    border_color = "#4caf50"
                else:
                    status_emoji = "❌"
                    border_color = "#f44336"
                
                with st.expander(f"{status_emoji} Command #{entry_num}: {entry['tool']} - {entry['timestamp']}", expanded=(i==1)):
                    st.markdown(f"**Tool:** `{entry['tool']}`")
                    st.markdown(f"**Status:** {status_emoji} {'SUCCESS' if entry['success'] else 'FAILED'}")
                    st.markdown(f"**Timestamp:** {entry['timestamp']}")
                    
                    st.markdown("**Parameters:**")
                    param_lines = []
                    for key, value in entry['parameters'].items():
                        param_lines.append(f"- `{key}`: {value}")
                    st.markdown("\n".join(param_lines))
                    
                    if entry['error']:
                        st.error(f"**Error:** {entry['error']}")
                    else:
                        st.markdown("**Result Preview:**")
                        st.code(entry['result_preview'], language='text')
    
    # TAB 6: Documentation
    with tab6:
        st.subheader("How to use BioQuery Local")
        
        st.markdown("""
        ## Quick Start
        
        1. **Ask a Question** (Tab 1)
           - Type what you want to do with your sequences
           - Example: "Translate ATGAAATTTCCC to protein"
           - The AI will automatically choose the right tool
        
        2. **Or Select a Tool** (Tab 2)
           - Choose a specific EMBOSS tool
           - Enter your sequence(s)
           - Click to run
        
        ## Available Tools
        
        - **Translate**: Convert DNA to protein sequence
        - **Reverse**: Get reverse complement of DNA
        - **ORF**: Find open reading frames
        - **Align**: Align two sequences
        - **Restriction**: Find restriction enzyme sites
        - **Info**: Get sequence statistics
        - **Six Frame**: Show all 6 translation frames
        
        ## Sequence Format
        
        Sequences can be:
        - Plain text (ATGAAATTT...)
        - FASTA format (>header\\nsequence)
        - Multiple lines (will be concatenated)
        
        ## Technology
        
        - **EMBOSS**: Industry-standard bioinformatics toolkit
        - **Gemma3**: Local AI model for natural language understanding
        - **Streamlit**: Interactive web interface
        
        ## Examples
        
        Try these queries:
        - "Translate DNA to protein: ATGAAATTTCCCGGGAAATTT"
        - "What's the reverse complement of GCTA?"
        - "Find ORFs in ATGAAATTTCCCGGGAAATTTAAAGGG"
        - "Align these sequences: ATGAAA and ATGCCC"
        """)
        
        st.markdown("---")
        st.caption("BioQuery NoLocal v1.0 - Built with EMBOSS, Google Gemini, and Streamlit")


if __name__ == "__main__":
    main()
