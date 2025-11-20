"""
BioQuery Local - Streamlit App
Interactive web interface for natural language bioinformatics analysis
"""

import streamlit as st
from emboss_wrapper import EMBOSSWrapper
from nlp_handler import NLPHandler
import traceback


# Page configuration
st.set_page_config(
    page_title="BioQuery Local",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

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
def initialize_tools():
    """Initialize EMBOSS and NLP tools (cached for performance)"""
    emboss = EMBOSSWrapper()
    nlp = NLPHandler()
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
    
    for line in fasta_content.split('\n'):
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_header is not None:
                sequences.append((current_header, ''.join(current_seq)))
            # Start new sequence
            current_header = line[1:]  # Remove '>'
            current_seq = []
        else:
            # Add to current sequence
            current_seq.append(line)
    
    # Don't forget the last sequence
    if current_header is not None:
        sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def main():
    """Main Streamlit application"""
    
    # Initialize session state for sequence storage
    if 'sequences' not in st.session_state:
        st.session_state.sequences = {}
    
    # Header
    col1, col2 = st.columns([1, 4])
    with col1:
        st.title("🧬")
    with col2:
        st.title("BioQuery NoLocal")
        st.caption("Cloud-powered natural language bioinformatics analysis tool")
    
    st.markdown("---")
    
    # Initialize tools
    with st.spinner("Initializing bioinformatics tools..."):
        try:
            emboss, nlp = initialize_tools()
            emboss_ready = emboss.check_emboss()
            nlp_ready = nlp.test_connection()
        except Exception as e:
            st.error(f"Failed to initialize tools: {str(e)}")
            return
    
    # Sidebar with status and info
    with st.sidebar:
        st.header("⚙️ System Status")
        
        # Status indicators
        col1, col2 = st.columns(2)
        with col1:
            if emboss_ready:
                st.success("✓ EMBOSS")
            else:
                st.error("✗ EMBOSS")
        with col2:
            if nlp_ready:
                st.success("✓ Gemini")
            else:
                st.error("✗ Gemini")
        
        st.markdown("---")
        
        # Available tools
        st.subheader("📊 Available Tools")
        tools = emboss.get_available_tools()
        for tool_name, description in tools.items():
            with st.expander(f"**{tool_name}**"):
                st.write(description)
        
        st.markdown("---")
        
        # About
        st.subheader("ℹ️ About")
        st.info(
            "BioQuery NoLocal uses:\n"
            "- **EMBOSS**: Bioinformatics analysis\n"
            "- **Google Gemini**: AI for natural language\n"
            "- **Streamlit**: Interactive interface"
        )
    
    # Main interface tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🤖 Natural Language Query",
        "🔧 Manual Tool Selection",
        "🌐 Genome Browser Query",
        "📁 FASTA Upload & Batch",
        "📚 Documentation"
    ])
    
    # TAB 1: Natural Language Query
    with tab1:
        st.subheader("Ask a bioinformatics question")
        st.write("Describe what you want to do with your sequences in natural language.")
        
        col1, col2 = st.columns([3, 1])
        with col1:
            user_query = st.text_area(
                "Your query:",
                placeholder="E.g., 'Translate this DNA sequence to protein: ATGAAATTTCCC' or 'What's the reverse complement of ATGAAA?'",
                height=100,
                key="nlp_query"
            )
        
        with col2:
            submit_button = st.button("🚀 Analyze", key="analyze_btn", use_container_width=True)
        
        if submit_button and user_query.strip():
            with st.spinner("Processing your query with AI..."):
                try:
                    # Parse the query
                    success, result = nlp.parse_user_query(user_query)
                    
                    if success:
                        # Check if this is multiple unrelated questions
                        if result.get('type') == 'multiple_questions':
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
                        st.markdown('<div class="error-box">', unsafe_allow_html=True)
                        st.error(f"Could not parse query: {result.get('error', 'Unknown error')}")
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
    
    # TAB 3: Genome Browser Query
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
        
        st.info("Upload a FASTA file to analyze multiple sequences at once. Select a tool to apply to all sequences.")
        
        # File upload
        uploaded_file = st.file_uploader(
            "Choose a FASTA file:",
            type=["fasta", "fa", "faa", "fna", "txt"],
            help="Upload a file containing sequences in FASTA format"
        )
        
        if uploaded_file:
            st.success(f"✓ File uploaded: {uploaded_file.name}")
            
            # Read and parse file
            try:
                content = uploaded_file.read().decode('utf-8')
                sequences = parse_fasta(content)
                
                if not sequences:
                    st.warning("No sequences found in file. Make sure it's in FASTA format (headers start with '>').")
                else:
                    st.info(f"Found {len(sequences)} sequence(s) in file")
                    
                    # Show preview
                    with st.expander("📄 Preview sequences"):
                        for i, (header, seq) in enumerate(sequences[:5], 1):  # Show first 5
                            st.text(f">{header}")
                            st.text(f"{seq[:100]}{'...' if len(seq) > 100 else ''}")
                            st.caption(f"Length: {len(seq)} bp/aa")
                            st.markdown("---")
                        if len(sequences) > 5:
                            st.caption(f"... and {len(sequences) - 5} more sequences")
                    
                    st.markdown("---")
                    
                    # Tool selection for batch processing
                    st.subheader("Select tool to apply")
                    
                    col1, col2 = st.columns([2, 1])
                    
                    with col1:
                        batch_tool = st.selectbox(
                            "Choose analysis tool:",
                            ["translate", "reverse", "gc", "info", "orf", "sixframe", "iep", "pepstats"],
                            key="batch_tool"
                        )
                    
                    with col2:
                        # Tool-specific parameters
                        params = {}
                        if batch_tool == "translate":
                            params['frame'] = st.slider("Reading frame:", 1, 3, 1, key="batch_frame")
                        elif batch_tool == "orf":
                            params['min_size'] = st.slider("Min ORF size:", 10, 500, 100, key="batch_orf")
                    
                    # Run analysis button
                    if st.button("🚀 Analyze All Sequences", type="primary", key="batch_run"):
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        results_container = st.container()
                        
                        all_results = []
                        
                        with results_container:
                            for i, (header, seq) in enumerate(sequences):
                                progress = (i + 1) / len(sequences)
                                progress_bar.progress(progress)
                                status_text.text(f"Processing {i+1}/{len(sequences)}: {header[:50]}...")
                                
                                try:
                                    # Run the selected tool
                                    if batch_tool == "translate":
                                        result = emboss.translate_sequence(seq, params.get('frame', 1))
                                    elif batch_tool == "reverse":
                                        result = emboss.reverse_complement(seq)
                                    elif batch_tool == "gc":
                                        result = emboss.calculate_gc_content(seq)
                                    elif batch_tool == "info":
                                        result = emboss.get_sequence_info(seq)
                                    elif batch_tool == "orf":
                                        result = emboss.find_orfs(seq, params.get('min_size', 100))
                                    elif batch_tool == "sixframe":
                                        result = emboss.get_six_frame_translation(seq)
                                    elif batch_tool == "iep":
                                        result = emboss.calculate_isoelectric_point(seq)
                                    elif batch_tool == "pepstats":
                                        result = emboss.get_protein_stats(seq)
                                    else:
                                        result = f"Tool {batch_tool} not yet implemented for batch"
                                    
                                    all_results.append((header, result))
                                    
                                except Exception as e:
                                    all_results.append((header, f"ERROR: {str(e)}"))
                            
                            status_text.text("✓ Analysis complete!")
                            progress_bar.progress(1.0)
                        
                        # Display results
                        st.markdown("---")
                        st.subheader("📊 Results")
                        
                        for header, result in all_results:
                            with st.expander(f"📄 {header}"):
                                st.code(result, language="text")
                        
                        # Download all results
                        combined_results = ""
                        for header, result in all_results:
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
            - **translate**: Convert DNA to protein
            - **reverse**: Reverse complement of DNA
            - **gc**: Calculate GC content percentage  
            - **info**: Sequence statistics (length, type, etc.)
            - **orf**: Find open reading frames
            - **sixframe**: Show all 6 translation frames
            - **iep**: Calculate isoelectric point (proteins)
            - **pepstats**: Protein statistics (molecular weight, amino acid composition)
            """)
    
    # TAB 5: Documentation
    with tab5:
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
