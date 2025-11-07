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


def main():
    """Main Streamlit application"""
    
    # Header
    col1, col2 = st.columns([1, 4])
    with col1:
        st.title("🧬")
    with col2:
        st.title("BioQuery Local")
        st.caption("Self-contained natural language bioinformatics analysis tool")
    
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
                st.success("✓ Ollama")
            else:
                st.error("✗ Ollama")
        
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
            "BioQuery Local uses:\n"
            "- **EMBOSS**: Bioinformatics analysis\n"
            "- **Ollama**: Local LLM for natural language\n"
            "- **Streamlit**: Interactive interface"
        )
    
    # Main interface tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "🤖 Natural Language Query",
        "🔧 Manual Tool Selection",
        "📊 Sequence Analysis",
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
                                analysis_result = emboss.run_tool(tool_name, **parameters)
                                
                                st.markdown('<div class="success-box">', unsafe_allow_html=True)
                                st.success("Analysis completed!")
                                st.subheader("Results:")
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
    
    # TAB 3: Batch Sequence Analysis
    with tab3:
        st.subheader("Analyze multiple sequences")
        
        st.info("Upload a FASTA file to analyze multiple sequences at once")
        
        uploaded_file = st.file_uploader("Choose a FASTA file:", type=["fasta", "fa", "faa", "txt"])
        
        if uploaded_file:
            st.success(f"File uploaded: {uploaded_file.name}")
            
            # Read file
            content = uploaded_file.read().decode('utf-8')
            st.text_area("File content:", content, height=200)
    
    # TAB 4: Documentation
    with tab4:
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
        st.caption("BioQuery Local v1.0 - Built with EMBOSS, Ollama, and Streamlit")


if __name__ == "__main__":
    main()
