"""
NLP Handler for BioQuery Local
Converts natural language requests to EMBOSS tool calls using Ollama
"""

import json
from typing import Dict, Optional, Tuple
from ollama import Client
import re


class NLPHandler:
    """Convert natural language queries to EMBOSS tool commands using Ollama"""
    
    def __init__(self, ollama_host: str = "http://192.168.128.1:11434", model: str = "gemma3:4b"):
        """Initialize the NLP handler with Ollama client
        
        Args:
            ollama_host: URL to Ollama server
            model: Model name to use (default: gemma3:4b)
        """
        self.ollama_host = ollama_host
        self.model = model
        self.client = Client(host=ollama_host)
        
        # Define the system prompt for the LLM
        self.system_prompt = """You are a bioinformatics assistant that helps users run EMBOSS analysis tools and query genomic databases.

When a user asks a question about DNA/protein sequences, genomic regions, genes, or tools, respond with a JSON object containing:

FOR SINGLE OPERATIONS:
1. "tool": the tool name (EMBOSS tool, 'genome_query', or 'gene_query')
2. "parameters": a dict with required parameters
3. "explanation": a brief explanation of what will be done

FOR MULTI-STEP OPERATIONS (when user says "then", "and then", "after that", "next", etc.):
1. "steps": a list of step objects, each containing "tool" and "parameters"
2. "explanation": overall explanation of the workflow
3. "use_previous_result": whether step N+1 should use output from step N

IMPORTANT MULTI-STEP CHAINING RULES:
- If Step 1 is gene_query and Step 2 is BLAST: DO NOT include "sequence" parameter in Step 2. The system will automatically fetch the sequence from the gene.
- If Step 1 is gene_query and Step 2 is translate/gc/analysis: DO NOT include "sequence" parameter in Step 2. The system will use the gene's sequence.
- Only provide "sequence" parameter if user explicitly gives you a raw sequence.
- For gene_name chaining: Just pass the gene_name in Step 1, and the gene_name (not sequence) in dependent steps - the system will resolve it.

IMPORTANT: When a user mentions a GENE SYMBOL (like ALKBH1, TP53, BRCA1), use the gene_name parameter, NOT sequence.

EMBOSS tools (use 'gene_name' for gene symbols, 'sequence' for raw DNA/protein sequences):
- translate: Translate DNA to protein. Needs "gene_name" OR "sequence", optional "frame" (1-3) and "transcript_variant" (e.g., "transcript variant 5")
- reverse: Reverse complement DNA. Needs "gene_name" OR "sequence"
- orf: Find open reading frames. Needs "gene_name" OR "sequence" and optional "min_size"
- align: Align two sequences. Needs "seq1" and "seq2"
- pattern: Search for patterns. Needs "gene_name" OR "sequence" and optional "pattern"
- restriction: Find restriction sites. Needs "gene_name" OR "sequence" and optional "enzyme"
- shuffle: Shuffle sequence. Needs "gene_name" OR "sequence"
- info: Get sequence info. Needs "gene_name" OR "sequence"
- sixframe: Show all 6 reading frames. Needs "gene_name" OR "sequence"
- gc: Calculate GC content. Needs "gene_name" OR "sequence"
- blast: Search NCBI databases for similar sequences. Needs "sequence" and optional "blast_type" (blastn/blastp/blastx), "database" (nt/nr), "max_results" (default 10), "expect_threshold" (default 10.0)
- blastn: DNA BLAST search. Needs "sequence"
- blastp: Protein BLAST search. Needs "sequence"
- blastx: DNA to protein BLAST search. Needs "sequence"
- search: General sequence search (uses blastn by default). Needs "sequence"

GENOME QUERY tool:
- genome_query: Query UCSC Genome Browser for genomic regions. Needs "genome", "chrom", "start", "end"

GENE QUERY tool:
- gene_query: Look up gene information (exons, CDS, transcript length). Needs "gene_name" and optional "genome" (default hg38) and "track" (default gencode)

Decision logic:
- If user mentions a GENE NAME/SYMBOL (like ALKBH1, TP53, BRCA1, etc.):
  - If asking about gene structure (exons, CDS, transcripts) -> use gene_query with gene_name
  - If asking to apply a tool (translate, gc, BLAST, etc.) to that gene -> use the tool with gene_name parameter and transcript_variant if specified
- If user provides raw DNA/RNA/protein sequences (like ATGCCC or MKLA...) or mentions "BLAST", "search", "similar sequences", "homologous" -> use the appropriate search/BLAST tool with sequence parameter
- If user mentions chromosome/genomic position -> use genome_query
- Otherwise use appropriate EMBOSS tool

MULTI-STEP EXAMPLE:
User: "Find gene info for ALKBH1, then calculate its GC content"
Response:
{
  "steps": [
    {"tool": "gene_query", "parameters": {"gene_name": "ALKBH1"}},
    {"tool": "gc", "parameters": {"gene_name": "ALKBH1"}}
  ],
  "explanation": "First query gene information for ALKBH1, then calculate GC content of its sequence",
  "use_previous_result": false
}

TRANSCRIPT VARIANT EXTRACTION:
- If user mentions "transcript variant N" or "primary variant" or "variant N", extract this as "transcript_variant": "transcript variant N"
- Examples: "transcript variant 5" -> "transcript_variant": "transcript variant 5"

Examples:
- "How many exons in ALKBH1?" -> gene_query with gene_name: ALKBH1
- "Find the gc content of ALKBH1" -> gc tool with gene_name: ALKBH1
- "Translate ATGCCC" -> translate tool with sequence: ATGCCC
- "Translate TP53" -> translate tool with gene_name: TP53
- "What's the reverse complement of BRCA1?" -> reverse tool with gene_name: BRCA1
- "Give the length of the protein made from transcript variant 5 of CARS1" -> translate tool with gene_name: CARS1, transcript_variant: "transcript variant 5"
- "Find gene info for ALKBH1, then calculate its GC content" -> multi-step with gene_query followed by gc
- "BLAST this sequence against the human genome: ATGCGATCG" -> blast tool with sequence: ATGCGATCG, blast_type: blastn
- "Find BRCA1 sequence, then BLAST it to find homologous genes" -> multi-step: gene_query for BRCA1 (NO sequence param in step 2), then blast (NO sequence param - system auto-fetches)
- "Search for similar proteins to MKLASELKD" -> blastp tool with sequence: MKLASELKD
- "Find homologous DNA sequences to this: ATGCCC" -> blastn tool with sequence: ATGCCC
- "Find TP53 gene info then translate it" -> multi-step: gene_query(gene_name: TP53), translate(gene_name: TP53, NOT sequence)

Always respond with ONLY valid JSON, no other text. Start with { and end with }"""

    def parse_user_query(self, query: str) -> Tuple[bool, Dict]:
        """Parse a user query and extract tool information
        
        Args:
            query: User's natural language query
        
        Returns:
            Tuple of (success: bool, result: dict with tool, parameters, explanation OR steps, explanation)
        """
        try:
            # Call Ollama with the system prompt
            response = self.client.generate(
                model=self.model,
                prompt=query,
                system=self.system_prompt,
                stream=False
            )
            
            # Extract the response text
            response_text = response['response'].strip()
            
            # Remove markdown code blocks if present
            if response_text.startswith('```'):
                response_text = response_text.split('```')[1]
                if response_text.startswith('json'):
                    response_text = response_text[4:]
            if response_text.endswith('```'):
                response_text = response_text[:-3]
            
            response_text = response_text.strip()
            
            # Parse JSON response
            result = json.loads(response_text)
            
            # Validate the response - can be either single tool or multi-step
            if 'steps' in result:
                # Multi-step query - clean up parameters
                if not isinstance(result['steps'], list) or len(result['steps']) == 0:
                    return False, {"error": "Invalid steps format"}
                
                # Smart parameter cleanup for multi-step workflows
                result['steps'] = self._cleanup_multistep_parameters(result['steps'])
                return True, result
            elif 'tool' in result and 'parameters' in result:
                # Single-step query
                return True, result
            else:
                return False, {"error": "Invalid response format from LLM"}
        
        except json.JSONDecodeError as e:
            return False, {"error": f"Failed to parse LLM response as JSON: {str(e)}"}
        except Exception as e:
            return False, {"error": f"Error processing query: {str(e)}"}
    
    def _cleanup_multistep_parameters(self, steps: list) -> list:
        """Clean up parameters in multi-step queries to enable proper chaining
        
        Rules:
        - If previous step is gene_query and current step is BLAST/translate/gc, remove 'sequence' param
        - This allows the system to automatically fetch sequences and chain them
        
        Args:
            steps: List of step dicts
        
        Returns:
            List of cleaned step dicts
        """
        cleaned_steps = []
        
        for idx, step in enumerate(steps):
            step_copy = step.copy()
            tool_name = step_copy.get('tool', '').lower()
            parameters = step_copy.get('parameters', {})
            if not isinstance(parameters, dict):
                parameters = {}
            parameters = parameters.copy()  # Make a copy to avoid modifying original
            
            # Check if previous step was gene_query
            if idx > 0:
                prev_step = steps[idx - 1]
                prev_tool = prev_step.get('tool', '').lower()
                
                # If previous was gene_query and current is a dependent tool
                if prev_tool == 'gene_query' and tool_name in ['blast', 'blastn', 'blastp', 'blastx', 'search', 'translate', 'gc', 'gc_content', 'reverse', 'orf', 'getorf']:
                    # Remove 'sequence' if it's just a gene name (not actual bases)
                    if 'sequence' in parameters:
                        seq_value = parameters['sequence']
                        # If sequence is short and all letters (like 'ALKBH1'), it's a gene name - remove it
                        if isinstance(seq_value, str) and len(seq_value) < 20 and seq_value.replace('_', '').isalpha():
                            print(f"🔧 [Cleanup] Removing sequence='{seq_value}' from step {idx+1} ({tool_name}) - detected as gene name")
                            del parameters['sequence']
            
            step_copy['parameters'] = parameters
            cleaned_steps.append(step_copy)
        
        return cleaned_steps
    
    def suggest_tools(self, query: str) -> list:
        """Suggest relevant EMBOSS tools based on a query
        
        Args:
            query: User's question
        
        Returns:
            List of suggested tool names
        """
        query_lower = query.lower()
        suggestions = []
        
        # Simple keyword matching for suggestions
        if any(word in query_lower for word in ['translate', 'protein', 'amino']):
            suggestions.append('translate')
        
        if any(word in query_lower for word in ['reverse', 'complement', 'rc']):
            suggestions.append('reverse')
        
        if any(word in query_lower for word in ['orf', 'open reading frame', 'gene', 'start codon']):
            suggestions.append('orf')
        
        if any(word in query_lower for word in ['align', 'alignment', 'compare', 'similarity']):
            suggestions.append('align')
        
        if any(word in query_lower for word in ['restrict', 'restriction', 'enzyme', 'cut']):
            suggestions.append('restriction')
        
        if any(word in query_lower for word in ['info', 'information', 'length', 'gc content']):
            suggestions.append('info')
        
        if any(word in query_lower for word in ['pattern', 'motif', 'search', 'find']):
            suggestions.append('pattern')
        
        if any(word in query_lower for word in ['six frame', 'sixframe', 'all frames']):
            suggestions.append('sixframe')
        
        if any(word in query_lower for word in ['shuffle', 'randomize']):
            suggestions.append('shuffle')
        
        return suggestions
    
    def validate_sequence(self, sequence: str) -> Tuple[bool, str]:
        """Validate if a sequence is in correct format
        
        Args:
            sequence: DNA or protein sequence
        
        Returns:
            Tuple of (is_valid: bool, message: str)
        """
        # Remove whitespace and convert to uppercase
        seq_clean = sequence.replace('\n', '').replace(' ', '').replace('\t', '').upper()
        
        if not seq_clean:
            return False, "Sequence is empty"
        
        # Check for valid DNA characters
        valid_dna = set('ATCGN')
        valid_protein = set('ACDEFGHIKLMNPQRSTVWY')
        
        seq_set = set(seq_clean)
        
        # Check if it looks like DNA or protein
        if seq_set.issubset(valid_dna):
            return True, f"Valid DNA sequence ({len(seq_clean)} bp)"
        elif seq_set.issubset(valid_protein):
            return True, f"Valid protein sequence ({len(seq_clean)} aa)"
        else:
            invalid_chars = seq_set - (valid_dna | valid_protein)
            return False, f"Invalid characters in sequence: {invalid_chars}"
    
    def extract_sequence_from_query(self, query: str) -> Optional[str]:
        """Try to extract a sequence from a user query
        
        Args:
            query: User query that might contain a sequence
        
        Returns:
            Extracted sequence or None
        """
        # Look for FASTA format
        fasta_match = re.search(r'>.*?\n([ATCGN]+)', query, re.IGNORECASE)
        if fasta_match:
            return fasta_match.group(1)
        
        # Look for plain sequence (continuous letters)
        seq_match = re.search(r'\b([ATCGN]{15,})\b', query, re.IGNORECASE)
        if seq_match:
            return seq_match.group(1)
        
        return None
    
    def format_response(self, success: bool, tool_info: Dict, result: str = None) -> str:
        """Format a nicely readable response
        
        Args:
            success: Whether the operation was successful
            tool_info: Information about the tool used
            result: The result from running the tool
        
        Returns:
            Formatted response string
        """
        if not success:
            return f"❌ Error: {tool_info.get('error', 'Unknown error')}"
        
        response = f"✓ Running {tool_info.get('tool', 'tool')}\n"
        response += f"📋 Explanation: {tool_info.get('explanation', 'N/A')}\n"
        
        if tool_info.get('parameters'):
            response += f"⚙️  Parameters:\n"
            for param, value in tool_info['parameters'].items():
                if isinstance(value, str) and len(value) > 50:
                    response += f"   - {param}: {value[:50]}...\n"
                else:
                    response += f"   - {param}: {value}\n"
        
        if result:
            response += f"\n📊 Results:\n{result}"
        
        return response
    
    def test_connection(self) -> bool:
        """Test if Ollama connection is working
        
        Returns:
            bool: True if connection successful
        """
        try:
            models = self.client.list()
            if models and len(models.models) > 0:
                print(f"✓ Connected to Ollama at {self.ollama_host}")
                print(f"✓ Available models: {[m.model for m in models.models]}")
                return True
            else:
                print("✗ No models available in Ollama")
                return False
        except Exception as e:
            print(f"✗ Failed to connect to Ollama: {str(e)}")
            return False


# Example usage
if __name__ == "__main__":
    print("=== BioQuery NLP Handler ===\n")
    
    # Initialize handler
    handler = NLPHandler()
    
    # Test connection
    print("Testing Ollama connection...")
    if handler.test_connection():
        print("\n" + "="*50 + "\n")
        
        # Test queries
        test_queries = [
            "Translate this DNA sequence to protein: ATGAAATTTCCCGGGAAATTT",
            "What is the reverse complement of ATGAAATTT?",
            "Find all open reading frames in ATGAAATTTCCCGGGAAATTTAAAGGG",
        ]
        
        for query in test_queries:
            print(f"Query: {query}")
            success, result = handler.parse_user_query(query)
            if success:
                print(f"✓ Tool: {result.get('tool')}")
                print(f"✓ Parameters: {result.get('parameters')}")
                print(f"✓ Explanation: {result.get('explanation')}")
            else:
                print(f"✗ Error: {result.get('error')}")
            print()
    else:
        print("Cannot proceed without Ollama connection")
