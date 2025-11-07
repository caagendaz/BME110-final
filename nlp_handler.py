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
        self.system_prompt = """You are a bioinformatics assistant that helps users run EMBOSS analysis tools.

When a user asks a question about DNA/protein sequences, respond with a JSON object containing:
1. "tool": the EMBOSS tool name (translate, reverse, orf, align, pattern, restriction, shuffle, info, sixframe)
2. "parameters": a dict with required parameters
3. "explanation": a brief explanation of what will be done

Available tools:
- translate: Translate DNA to protein. Needs "sequence" and optional "frame" (1-3)
- reverse: Reverse complement DNA. Needs "sequence"
- orf: Find open reading frames. Needs "sequence" and optional "min_size"
- align: Align two sequences. Needs "seq1" and "seq2"
- pattern: Search for patterns. Needs "sequence" and optional "pattern"
- restriction: Find restriction sites. Needs "sequence" and optional "enzyme"
- shuffle: Shuffle sequence. Needs "sequence"
- info: Get sequence info. Needs "sequence"
- sixframe: Show all 6 reading frames. Needs "sequence"

Always respond with ONLY valid JSON, no other text. Start with { and end with }"""

    def parse_user_query(self, query: str) -> Tuple[bool, Dict]:
        """Parse a user query and extract tool information
        
        Args:
            query: User's natural language query
        
        Returns:
            Tuple of (success: bool, result: dict with tool, parameters, explanation)
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
            
            # Validate the response
            if 'tool' not in result or 'parameters' not in result:
                return False, {"error": "Invalid response format from LLM"}
            
            return True, result
        
        except json.JSONDecodeError as e:
            return False, {"error": f"Failed to parse LLM response as JSON: {str(e)}"}
        except Exception as e:
            return False, {"error": f"Error processing query: {str(e)}"}
    
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
