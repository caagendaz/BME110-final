"""
EMBOSS Wrapper for BioQuery Local
Maps natural language to EMBOSS bioinformatics tools
"""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class EMBOSSWrapper:
    """Wrapper for common EMBOSS bioinformatics tools"""
    
    def __init__(self):
        """Initialize the EMBOSS wrapper and verify installation"""
        # Map natural language to EMBOSS commands
        self.tool_map = {
            'translate': 'transeq',
            'reverse': 'revseq',
            'orf': 'getorf',
            'align': 'needle',
            'pattern': 'fuzznuc',
            'restriction': 'restrict',
            'shuffle': 'shuffleseq',
            'info': 'infoseq',
            'sixframe': 'sixpack',
            'palindrome': 'palindrome',
            'fuzzy': 'fuzzpro',
            'dotplot': 'dotmatcher',
            'consensus': 'cons'
        }
        
        # Tool descriptions for user guidance
        self.tool_descriptions = {
            'translate': 'Translate a DNA sequence into protein',
            'reverse': 'Reverse and complement a DNA sequence',
            'orf': 'Find open reading frames in a sequence',
            'align': 'Global sequence alignment using Needleman-Wunsch',
            'pattern': 'Search for patterns in sequences',
            'restriction': 'Find restriction enzyme sites',
            'shuffle': 'Shuffle a sequence',
            'info': 'Get information about sequences',
            'sixframe': 'Display all six translation frames',
            'palindrome': 'Find palindromic sequences',
            'fuzzy': 'Fuzzy pattern matching',
            'dotplot': 'Create dot plots between sequences',
            'consensus': 'Generate consensus sequence'
        }
        
        # Verify EMBOSS installation
        self.check_emboss()
    
    def check_emboss(self) -> bool:
        """Verify EMBOSS tools are available
        
        Returns:
            bool: True if EMBOSS is installed, False otherwise
        """
        try:
            result = subprocess.run(
                ['embossversion'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                version_info = result.stdout.strip()
                print(f"✓ EMBOSS found: {version_info}")
                return True
            else:
                print("✗ EMBOSS not found or error occurred")
                return False
        except FileNotFoundError:
            print("✗ EMBOSS not found. Please install with: conda install -c bioconda emboss")
            return False
        except subprocess.TimeoutExpired:
            print("✗ EMBOSS check timed out")
            return False
    
    def get_available_tools(self) -> Dict[str, str]:
        """Get list of available tools and their descriptions
        
        Returns:
            dict: Mapping of tool names to descriptions
        """
        return self.tool_descriptions.copy()
    
    def translate_sequence(self, sequence: str, frame: int = 1) -> str:
        """Translate DNA sequence to protein
        
        Args:
            sequence: DNA sequence string
            frame: Reading frame (1, 2, or 3)
        
        Returns:
            str: Translated protein sequence
        """
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">input_seq\n{sequence}\n")
                input_file = f.name
            
            output_file = input_file.replace('.fasta', '_translated.fasta')
            
            result = subprocess.run(
                ['transeq', '-sequence', input_file, '-outseq', output_file, '-frame', str(frame)],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0 and os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    translation = f.read()
                # Cleanup
                os.remove(input_file)
                os.remove(output_file)
                return translation
            else:
                return f"Error: {result.stderr}"
        
        except Exception as e:
            return f"Error during translation: {str(e)}"
    
    def reverse_complement(self, sequence: str) -> str:
        """Reverse complement a DNA sequence
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            str: Reverse complemented sequence
        """
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">input_seq\n{sequence}\n")
                input_file = f.name
            
            output_file = input_file.replace('.fasta', '_rc.fasta')
            
            result = subprocess.run(
                ['revseq', '-sequence', input_file, '-outseq', output_file],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0 and os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    rc_seq = f.read()
                # Cleanup
                os.remove(input_file)
                os.remove(output_file)
                return rc_seq
            else:
                return f"Error: {result.stderr}"
        
        except Exception as e:
            return f"Error during reverse complement: {str(e)}"
    
    def find_orfs(self, sequence: str, min_size: int = 100) -> str:
        """Find open reading frames in a sequence
        
        Args:
            sequence: DNA sequence string
            min_size: Minimum ORF size in base pairs
        
        Returns:
            str: ORFs found in the sequence
        """
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">input_seq\n{sequence}\n")
                input_file = f.name
            
            output_file = input_file.replace('.fasta', '_orfs.fasta')
            
            result = subprocess.run(
                ['getorf', '-sequence', input_file, '-outseq', output_file, '-minsize', str(min_size)],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0 and os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    orfs = f.read()
                # Cleanup
                os.remove(input_file)
                os.remove(output_file)
                return orfs
            else:
                return f"Error: {result.stderr}"
        
        except Exception as e:
            return f"Error finding ORFs: {str(e)}"
    
    def get_sequence_info(self, sequence: str) -> str:
        """Get information about a sequence
        
        Args:
            sequence: DNA or protein sequence string
        
        Returns:
            str: Sequence information (length, GC content, etc.)
        """
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">input_seq\n{sequence}\n")
                input_file = f.name
            
            result = subprocess.run(
                ['infoseq', '-sequence', input_file],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            # Cleanup
            os.remove(input_file)
            
            if result.returncode == 0:
                return result.stdout
            else:
                return f"Error: {result.stderr}"
        
        except Exception as e:
            return f"Error getting sequence info: {str(e)}"
    
    def find_restriction_sites(self, sequence: str, enzyme: str = None) -> str:
        """Find restriction enzyme sites in a sequence
        
        Args:
            sequence: DNA sequence string
            enzyme: Specific enzyme to search for (optional)
        
        Returns:
            str: Restriction sites found
        """
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">input_seq\n{sequence}\n")
                input_file = f.name
            
            output_file = input_file.replace('.fasta', '_restrict.txt')
            
            cmd = ['restrict', '-sequence', input_file, '-outfile', output_file]
            if enzyme:
                cmd.extend(['-enzyme', enzyme])
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0 and os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    sites = f.read()
                # Cleanup
                os.remove(input_file)
                os.remove(output_file)
                return sites
            else:
                return f"Error: {result.stderr}"
        
        except Exception as e:
            return f"Error finding restriction sites: {str(e)}"
    
    def align_sequences(self, seq1: str, seq2: str) -> str:
        """Align two sequences using Needleman-Wunsch algorithm
        
        Args:
            seq1: First sequence
            seq2: Second sequence
        
        Returns:
            str: Alignment result
        """
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f1:
                f1.write(f">seq1\n{seq1}\n")
                file1 = f1.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f2:
                f2.write(f">seq2\n{seq2}\n")
                file2 = f2.name
            
            output_file = tempfile.NamedTemporaryFile(suffix='.txt', delete=False).name
            
            result = subprocess.run(
                ['needle', '-asequence', file1, '-bsequence', file2, '-outfile', output_file],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0 and os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    alignment = f.read()
                # Cleanup
                os.remove(file1)
                os.remove(file2)
                os.remove(output_file)
                return alignment
            else:
                return f"Error: {result.stderr}"
        
        except Exception as e:
            return f"Error during alignment: {str(e)}"
    
    def get_six_frame_translation(self, sequence: str) -> str:
        """Get all six reading frame translations
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            str: Six frame translation output
        """
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">input_seq\n{sequence}\n")
                input_file = f.name
            
            output_file = input_file.replace('.fasta', '_sixpack.fasta')
            
            result = subprocess.run(
                ['sixpack', '-sequence', input_file, '-outfile', output_file],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0 and os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    translation = f.read()
                # Cleanup
                os.remove(input_file)
                os.remove(output_file)
                return translation
            else:
                return f"Error: {result.stderr}"
        
        except Exception as e:
            return f"Error during six frame translation: {str(e)}"
    
    def run_tool(self, tool_name: str, **kwargs) -> str:
        """Generic method to run any EMBOSS tool
        
        Args:
            tool_name: Name of the tool to run (natural language or EMBOSS name)
            **kwargs: Tool-specific parameters
        
        Returns:
            str: Tool output or error message
        """
        # Map natural language name to EMBOSS command if needed
        emboss_name = self.tool_map.get(tool_name.lower(), tool_name.lower())
        
        # Route to specific methods
        if emboss_name == 'transeq':
            return self.translate_sequence(kwargs.get('sequence', ''), kwargs.get('frame', 1))
        elif emboss_name == 'revseq':
            return self.reverse_complement(kwargs.get('sequence', ''))
        elif emboss_name == 'getorf':
            return self.find_orfs(kwargs.get('sequence', ''), kwargs.get('min_size', 100))
        elif emboss_name == 'infoseq':
            return self.get_sequence_info(kwargs.get('sequence', ''))
        elif emboss_name == 'restrict':
            return self.find_restriction_sites(kwargs.get('sequence', ''), kwargs.get('enzyme', None))
        elif emboss_name == 'needle':
            return self.align_sequences(kwargs.get('seq1', ''), kwargs.get('seq2', ''))
        elif emboss_name == 'sixpack':
            return self.get_six_frame_translation(kwargs.get('sequence', ''))
        else:
            return f"Tool '{tool_name}' not yet implemented"


if __name__ == "__main__":
    # Test the wrapper
    wrapper = EMBOSSWrapper()
    
    print("\n=== Available EMBOSS Tools ===")
    for tool, desc in wrapper.get_available_tools().items():
        print(f"  {tool}: {desc}")
    
    # Test with a sample DNA sequence
    test_seq = "ATGAAATTTCCCGGGAAATTT"
    print(f"\n=== Testing with sequence: {test_seq} ===")
    
    print("\nSequence Info:")
    print(wrapper.get_sequence_info(test_seq))
    
    print("\nReverse Complement:")
    print(wrapper.reverse_complement(test_seq))
