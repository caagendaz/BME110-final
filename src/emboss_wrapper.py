"""
EMBOSS Wrapper for BioQuery NoLocal
Maps natural language to EMBOSS bioinformatics tools
"""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio import Entrez
import requests
import json
from datetime import datetime
import traceback


class EMBOSSWrapper:
    """Wrapper for common EMBOSS bioinformatics tools"""
    
    def __init__(self):
        """Initialize the EMBOSS wrapper and verify installation"""
        # Set NCBI Entrez email for downloads
        Entrez.email = "user@bioquery.local"
        
        # Initialize command log
        self.command_log = []
        
        # Map natural language to EMBOSS commands (shortcuts for common operations)
        self.tool_map = {
            'translate': 'transeq',
            'reverse': 'revseq',
            'orf': 'getorf',
            'align': 'needle',
            'water': 'water',
            'waterman': 'water',
            'local_align': 'water',
            'smith_waterman': 'water',
            'pattern': 'fuzznuc',
            'restriction': 'restrict',
            'shuffle': 'shuffleseq',
            'info': 'infoseq',
            'sixframe': 'sixpack',
            'palindrome': 'palindrome',
            'fuzzy': 'fuzzpro',
            'dotplot': 'dotmatcher',
            'consensus': 'cons',
            'gc': 'gc_content',
            'download': 'download_sequence',
            'geecee': 'gc_content',
            'isoelectric': 'iep',
            'charge': 'iep',
            'gene_query': 'gene_query',
            'gene': 'gene_query',
            'blast': 'blast',
            'blastn': 'blastn',
            'blastp': 'blastp',
            'blastx': 'blastx',
            'search': 'blast',
            'cusp': 'cusp',
            'codon_usage': 'cusp',
            'pepstats': 'pepstats',
            'protein_stats': 'pepstats',
            'molecular_weight': 'pepstats',
            'weight': 'pepstats',
            'mass': 'pepstats',
            'wordcount': 'wordcount',
            'oligonucleotide': 'wordcount',
            'bedtools': 'bedtools_intersect',
            'intersect': 'bedtools_intersect',
            'blat': 'blat_search',
            'ucsc_gene': 'ucsc_gene_info',
            'gtex': 'gtex_expression',
            'expression': 'gtex_expression',
            'tissue_expression': 'gtex_expression',
            'ucsc_table': 'ucsc_table_query',
            'table_browser': 'ucsc_table_query',
            'pubmed': 'pubmed_search',
            'literature': 'pubmed_search',
            'articles': 'pubmed_search',
            'overlap': 'track_intersection',
            'intersect_tracks': 'track_intersection',
            'track_overlap': 'track_intersection',
            'neighboring_genes': 'find_neighboring_genes',
            'flanking_genes': 'find_neighboring_genes',
            'genes_near': 'find_neighboring_genes',
            'neighbor_genes': 'find_neighboring_genes',
            # Biopython-powered tools
            'msa': 'multiple_sequence_alignment',
            'multiple_align': 'multiple_sequence_alignment',
            'pdb': 'protein_structure',
            'structure': 'protein_structure',
            'protein_structure': 'protein_structure',
            'biopython_mw': 'molecular_weight_biopython',
            'protparam': 'molecular_weight_biopython',
            # Phylogenetics
            'phylo': 'phylogenetic_tree',
            'phylogenetic_tree': 'phylogenetic_tree',
            'tree': 'phylogenetic_tree',
            'phylogeny': 'phylogenetic_tree',
            # Motif analysis
            'motif': 'find_motifs',
            'motifs': 'find_motifs',
            'find_motif': 'find_motifs',
            'consensus_motif': 'motif_consensus',
            # Primer design
            'primer': 'primer_analysis',
            'primers': 'primer_analysis',
            'primer_tm': 'primer_analysis',
            'melting_temp': 'primer_analysis',
            # Restriction enzymes (advanced)
            'restriction_batch': 'restriction_batch_analysis',
            'restriction_map': 'restriction_map',
            'digest': 'restriction_batch_analysis',
            # Sequence manipulation
            'extract_features': 'extract_sequence_features',
            'features': 'extract_sequence_features',
            'transcribe': 'transcribe_sequence',
            'back_transcribe': 'back_transcribe',
            # Codon optimization
            'codon_optimize': 'codon_optimize',
            'optimize_codons': 'codon_optimize',
            # Sequence comparison
            'pairwise_compare': 'pairwise_comparison',
            'compare_sequences': 'pairwise_comparison',
            # Protein analysis
            'secondary_structure': 'predict_secondary_structure',
            'hydrophobicity': 'hydrophobicity_plot',
            # Sequence statistics
            'entropy': 'sequence_entropy',
            'complexity': 'sequence_complexity'
        }
        
        # Tool descriptions for user guidance
        self.tool_descriptions = {
            'translate': 'Translate a DNA sequence into protein',
            'reverse': 'Reverse and complement a DNA sequence',
            'orf': 'Find open reading frames in a sequence',
            'align': 'Global sequence alignment using Needleman-Wunsch',
            'water': 'Local sequence alignment using Smith-Waterman',
            'pattern': 'Search for patterns in sequences',
            'restriction': 'Find restriction enzyme sites',
            'shuffle': 'Shuffle a sequence',
            'info': 'Get information about sequences',
            'sixframe': 'Display all six translation frames',
            'palindrome': 'Find palindromic sequences',
            'fuzzy': 'Fuzzy pattern matching',
            'dotplot': 'Create dot plots between sequences',
            'consensus': 'Generate consensus sequence',
            'gc': 'Calculate GC content of a sequence',
            'download': 'Download sequence from NCBI database',
            'isoelectric': 'Calculate isoelectric point of protein',
            'iep': 'Calculate isoelectric point of protein',
            'cusp': 'Calculate codon usage statistics',
            'pepstats': 'Calculate protein statistics (MW, pI, composition)',
            'wordcount': 'Count oligonucleotide word frequencies',
            'bedtools': 'Find overlaps between genomic regions (BED files)',
            'blat': 'Search sequences against genome (UCSC BLAT)',
            'ucsc_gene': 'Get gene information from UCSC Genome Browser',
            'gtex': 'Get gene expression data across tissues (GTEx)',
            # Biopython tools
            'msa': 'Multiple sequence alignment (Biopython)',
            'pdb': 'Get protein 3D structure information from PDB',
            'protparam': 'Advanced protein analysis (Biopython ProtParam)',
            'phylo': 'Build phylogenetic tree from sequences',
            'motif': 'Find and analyze sequence motifs',
            'primer': 'Primer design and Tm calculation',
            'restriction_batch': 'Batch restriction enzyme analysis',
            'restriction_map': 'Generate restriction map',
            'extract_features': 'Extract sequence features and annotations',
            'transcribe': 'Transcribe DNA to RNA',
            'back_transcribe': 'Back-transcribe RNA to DNA',
            'codon_optimize': 'Optimize codon usage for expression',
            'pairwise_compare': 'Detailed pairwise sequence comparison',
            'secondary_structure': 'Predict protein secondary structure',
            'hydrophobicity': 'Calculate hydrophobicity profile',
            'entropy': 'Calculate sequence entropy/complexity',
            'complexity': 'Analyze sequence complexity'
        }
        
        # Cache for available EMBOSS tools
        self.available_emboss_tools = None
        
        # Cache for AI-resolved tool names
        self.tool_resolution_cache = {}
        
        # Verify EMBOSS installation
        self.check_emboss()
    
    def _resolve_tool_with_ai(self, tool_name: str) -> str:
        """Use Gemini API to resolve unknown tool names to EMBOSS tool names
        
        Args:
            tool_name: User's requested tool name
            
        Returns:
            str: EMBOSS tool name, or original name if not found
        """
        # Check cache first
        if tool_name in self.tool_resolution_cache:
            return self.tool_resolution_cache[tool_name]
        
        try:
            import google.generativeai as genai
            import os
            
            # Configure Gemini
            api_key = os.getenv('GEMINI_API_KEY')
            if not api_key:
                return tool_name  # Can't resolve without API key
            
            genai.configure(api_key=api_key)
            model = genai.GenerativeModel('gemini-2.0-flash-exp')
            
            prompt = f"""You are an EMBOSS bioinformatics tool expert. 

The user requested: "{tool_name}"

What is the exact EMBOSS command-line tool name for this? 

Common EMBOSS tools include:
- needle (global alignment)
- water (local alignment) 
- stretcher (global alignment for longer sequences)
- matcher (local alignment, shows all matches)
- dotmatcher (dot plot)
- transeq (translate DNA to protein)
- revseq (reverse complement)
- getorf (find ORFs)
- restrict (restriction sites)
- pepstats (protein statistics)
- seqret (sequence format conversion)
- trimest, trimseq (sequence trimming)

Respond with ONLY the exact EMBOSS tool name, nothing else. If you don't know, respond with "UNKNOWN"."""

            response = model.generate_content(prompt)
            resolved = response.text.strip().lower()
            
            # Validate response
            if resolved and resolved != "unknown" and len(resolved) < 30 and resolved.isalnum():
                self.tool_resolution_cache[tool_name] = resolved
                return resolved
            else:
                return tool_name
                
        except Exception as e:
            print(f"AI tool resolution failed: {e}")
            return tool_name
    
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
    
    def _log_command(self, tool: str, parameters: dict, result: str, error: Optional[str] = None):
        """Log a command execution for debugging
        
        Args:
            tool: Tool name that was executed
            parameters: Parameters passed to the tool
            result: Result summary or first 500 chars
            error: Error message if command failed
        """
        # Truncate long sequences in parameters for readability
        logged_params = {}
        for key, value in parameters.items():
            if isinstance(value, str) and len(value) > 100:
                logged_params[key] = f"{value[:100]}... ({len(value)} chars total)"
            else:
                logged_params[key] = value
        
        # Truncate result
        logged_result = result[:500] + "..." if len(result) > 500 else result
        
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "tool": tool,
            "parameters": logged_params,
            "result_preview": logged_result,
            "success": error is None,
            "error": error
        }
        
        self.command_log.append(log_entry)
    
    def get_command_log(self) -> List[Dict]:
        """Get the full command log
        
        Returns:
            List of command log entries
        """
        return self.command_log
    
    def get_formatted_log(self) -> str:
        """Get formatted command log as a readable string
        
        Returns:
            Formatted log string
        """
        if not self.command_log:
            return "No commands executed yet."
        
        log_lines = ["=" * 80, "COMMAND EXECUTION LOG", "=" * 80, ""]
        
        for i, entry in enumerate(self.command_log, 1):
            log_lines.append(f"Command #{i} - {entry['timestamp']}")
            log_lines.append(f"Tool: {entry['tool']}")
            log_lines.append(f"Status: {'✓ SUCCESS' if entry['success'] else '✗ FAILED'}")
            log_lines.append("Parameters:")
            for key, value in entry['parameters'].items():
                log_lines.append(f"  - {key}: {value}")
            if entry['error']:
                log_lines.append(f"Error: {entry['error']}")
            else:
                log_lines.append(f"Result Preview: {entry['result_preview']}")
            log_lines.append("")
            log_lines.append("-" * 80)
            log_lines.append("")
        
        return "\n".join(log_lines)
    
    def clear_log(self):
        """Clear the command log"""
        self.command_log = []
    
    def discover_all_emboss_tools(self) -> List[str]:
        """Discover all available EMBOSS tools by checking which have the EMBOSS help signature
        
        Returns:
            list: Names of all available EMBOSS tools (deduplicated)
        """
        if self.available_emboss_tools is not None:
            return self.available_emboss_tools
        
        try:
            # Get the bin directory where EMBOSS tools are located
            result = subprocess.run(
                ['which', 'transeq'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode != 0:
                return []
            
            bin_dir = os.path.dirname(result.stdout.strip())
            tools = set()  # Use set to avoid duplicates
            
            # Check each executable in the bin directory
            for filename in os.listdir(bin_dir):
                filepath = os.path.join(bin_dir, filename)
                if os.path.isfile(filepath) and os.access(filepath, os.X_OK):
                    # Quick check: see if it responds to -help with EMBOSS signature
                    try:
                        result = subprocess.run(
                            [filepath, '-help'],
                            capture_output=True,
                            timeout=0.5,
                            text=True
                        )
                        # EMBOSS tools show "EMBOSS" or the tool version in help
                        if 'EMBOSS' in result.stdout or 'EMBOSS' in result.stderr or 'Version:' in result.stdout:
                            # Remove leading underscore if present (conda adds these)
                            tool_name = filename.lstrip('_')
                            tools.add(tool_name)
                    except (subprocess.TimeoutExpired, Exception):
                        pass
            
            self.available_emboss_tools = sorted(list(tools))
            return self.available_emboss_tools
        
        except Exception as e:
            print(f"Error discovering EMBOSS tools: {e}")
            return []
    
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
    
    def dna_to_rna(self, sequence: str) -> str:
        """Convert DNA sequence to RNA (T → U substitution)
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            str: RNA sequence (with U instead of T)
        """
        try:
            # Simple string replacement is most efficient for T→U
            # But we can also use biosed if needed for more complex replacements
            rna_sequence = sequence.upper().replace('T', 'U')
            return f">rna_sequence\n{rna_sequence}\n"
        
        except Exception as e:
            return f"Error during DNA to RNA conversion: {str(e)}"
    
    def rna_to_dna(self, sequence: str) -> str:
        """Convert RNA sequence to DNA (U → T substitution)
        
        Args:
            sequence: RNA sequence string
        
        Returns:
            str: DNA sequence (with T instead of U)
        """
        try:
            dna_sequence = sequence.upper().replace('U', 'T')
            return f">dna_sequence\n{dna_sequence}\n"
        
        except Exception as e:
            return f"Error during RNA to DNA conversion: {str(e)}"
    
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
    
    def calculate_gc_content(self, sequence: str) -> str:
        """Calculate GC content of a sequence
        
        Args:
            sequence: DNA sequence string
        
        Returns:
            str: GC content percentage and analysis
        """
        try:
            # Remove non-nucleotide characters
            clean_seq = ''.join(c.upper() for c in sequence if c in 'ATGCN')
            
            if len(clean_seq) == 0:
                return "Error: No valid nucleotides found in sequence"
            
            gc_count = clean_seq.count('G') + clean_seq.count('C')
            gc_percent = (gc_count / len(clean_seq)) * 100
            at_count = clean_seq.count('A') + clean_seq.count('T')
            at_percent = (at_count / len(clean_seq)) * 100
            
            result = f"""GC Content Analysis
=====================
Sequence Length: {len(clean_seq)} bp
GC Count: {gc_count}
GC Percentage: {gc_percent:.2f}%
AT Count: {at_count}
AT Percentage: {at_percent:.2f}%

Interpretation:
- Typical bacterial GC content: 30-70%
- Typical human GC content: ~40%
- GC content affects DNA stability and PCR
"""
            return result
        
        except Exception as e:
            return f"Error calculating GC content: {str(e)}"
    
    def download_sequence_from_ncbi(self, accession: str, start: int = None, end: int = None, db: str = 'nucleotide') -> str:
        """Download sequence region from NCBI using accession number and optional coordinates
        
        Args:
            accession: NCBI accession number (e.g., 'NC_000001' for human chr1, 'NM_001234' for mRNA)
            start: Start position (optional, if not specified gets full sequence)
            end: End position (optional, if not specified gets full sequence)
            db: Database to search ('nucleotide' or 'protein')
        
        Returns:
            str: FASTA sequence or error message
        
        Examples:
            - download_sequence_from_ncbi('NC_000001', start=1000, end=2000) -> 1000 bp region
            - download_sequence_from_ncbi('NM_001234') -> full mRNA sequence
        """
        try:
            # Fetch sequence with optional range
            if start and end:
                handle = Entrez.efetch(
                    db=db, 
                    id=accession, 
                    rettype="fasta", 
                    retmode="text",
                    seq_start=start,
                    seq_stop=end
                )
            else:
                handle = Entrez.efetch(
                    db=db, 
                    id=accession, 
                    rettype="fasta", 
                    retmode="text"
                )
            
            record = handle.read()
            handle.close()
            
            if record:
                if start and end:
                    return f"Successfully downloaded {accession} region {start}-{end}:\n\n{record}"
                else:
                    return f"Successfully downloaded {accession}:\n\n{record}"
            else:
                return f"No sequence found for accession: {accession}"
        
        except Exception as e:
            return f"Error downloading from NCBI: {str(e)}. Make sure accession number is valid."
    
    def query_ucsc_genome(self, genome: str, chrom: str, start: int, end: int) -> str:
        """Query UCSC Genome Browser API for sequence data (no storage needed)
        
        Args:
            genome: Genome assembly (e.g., 'hg38' for human, 'mm10' for mouse, 'dm6' for fly)
            chrom: Chromosome (e.g., 'chr1', 'chr2', 'chrX')
            start: Start position (0-based)
            end: End position (0-based, exclusive)
        
        Returns:
            str: FASTA sequence or error message
        
        Examples:
            - query_ucsc_genome('hg38', 'chr1', 1000, 2000) -> 1000 bp from human chr1
            - query_ucsc_genome('mm10', 'chrX', 50000, 51000) -> 1000 bp from mouse chrX
        
        Available genomes: hg38, hg37, mm10, mm9, dm6, dm3, ce10, sacCer3
        """
        try:
            # If single position provided, create small range around it
            if start == end or end == 0:
                end = start + 100  # Get 100bp window
            
            # UCSC DAS server endpoint (no authentication needed, no storage)
            base_url = "https://genome.ucsc.edu/cgi-bin/das"
            
            # Format: /database/dna?segment=chrom:start,end
            url = f"{base_url}/{genome}/dna?segment={chrom}:{start},{end}"
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Parse DAS XML response
            if response.status_code == 200:
                # Extract sequence from DAS format
                import xml.etree.ElementTree as ET
                root = ET.fromstring(response.text)
                
                # Find sequence in XML
                for dna_elem in root.findall('.//{http://www.biodas.org/xmlns/das}DNA'):
                    seq_text = dna_elem.text
                    if seq_text:
                        # Clean whitespace
                        sequence = ''.join(seq_text.split())
                        
                        region_info = f"{genome}:{chrom}:{start}-{end}"
                        return f">UCSC_{region_info}\\n{sequence}"
                
                return f"No sequence data returned for {genome}:{chrom}:{start}-{end}"
            else:
                return f"Error: UCSC API returned status {response.status_code}"
        
        except requests.exceptions.Timeout:
            return "Error: UCSC API request timed out. Try smaller region."
        except Exception as e:
            return f"Error querying UCSC: {str(e)}. Check genome/chromosome names."
    
    def query_gene_info(self, gene_name: str, genome: str = 'hg38', track: str = 'gencode') -> str:
        """Query gene information from Ensembl BioMart
        
        Args:
            gene_name: Gene symbol (e.g., 'ALKBH1', 'TP53', 'BRCA1')
            genome: Genome assembly ('hg38', 'hg37', 'mm10', etc.)
            track: Track/database ('gencode', 'ensembl', 'refseq')
        
        Returns:
            str: Gene information including exons, CDS, transcript length
        
        Examples:
            - query_gene_info('ALKBH1') -> Full gene annotation
            - query_gene_info('TP53', genome='hg38') -> TP53 with all transcripts
        """
        try:
            # Use Ensembl REST API to get gene information
            # First, get the gene ID from gene symbol
            base_url = "https://rest.ensembl.org"
            
            # Search for gene by symbol
            search_url = f"{base_url}/xrefs/symbol/homo_sapiens/{gene_name}?external_db=HGNC"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(search_url, headers=headers, timeout=30)
            response.raise_for_status()
            
            gene_data = response.json()
            
            if not gene_data or len(gene_data) == 0:
                return f"Gene '{gene_name}' not found in Ensembl. Try HGNC symbol (e.g., ALKBH1)"
            
            # Get the Ensembl gene ID (stable ID)
            gene_id = gene_data[0].get('id')
            
            if not gene_id:
                return f"Could not find Ensembl gene ID for {gene_name}"
            
            # Now get detailed information about the gene
            gene_url = f"{base_url}/lookup/id/{gene_id}?expand=1"
            gene_response = requests.get(gene_url, headers=headers, timeout=30)
            gene_response.raise_for_status()
            
            gene_info = gene_response.json()
            
            result = f"Gene Information: {gene_name}\n"
            result += "=" * 70 + "\n"
            result += f"Ensembl Gene ID: {gene_id}\n"
            result += f"Description: {gene_info.get('description', 'N/A')}\n"
            result += f"Biotype: {gene_info.get('biotype', 'N/A')}\n"
            result += f"Location: {gene_info.get('seq_region_name')}:{gene_info.get('start')}-{gene_info.get('end')} ({gene_info.get('strand', '?')})\n"
            result += "\n"
            
            # Process transcripts
            transcripts = gene_info.get('Transcript', [])
            result += f"Total Transcripts: {len(transcripts)}\n\n"
            
            if transcripts:
                for i, transcript in enumerate(transcripts, 1):
                    tx_id = transcript.get('id', 'N/A')
                    tx_biotype = transcript.get('biotype', 'N/A')
                    
                    # Get exons
                    exons = transcript.get('Exon', [])
                    exon_count = len(exons)
                    
                    # Calculate lengths
                    tx_start = transcript.get('start', 0)
                    tx_end = transcript.get('end', 0)
                    tx_length = tx_end - tx_start + 1  # Inclusive
                    
                    # CDS coordinates
                    cds_start = transcript.get('Translation', {}).get('start', tx_start) if transcript.get('Translation') else tx_start
                    cds_end = transcript.get('Translation', {}).get('end', tx_end) if transcript.get('Translation') else tx_end
                    cds_length = cds_end - cds_start + 1 if (cds_start and cds_end) else 0
                    
                    result += f"Transcript {i}: {tx_id}\n"
                    result += f"  Biotype: {tx_biotype}\n"
                    result += f"  Position: {transcript.get('seq_region_name')}:{tx_start}-{tx_end}\n"
                    result += f"  Number of exons: {exon_count}\n"
                    result += f"  Full transcript length (with UTRs): {tx_length:,} bp\n"
                    result += f"  CDS (coding region) length: {cds_length:,} bp\n"
                    
                    if exons:
                        result += f"  Exon coordinates:\n"
                        for j, exon in enumerate(exons, 1):
                            ex_start = exon.get('start')
                            ex_end = exon.get('end')
                            ex_length = ex_end - ex_start + 1
                            result += f"    Exon {j}: {ex_start:,}-{ex_end:,} ({ex_length:,} bp)\n"
                    
                    result += "\n"
            
            return result
        
        except requests.exceptions.Timeout:
            return f"Error: Ensembl API request timed out searching for {gene_name}"
        except requests.exceptions.ConnectionError:
            return f"Error: Could not connect to Ensembl. Check internet connection."
        except Exception as e:
            return f"Error querying gene info from Ensembl: {str(e)}\nExample: query_gene_info('ALKBH1')"

    def get_best_transcript_for_gene(self, gene_name: str, species: str = 'homo_sapiens') -> Optional[str]:
        """Return the best transcript ID for a gene (prefer protein_coding, longest transcript)

        Args:
            gene_name: Gene symbol
            species: species for Ensembl (default homo_sapiens)

        Returns:
            transcript_id (str) or None
        """
        try:
            base_url = "https://rest.ensembl.org"
            headers = {"Content-Type": "application/json"}

            # Get gene stable id via xrefs by symbol
            search_url = f"{base_url}/xrefs/symbol/{species}/{gene_name}?external_db=HGNC"
            r = requests.get(search_url, headers=headers, timeout=30)
            r.raise_for_status()
            xrefs = r.json()
            if not xrefs:
                return None
            gene_id = xrefs[0].get('id')

            # Lookup gene with transcripts expanded
            gene_url = f"{base_url}/lookup/id/{gene_id}?expand=1"
            gr = requests.get(gene_url, headers=headers, timeout=30)
            gr.raise_for_status()
            gene_info = gr.json()

            transcripts = gene_info.get('Transcript', [])
            # Filter protein_coding transcripts
            prot_trans = [t for t in transcripts if t.get('biotype') == 'protein_coding']
            candidates = prot_trans if prot_trans else transcripts

            # Choose longest transcript (by length)
            best = None
            best_len = 0
            for t in candidates:
                start = t.get('start', 0)
                end = t.get('end', 0)
                length = end - start + 1
                if length > best_len:
                    best_len = length
                    best = t.get('id')

            return best
        except Exception:
            return None

    def get_transcript_sequence(self, transcript_id: str, seq_type: str = 'cdna') -> str:
        """Fetch sequence for a transcript from Ensembl

        Args:
            transcript_id: Ensembl transcript ID (e.g., ENST...)
            seq_type: 'cdna', 'cds', or 'genomic'

        Returns:
            sequence string or error message starting with 'Error'
        """
        try:
            base_url = "https://rest.ensembl.org"
            headers = {"Content-Type": "text/plain"}
            seq_url = f"{base_url}/sequence/id/{transcript_id}?type={seq_type}"
            r = requests.get(seq_url, headers=headers, timeout=30)
            r.raise_for_status()
            seq = r.text.strip()
            if not seq:
                return f"Error: Empty sequence returned for {transcript_id}"
            return seq
        except requests.exceptions.Timeout:
            return "Error: Ensembl sequence request timed out"
        except Exception as e:
            return f"Error fetching sequence for {transcript_id}: {str(e)}"
    
    def summarize_gene_info(self, gene_data: str) -> str:
        """Summarize detailed gene information into natural language
        
        Args:
            gene_data: Raw gene information output from query_gene_info()
        
        Returns:
            str: Natural language summary of the gene information
        """
        try:
            # Use Ollama to generate a natural language summary
            from ollama import Client
            
            client = Client(host='http://192.168.128.1:11434')
            
            prompt = f"""You are a helpful bioinformatics assistant. 
            
Summarize this gene information in 2-3 sentences, focusing on the main transcript (usually the first one listed). 
Answer the specific question briefly and naturally.

Gene data:
{gene_data}

Provide a concise, natural language summary suitable for a user query answer. 
Focus on: number of exons, CDS length, and transcript length for the main transcript.
Keep it conversational and friendly."""

            response = client.generate(
                model='gemma3:4b',
                prompt=prompt,
                stream=False
            )
            
            summary = response['response'].strip()
            return summary
        
        except Exception as e:
            # If summarization fails, fall back to extracting key info manually
            lines = gene_data.split('\n')
            
            # Try to extract key information
            summary_lines = []
            for i, line in enumerate(lines):
                if 'Total Transcripts:' in line:
                    summary_lines.append(line.strip())
                if 'Number of exons:' in line and i < 50:  # First occurrence (main transcript)
                    summary_lines.append(line.strip())
                if 'CDS (coding region) length:' in line and i < 50:
                    summary_lines.append(line.strip())
                if 'Full transcript length' in line and i < 50:
                    summary_lines.append(line.strip())
            
            if summary_lines:
                return " | ".join(summary_lines)
            else:
                return "Gene information retrieved. See detailed results below."
    
    def _resolve_sequence_from_gene(self, gene_name: str, seq_type: str = 'cdna') -> Optional[str]:
        """Helper: Given a gene name, return its sequence (or None if failed)
        
        Args:
            gene_name: Gene symbol (e.g., 'ALKBH1')
            seq_type: 'cdna' (default), 'cds', or 'genomic'
        
        Returns:
            sequence string or None
        """
        try:
            transcript_id = self.get_best_transcript_for_gene(gene_name)
            if not transcript_id:
                return None
            seq = self.get_transcript_sequence(transcript_id, seq_type=seq_type)
            if seq.startswith('Error'):
                return None
            return seq
        except Exception:
            return None
    
    def run_blast(self, sequence: str, blast_type: str = 'blastn', database: str = 'nt', 
                  max_results: int = 10, expect_threshold: float = 10.0, exclude_taxa: str = None,
                  word_size: int = None, organism: str = None) -> str:
        """Run BLAST search against NCBI databases using remote Entrez service
        
        Args:
            sequence: Query sequence (DNA, RNA, or protein) or gene name to resolve
            blast_type: 'blastn' (DNA), 'blastp' (protein), 'blastx' (DNA->protein), 'tblastn' (protein->DNA)
            database: NCBI database ('nt' for nucleotide, 'nr' for protein, etc.)
            max_results: Maximum number of results to return
            expect_threshold: E-value threshold for reporting
            exclude_taxa: Taxa to exclude (e.g., 'primates')
            word_size: Word size for BLAST (default varies by program)
            organism: Limit search to specific organism or taxonomic group (e.g., 'Archaea', 'Bacteria')
        
        Returns:
            str: Formatted BLAST results with hits, E-values, and identity percentages
        """
        try:
            from Bio.Blast import NCBIXML, NCBIWWW
            
            # Check if sequence is a gene name (short, all letters) - if so, resolve it
            if sequence and len(sequence) < 20 and sequence.isalpha() and sequence.isupper():
                print(f"Resolving gene {sequence} to sequence...")
                gene_info = self.query_gene_info(sequence)
                if gene_info and not gene_info.startswith('Error'):
                    # Extract sequence from gene info (it's in FASTA format)
                    lines = gene_info.split('\n')
                    seq_lines = [l.strip() for l in lines if l.strip() and not l.startswith('>')]
                    sequence = ''.join(seq_lines)
                    if not sequence:
                        return f"Error: Could not extract sequence from gene {sequence}"
                else:
                    return f"Error: Could not resolve gene {sequence}: {gene_info}"
            
            # Auto-detect protein sequences (contain amino acids not in DNA/RNA)
            # Only convert if it's CLEARLY a protein sequence
            upper_seq = sequence.upper()
            dna_nucleotides = set('ATGCNU')  # DNA/RNA nucleotides
            protein_only_aa = set('EFILPQZ')  # These amino acids don't appear in DNA/RNA
            
            # Count nucleotides vs protein-only amino acids
            nucleotide_count = sum(1 for c in upper_seq if c in dna_nucleotides)
            protein_only_count = sum(1 for c in upper_seq if c in protein_only_aa)
            nucleotide_ratio = nucleotide_count / len(sequence) if len(sequence) > 0 else 0
            
            # Only convert to protein BLAST if:
            # 1. Contains protein-only amino acids AND
            # 2. Less than 95% of sequence is DNA nucleotides
            if protein_only_count > 0 and nucleotide_ratio < 0.95:
                # Definitely a protein sequence
                if blast_type == 'blastn':
                    blast_type = 'blastp'
                    if database == 'nt':
                        database = 'nr'
                    print(f"Auto-detected protein sequence (protein-only AAs: {protein_only_count}, nucleotide ratio: {nucleotide_ratio:.2%}), using blastp")
            elif blast_type == 'blastn':
                print(f"Using blastn for nucleotide sequence (nucleotide ratio: {nucleotide_ratio:.2%})")
            
            print(f"Submitting {blast_type} search to NCBI (this may take 10-30 seconds)...")
            print(f"Query sequence length: {len(sequence)} {'aa' if blast_type == 'blastp' else 'bp'}")
            
            # Normalize database names for NCBI
            db_mapping = {
                'nucleotide': 'nt',
                'nucleotides': 'nt',
                'protein': 'nr',
                'proteins': 'nr',
                'core_nt': 'nt',
                'core_nr': 'nr'
            }
            database = db_mapping.get(database.lower(), database)
            
            # Build Entrez query for taxonomic filtering
            entrez_query = None
            
            # Handle organism parameter (Archaea, Bacteria, specific taxa)
            if organism:
                organism_lower = organism.lower()
                if organism_lower == 'archaea':
                    entrez_query = "txid2157[Organism:exp]"  # Archaea domain
                    print(f"Filtering: Archaea domain only")
                elif organism_lower == 'bacteria':
                    entrez_query = "txid2[Organism:exp]"  # Bacteria domain
                    print(f"Filtering: Bacteria domain only")
                else:
                    # Custom organism filter
                    entrez_query = f"{organism}[Organism]"
                    print(f"Filtering: {organism}")
            
            # Handle exclude_taxa (can combine with organism filter)
            if exclude_taxa and 'primate' in exclude_taxa.lower():
                # Use NCBI taxonomy to exclude all primates (taxid:9443)
                # If no organism filter yet, search all BUT exclude primates
                if entrez_query:
                    entrez_query += " NOT txid9443[Organism:exp]"
                else:
                    # No organism specified, just exclude primates from all results
                    entrez_query = "NOT txid9443[Organism:exp]"
                print(f"Filtering: Excluding primates (NCBI taxid:9443)")
            
            # Request more results when filtering to ensure good hits
            if entrez_query:
                actual_hitlist_size = max(max_results * 10, 200)
            else:
                actual_hitlist_size = max_results
            
            # Build BLAST parameters
            blast_params = {
                'program': blast_type,
                'database': database,
                'sequence': sequence,
                'hitlist_size': actual_hitlist_size,
                'expect': expect_threshold,
                'format_type': 'XML'
            }
            
            # Use megablast for blastn (default, optimized for highly similar sequences)
            if blast_type == 'blastn':
                blast_params['megablast'] = True
                print(f"Using megablast (default for blastn)")
            
            # Add optional parameters
            if entrez_query:
                blast_params['entrez_query'] = entrez_query
            if word_size:
                blast_params['word_size'] = word_size
                print(f"Using word size: {word_size}")
            
            # Submit BLAST query using NCBIWWW (web interface)
            result_handle = NCBIWWW.qblast(**blast_params)
            
            # Parse results
            blast_records = list(NCBIXML.parse(result_handle))
            
            output = []
            output.append(f"BLAST {blast_type.upper()} Results")
            output.append("=" * 60)
            
            if len(blast_records) == 0:
                output.append("No results returned from BLAST search")
                return "\n".join(output)
            
            record = blast_records[0]
            output.append(f"\nQuery: {record.query[:50]}...")
            output.append(f"Query Length: {record.query_length} bp")
            
            if len(record.alignments) == 0:
                output.append(f"Number of Alignments: 0")
                output.append("-" * 60)
                output.append("No matches found")
            else:
                # Note: If entrez_query was used, NCBI already filtered the results
                if entrez_query:
                    output.append(f"Number of Alignments: {len(record.alignments)} (pre-filtered by NCBI taxonomy)")
                    output.append(f"Showing top {min(max_results, len(record.alignments))} results")
                else:
                    output.append(f"Number of Alignments: {len(record.alignments)}")
                output.append("-" * 60)
                
                # Sort by bit score (BLAST default) to match NCBI web interface
                # IMPORTANT: Sort BEFORE taking top N results, not after!
                if entrez_query and len(record.alignments) > 0:
                    def get_bit_score(alignment):
                        if alignment.hsps:
                            return alignment.hsps[0].score
                        return 0
                    
                    # Sort ALL alignments by bit score (BLAST default)
                    all_alignments_sorted = sorted(record.alignments, key=get_bit_score, reverse=True)
                    alignments_to_show = all_alignments_sorted[:max_results]
                    output.append("Results sorted by bit score (BLAST default)")
                    output.append("")
                else:
                    # No taxonomy filtering - use NCBI's default sorting (already by bit score)
                    alignments_to_show = record.alignments[:max_results]
                
                # Debug: Print organisms sorted by identity to console
                if entrez_query and len(record.alignments) > 0:
                    print(f"\nDEBUG: Received {len(record.alignments)} total results from NCBI")
                    # Sort by identity for debug display
                    def get_identity_pct_debug(aln):
                        if aln.hsps:
                            hsp = aln.hsps[0]
                            return (100 * hsp.identities / hsp.align_length) if hsp.align_length > 0 else 0
                        return 0
                    sorted_alns = sorted(record.alignments, key=get_identity_pct_debug, reverse=True)
                    print(f"Top {min(10, len(sorted_alns))} organisms by identity:")
                    for i, aln in enumerate(sorted_alns[:10], 1):
                        # Extract organism name from title
                        title_parts = aln.title.split('|')
                        if len(title_parts) > 4:
                            organism_part = title_parts[4].strip()
                        else:
                            organism_part = aln.title
                        hsp = aln.hsps[0] if aln.hsps else None
                        if hsp:
                            identity_pct = 100 * hsp.identities / hsp.align_length if hsp.align_length > 0 else 0
                            print(f"  {i}. {organism_part[:80]} ({identity_pct:.1f}% identity)")
                    
                    # Check all Cynocephalus volans hits
                    print(f"\nAll Cynocephalus volans hits:")
                    for i, aln in enumerate(record.alignments, 1):
                        if 'Cynocephalus volans' in aln.title:
                            hsp = aln.hsps[0] if aln.hsps else None
                            if hsp:
                                identity_pct = 100 * hsp.identities / hsp.align_length if hsp.align_length > 0 else 0
                                print(f"  {aln.accession}: {identity_pct:.2f}% identity, E-value: {hsp.expect:.2e}")
                
                # Show filtered results
                for alignment_idx, alignment in enumerate(alignments_to_show, 1):
                    output.append(f"\nHit {alignment_idx}: {alignment.title}")
                    output.append(f"  Accession: {alignment.accession}")
                    output.append(f"  Length: {alignment.length} bp")
                    
                    # Show top HSP (High-Scoring Pair)
                    if alignment.hsps:
                        hsp = alignment.hsps[0]
                        output.append(f"  E-value: {hsp.expect:.2e}")
                        output.append(f"  Score: {hsp.score} bits")
                        identity_pct = 100 * hsp.identities // hsp.align_length if hsp.align_length > 0 else 0
                        output.append(f"  Identity: {hsp.identities}/{hsp.align_length} ({identity_pct}%)")
                        output.append(f"  Query alignment: {hsp.query[:70]}")
                        output.append(f"  Subject alignment: {hsp.sbjct[:70]}")
            
            return "\n".join(output)
            
        except Exception as e:
            import traceback
            error_detail = traceback.format_exc()
            return f"BLAST search failed: {str(e)}\n\nDetails: {error_detail}\n\nMake sure you have an internet connection and BioPython is properly installed."

    def run_tool(self, tool_name: str, **kwargs) -> str:
        """Generic method to run any EMBOSS tool or BLAST
        
        Supports:
        - Specific hardcoded implementations for common tools (transeq, revseq, etc.)
        - Generic fallback for any other EMBOSS tool (iep, charge, etc.)
        - Gene-based access: any tool can accept gene_name and auto-resolve to sequence
        - Transcript variant selection: specify which transcript variant to use
        - BLAST searches via NCBI
        
        Args:
            tool_name: Name of the tool to run (natural language or EMBOSS name)
            **kwargs: Tool-specific parameters (can include 'gene_name', 'sequence', 'transcript_variant')
        
        Returns:
            str: Tool output or error message
        """
        # Log the command execution
        try:
            result = self._run_tool_internal(tool_name, **kwargs)
            self._log_command(tool_name, kwargs, result)
            return result
        except RecursionError:
            error_msg = f"Error: Maximum recursion depth exceeded while running {tool_name}"
            self._log_command(tool_name, kwargs, "", error=error_msg)
            return error_msg
        except Exception as e:
            error_msg = f"Error: {str(e)}"
            self._log_command(tool_name, kwargs, "", error=error_msg)
            return error_msg
    
    def _run_tool_internal(self, tool_name: str, **kwargs) -> str:
        """Internal implementation of run_tool (without logging wrapper)
        
        Args:
            tool_name: Name of the tool to run
            **kwargs: Tool-specific parameters
        
        Returns:
            str: Tool output or error message
        """
        # Map natural language name to EMBOSS command if needed
        emboss_name = self.tool_map.get(tool_name.lower(), tool_name.lower())
        
        # If tool not in map, try asking Gemini for the correct EMBOSS tool name
        if emboss_name == tool_name.lower() and emboss_name not in self.tool_map.values():
            resolved_tool = self._resolve_tool_with_ai(tool_name)
            if resolved_tool and resolved_tool != tool_name:
                print(f"AI resolved '{tool_name}' → '{resolved_tool}'")
                emboss_name = resolved_tool
        
        # If a gene_name is provided but no sequence, try to resolve the sequence first
        if ('gene_name' in kwargs or 'gene' in kwargs) and 'sequence' not in kwargs:
            gene = kwargs.get('gene_name') or kwargs.get('gene')
            transcript_variant = kwargs.get('transcript_variant', None)
            
            # If a specific transcript variant is requested, fetch that one
            if transcript_variant:
                # Extract variant number from "transcript variant N"
                import re
                match = re.search(r'variant\s+(\d+)', transcript_variant.lower())
                if match:
                    variant_num = int(match.group(1))
                    # Query gene info and get specific transcript
                    gene_info = self.query_gene_info(gene)
                    if gene_info and not gene_info.startswith('Error'):
                        # Extract transcript from the numbered position
                        lines = gene_info.split('\n')
                        transcript_id = None
                        for i, line in enumerate(lines):
                            # Check if this line has "Transcript N: ENSTNNN..."
                            if f'Transcript {variant_num}:' in line:
                                # The transcript ID might be on this line
                                parts = line.split()
                                for part in parts:
                                    if part.startswith('ENST'):
                                        transcript_id = part
                                        break
                                break
                        
                        if transcript_id:
                            seq = self.get_transcript_sequence(transcript_id, seq_type='cds')
                            if seq and not seq.startswith('Error'):
                                kwargs['sequence'] = seq
                                kwargs['_gene_name'] = gene
                                kwargs['_transcript_variant'] = transcript_variant
                            else:
                                return f"Could not fetch sequence for {gene} {transcript_variant}"
                        else:
                            return f"Could not find {transcript_variant} in {gene} gene info"
            else:
                # No specific variant requested, use best transcript
                seq = self._resolve_sequence_from_gene(gene)
                if seq:
                    kwargs['sequence'] = seq
                    # Also keep gene name for context in output
                    kwargs['_gene_name'] = gene
                else:
                    return f"Could not resolve sequence for gene: {gene}"
        
        # Route to specific hardcoded methods (optimized implementations)
        if emboss_name == 'transeq':
            result = self.translate_sequence(kwargs.get('sequence', ''), kwargs.get('frame', 1))
            # Add gene/transcript context if available
            if '_gene_name' in kwargs:
                gene = kwargs['_gene_name']
                variant_info = f" from {kwargs['_transcript_variant']}" if '_transcript_variant' in kwargs else ""
                result = f"Translation of {gene}{variant_info}:\n\n{result}"
                # Count amino acids
                lines = result.split('\n')
                aa_count = len([c for c in result if c.isalpha() and c.islower()])
                if aa_count == 0:
                    for line in lines:
                        if line and not line.startswith('>'):
                            aa_count += len([c for c in line if c.isalpha()])
                            break
            return result
        elif emboss_name == 'revseq':
            return self.reverse_complement(kwargs.get('sequence', ''))
        elif emboss_name in ['dna_to_rna', 'dna2rna', 'torna', 'rna']:
            return self.dna_to_rna(kwargs.get('sequence', ''))
        elif emboss_name in ['rna_to_dna', 'rna2dna', 'todna']:
            return self.rna_to_dna(kwargs.get('sequence', ''))
        elif emboss_name == 'getorf':
            return self.find_orfs(kwargs.get('sequence', ''), kwargs.get('min_size', 100))
        elif emboss_name == 'infoseq':
            return self.get_sequence_info(kwargs.get('sequence', ''))
        elif emboss_name == 'restrict':
            return self.find_restriction_sites(kwargs.get('sequence', ''), kwargs.get('enzyme', None))
        elif emboss_name == 'sixpack':
            return self.get_six_frame_translation(kwargs.get('sequence', ''))
        elif emboss_name == 'gc_content':
            result = self.calculate_gc_content(kwargs.get('sequence', ''))
            # If we resolved from a gene, add context
            if '_gene_name' in kwargs:
                gene = kwargs['_gene_name']
                result = f"Gene {gene} GC Content:\n\n{result}"
            return result
        elif emboss_name == 'download_sequence':
            # Check if using UCSC API (genome-based)
            if 'genome' in kwargs:
                return self.query_ucsc_genome(
                    kwargs.get('genome', ''),
                    kwargs.get('chrom', ''),
                    kwargs.get('start', 0),
                    kwargs.get('end', 0)
                )
            else:
                # Fall back to NCBI
                return self.download_sequence_from_ncbi(
                    kwargs.get('accession', ''),
                    kwargs.get('start', None),
                    kwargs.get('end', None),
                    kwargs.get('db', 'nucleotide')
                )
        elif emboss_name == 'gene_query':
            # Query gene information from Ensembl
            return self.query_gene_info(kwargs.get('gene_name', ''), kwargs.get('genome', 'hg38'), kwargs.get('track', 'gencode'))
        elif emboss_name in ['blast', 'blastn', 'blastp', 'blastx', 'search']:
            # Route BLAST requests to remote NCBI service
            blast_type = 'blastn' if emboss_name in ['blast', 'search'] else emboss_name
            
            # Get sequence - may be provided directly or as a gene name to resolve
            sequence = kwargs.get('sequence', '')
            if not sequence and ('gene_name' in kwargs or 'gene' in kwargs):
                # Resolve gene to sequence first
                gene = kwargs.get('gene_name') or kwargs.get('gene')
                gene_info = self.query_gene_info(gene)
                if gene_info and not gene_info.startswith('Error'):
                    # Extract sequence from gene info (it's typically in FASTA format)
                    lines = gene_info.split('\n')
                    seq_lines = [l for l in lines if l and not l.startswith('>')]
                    sequence = ''.join(seq_lines)
            
            return self.run_blast(
                sequence,
                blast_type=blast_type,
                database=kwargs.get('database', 'nt'),
                max_results=int(kwargs.get('max_results', 10)),
                expect_threshold=float(kwargs.get('expect_threshold', 10.0)),
                exclude_taxa=kwargs.get('exclude_taxa', None),
                word_size=kwargs.get('word_size', None),
                organism=kwargs.get('organism', None)
            )
        elif emboss_name == 'bedtools_intersect':
            # BEDTools intersect for genomic overlap analysis
            result = self.bedtools_intersect(
                kwargs.get('file_a', ''),
                kwargs.get('file_b', ''),
                kwargs.get('output_format', 'bed')
            )
            if result.get('success'):
                return result.get('output', '') + f"\n\n{result.get('summary', '')}"
            else:
                return f"BEDTools error: {result.get('error', 'Unknown error')}"
        elif emboss_name == 'blat_search':
            # BLAT search via UCSC
            result = self.blat_search(
                kwargs.get('sequence', ''),
                kwargs.get('database', 'hg38'),
                kwargs.get('output_format', 'psl')
            )
            if result.get('success'):
                return result.get('output', '') + f"\n\n{result.get('summary', '')}\n{result.get('note', '')}"
            else:
                return f"BLAT error: {result.get('error', 'Unknown error')}"
        elif emboss_name == 'ucsc_gene_info':
            # UCSC gene information
            result = self.ucsc_gene_info(
                kwargs.get('gene_name', ''),
                kwargs.get('database', 'hg38')
            )
            if result.get('success'):
                output = f"Gene: {result.get('gene', '')}\n"
                output += f"Database: {result.get('database', '')}\n"
                output += f"Note: {result.get('note', '')}\n"
                output += f"Browser link: {result.get('api_endpoint', '')}"
                return output
            else:
                return f"UCSC error: {result.get('error', 'Unknown error')}"
        elif emboss_name == 'gtex_expression':
            # GTEx tissue expression data
            result = self.gtex_expression(
                kwargs.get('gene_name', ''),
                kwargs.get('top_n', 10)
            )
            if result.get('success'):
                return result.get('output', '')
            else:
                return f"GTEx error: {result.get('error', 'Unknown error')}"
        elif emboss_name == 'ucsc_table_query':
            # UCSC Table Browser query
            return self.ucsc_table_query(
                genome=kwargs.get('genome', 'hg38'),
                track=kwargs.get('track', ''),
                table=kwargs.get('table', None),
                chrom=kwargs.get('chrom', None),
                start=kwargs.get('start', None),
                end=kwargs.get('end', None),
                filter_field=kwargs.get('filter_field', None),
                filter_value=kwargs.get('filter_value', None)
            )
        elif emboss_name == 'pubmed_search':
            # PubMed literature search
            return self.pubmed_search(
                query=kwargs.get('query', ''),
                max_results=int(kwargs.get('max_results', 10)),
                year=kwargs.get('year', None)
            )
        elif emboss_name == 'track_intersection':
            # Track intersection/overlap analysis
            return self.track_intersection(
                genome=kwargs.get('genome', 'hg38'),
                track1=kwargs.get('track1', ''),
                track2=kwargs.get('track2', ''),
                chrom=kwargs.get('chrom', None),
                max_results=int(kwargs.get('max_results', 10))
            )
        elif emboss_name == 'find_neighboring_genes':
            # Find genes neighboring a target gene
            return self.find_neighboring_genes(
                gene_name=kwargs.get('gene_name', ''),
                genome=kwargs.get('genome', 'hg38'),
                track=kwargs.get('track', 'knownGene')
            )
        # Biopython-powered tools
        elif emboss_name == 'multiple_sequence_alignment':
            sequences = kwargs.get('sequences', [])
            if not sequences:
                return "Error: Provide 'sequences' as list of (name, seq) tuples"
            return self.multiple_sequence_alignment(sequences)
        elif emboss_name == 'protein_structure':
            pdb_id = kwargs.get('pdb_id', kwargs.get('pdb', ''))
            if not pdb_id:
                return "Error: Provide 'pdb_id' (e.g., '1ABC')"
            return self.get_protein_structure_info(pdb_id)
        elif emboss_name == 'molecular_weight_biopython':
            sequence = kwargs.get('sequence', '')
            seq_type = kwargs.get('seq_type', 'protein')
            return self.calculate_molecular_weight_biopython(sequence, seq_type)
        # Phylogenetics
        elif emboss_name == 'phylogenetic_tree':
            sequences = kwargs.get('sequences', [])
            if not sequences:
                return "Error: Provide 'sequences' as list of (name, seq) tuples"
            return self.phylogenetic_tree(sequences)
        # Motif analysis
        elif emboss_name == 'find_motifs':
            sequences = kwargs.get('sequences', [])
            if not sequences:
                return "Error: Provide 'sequences' as list of (name, seq) tuples"
            return self.find_motifs(sequences)
        elif emboss_name == 'motif_consensus':
            sequences = kwargs.get('sequences', [])
            if not sequences:
                return "Error: Provide 'sequences' as list of (name, seq) tuples"
            return self.motif_consensus(sequences)
        # Primer analysis
        elif emboss_name == 'primer_analysis':
            sequence = kwargs.get('sequence', '')
            if not sequence:
                return "Error: Provide primer 'sequence'"
            return self.primer_analysis(sequence)
        # Restriction enzymes
        elif emboss_name == 'restriction_batch_analysis':
            sequence = kwargs.get('sequence', '')
            enzymes = kwargs.get('enzymes', None)
            if not sequence:
                return "Error: Provide DNA 'sequence'"
            return self.restriction_batch_analysis(sequence, enzymes)
        elif emboss_name == 'restriction_map':
            sequence = kwargs.get('sequence', '')
            if not sequence:
                return "Error: Provide DNA 'sequence'"
            return self.restriction_map(sequence)
        # Sequence manipulation
        elif emboss_name == 'extract_sequence_features':
            sequence = kwargs.get('sequence', '')
            seq_type = kwargs.get('seq_type', 'dna')
            return self.extract_sequence_features(sequence, seq_type)
        elif emboss_name == 'transcribe_sequence':
            sequence = kwargs.get('sequence', '')
            return self.transcribe_sequence(sequence)
        elif emboss_name == 'back_transcribe':
            sequence = kwargs.get('sequence', '')
            return self.back_transcribe(sequence)
        elif emboss_name == 'codon_optimize':
            sequence = kwargs.get('sequence', '')
            organism = kwargs.get('organism', 'human')
            return self.codon_optimize(sequence, organism)
        # Sequence comparison
        elif emboss_name == 'pairwise_comparison':
            seq1 = kwargs.get('seq1', '')
            seq2 = kwargs.get('seq2', '')
            if not (seq1 and seq2):
                return "Error: Provide 'seq1' and 'seq2'"
            return self.pairwise_comparison(seq1, seq2)
        # Protein analysis
        elif emboss_name == 'predict_secondary_structure':
            sequence = kwargs.get('sequence', '')
            return self.predict_secondary_structure(sequence)
        elif emboss_name == 'hydrophobicity_plot':
            sequence = kwargs.get('sequence', '')
            return self.hydrophobicity_plot(sequence)
        # Sequence statistics
        elif emboss_name == 'sequence_entropy':
            sequence = kwargs.get('sequence', '')
            return self.sequence_entropy(sequence)
        elif emboss_name == 'sequence_complexity':
            sequence = kwargs.get('sequence', '')
            return self.sequence_complexity(sequence)
        else:
            # Generic EMBOSS tool fallback
            return self._run_generic_emboss_tool(emboss_name, **kwargs)
    
    def _run_generic_emboss_tool(self, tool_name: str, **kwargs) -> str:
        """Generic runner for any EMBOSS tool that accepts sequence input
        
        Supports both single-sequence and two-sequence tools:
        - Single sequence: provide 'sequence' parameter
        - Two sequences: provide 'seq1' and 'seq2' parameters
        
        Args:
            tool_name: EMBOSS tool name (e.g., 'iep', 'charge', 'dotmatcher')
            **kwargs: Parameters including 'sequence' OR ('seq1' and 'seq2')
        
        Returns:
            str: Tool output or error message
        """
        # Check if this is a two-sequence tool (seq1 and seq2 provided)
        seq1 = kwargs.get('seq1', '')
        seq2 = kwargs.get('seq2', '')
        sequence = kwargs.get('sequence', '')
        
        # Determine if single or two-sequence tool
        is_two_sequence = bool(seq1 and seq2)
        
        if not is_two_sequence and not sequence:
            return f"Error: No sequence provided for tool '{tool_name}'. Provide 'sequence' for single-sequence tools or 'seq1' and 'seq2' for two-sequence tools."
        
        try:
            if is_two_sequence:
                # Two-sequence tool (dotplot, needle, water, etc.)
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f1:
                    f1.write(f">seq1\n{seq1}\n")
                    file1 = f1.name
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f2:
                    f2.write(f">seq2\n{seq2}\n")
                    file2 = f2.name
                
                # Check if this is a graphics tool (dotmatcher, dotpath, dottup, polydot)
                is_graphics_tool = tool_name in ['dotmatcher', 'dotpath', 'dottup', 'polydot']
                
                if is_graphics_tool:
                    # Graphics tools use -graph instead of -outfile
                    output_file = tempfile.NamedTemporaryFile(suffix='.png', delete=False).name
                    cmd = [tool_name, '-asequence', file1, '-bsequence', file2, '-graph', 'png', '-goutfile', output_file]
                else:
                    # Text output tools use -outfile
                    output_file = tempfile.NamedTemporaryFile(suffix=f'_{tool_name}_output.txt', delete=False).name
                    cmd = [tool_name, '-asequence', file1, '-bsequence', file2, '-outfile', output_file]
                
                # Add any additional parameters
                skip_keys = {'seq1', 'seq2', '_gene_name', 'gene_name', 'gene'}
                for key, value in kwargs.items():
                    if key not in skip_keys and value is not None:
                        param_name = f'-{key.replace("_", "")}'
                        cmd.extend([param_name, str(value)])
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                
                if is_graphics_tool:
                    # Graphics tools like dotmatcher append .1.png to the filename
                    # Check for the actual generated file
                    actual_output_file = None
                    
                    # Try common patterns
                    possible_files = [
                        output_file,  # exact name
                        f"{output_file}.1.png",  # dotmatcher pattern
                        output_file.replace('.png', '.1.png')  # another pattern
                    ]
                    
                    for pfile in possible_files:
                        if os.path.exists(pfile) and os.path.getsize(pfile) > 0:
                            actual_output_file = pfile
                            break
                    
                    if actual_output_file:
                        # Cleanup temp sequence files only
                        try:
                            os.remove(file1)
                            os.remove(file2)
                        except:
                            pass
                        return f"IMAGE_FILE:{actual_output_file}"
                    else:
                        error_msg = result.stderr if result.stderr else f"Graphics file not generated. Checked: {possible_files}"
                        return f"Error running {tool_name}: {error_msg}"
                
                elif result.returncode == 0 and os.path.exists(output_file):
                    # Text output tools
                    with open(output_file, 'r') as f:
                        output = f.read()
                    # Cleanup
                    try:
                        os.remove(file1)
                        os.remove(file2)
                        os.remove(output_file)
                    except:
                        pass
                    return output
                else:
                    error_msg = result.stderr if result.stderr else "Tool execution failed"
                    return f"Error running {tool_name}: {error_msg}"
            
            else:
                # Single-sequence tool (iep, charge, gc, etc.)
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                    f.write(f">input_seq\n{sequence}\n")
                    input_file = f.name
                
                output_file = input_file.replace('.fasta', f'_{tool_name}_output.txt')
                
                # Build command with the tool and basic parameters
                cmd = [tool_name, '-sequence', input_file, '-outfile', output_file]
                
                # Add any additional parameters from kwargs (excluding internal/non-EMBOSS keys)
                skip_keys = {'sequence', '_gene_name', 'gene_name', 'gene'}
                for key, value in kwargs.items():
                    if key not in skip_keys and value is not None:
                        # Convert kwargs to EMBOSS format (-param value)
                        param_name = f'-{key.replace("_", "")}'
                        cmd.extend([param_name, str(value)])
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=30
                )
                
                if result.returncode == 0 and os.path.exists(output_file):
                    with open(output_file, 'r') as f:
                        output = f.read()
                    
                    # Add gene context if available
                    if '_gene_name' in kwargs:
                        gene = kwargs['_gene_name']
                        output = f"Gene {gene} - {tool_name.upper()}:\n\n{output}"
                    
                    # Cleanup
                    try:
                        os.remove(input_file)
                        os.remove(output_file)
                    except:
                        pass
                    
                    return output
                else:
                    # If no output file, return stdout/stderr
                    error_msg = result.stderr if result.stderr else "Tool execution failed"
                    return f"Error running {tool_name}: {error_msg}"
        
        except Exception as e:
            return f"Error: Failed to run {tool_name}: {str(e)}"
    
    def execute_multi_step(self, steps: List[Dict], use_previous_result: bool = False) -> Tuple[bool, List[Dict]]:
        """Execute multiple analysis steps in sequence with caching
        
        Args:
            steps: List of step dicts containing "tool" and "parameters"
            use_previous_result: If True, feed output from step N to step N+1
        
        Returns:
            Tuple of (success: bool, results: list of result dicts with step info and output)
        """
        results = []
        cached_result = None
        cached_sequence = None
        previous_tool = None
        previous_gene_name = None  # Track gene name from gene_query
        previous_sequence_type = 'cdna'  # Default sequence type
        
        try:
            for idx, step in enumerate(steps):
                tool_name = step.get('tool')
                parameters = step.get('parameters', {}).copy()
                
                # Track gene name and sequence type from gene_query step
                if tool_name == 'gene_query':
                    previous_gene_name = parameters.get('gene_name') or parameters.get('gene')
                    # Store the sequence type requested (mRNA, protein, etc.)
                    seq_type_param = parameters.get('sequence_type', 'mRNA').lower()
                    if seq_type_param in ['mrna', 'transcript', 'cdna']:
                        previous_sequence_type = 'cdna'
                    elif seq_type_param in ['protein', 'peptide']:
                        previous_sequence_type = 'protein'
                    else:
                        # Default: if user asks for both mRNA and protein, prioritize mRNA for BLAST
                        previous_sequence_type = 'cdna'
                    print(f"[Step {idx+1}/{len(steps)}] Tracked gene_name={previous_gene_name}, sequence_type={previous_sequence_type}")
                
                # Smart chaining: if previous step was gene_query and this step needs sequence
                if previous_tool == 'gene_query' and tool_name in ['blast', 'blastn', 'blastp', 'blastx', 'search', 'gc', 'translate', 'reverse', 'orf', 'pepstats', 'iep', 'cusp']:
                    seq_type = previous_sequence_type  # Use the sequence type from gene_query
                    
                    # Check if sequence parameter indicates it should come from previous step
                    sequence_param = parameters.get('sequence', '')
                    needs_fetch = False
                    
                    if not sequence_param or 'sequence' not in parameters:
                        needs_fetch = True
                    elif isinstance(sequence_param, str):
                        # Check for various ways Gemini might say "use previous result"
                        lower_seq = sequence_param.lower()
                        if any(phrase in lower_seq for phrase in [
                            'use_previous_result', 'previous step', 'from previous', 
                            'use previous', 'from step', 'retrieved', 'from gene_query'
                        ]):
                            needs_fetch = True
                    
                    if needs_fetch and previous_gene_name:
                        print(f"[Step {idx+1}/{len(steps)}] Fetching {seq_type} sequence for {previous_gene_name}...")
                        # Use _resolve_sequence_from_gene which handles gene name -> transcript -> sequence
                        sequence = self._resolve_sequence_from_gene(previous_gene_name, seq_type=seq_type)
                        if sequence:
                            # Clean up the sequence (remove formatting)
                            sequence = ''.join(c for c in sequence if c.isalpha())
                            parameters['sequence'] = sequence
                            print(f"[Step {idx+1}/{len(steps)}] Using {seq_type} sequence ({len(sequence)} {'aa' if seq_type == 'protein' else 'bp'}) for {tool_name}")
                        else:
                            return False, [{
                                'step': idx + 1,
                                'error': f"Could not fetch {seq_type} sequence for {previous_gene_name}",
                                'details': "Failed to resolve gene name to transcript sequence"
                            }]
                
                # If using previous result, inject it into current step
                if use_previous_result and cached_result and idx > 0:
                    # Try to use the sequence output from previous step
                    if 'sequence' not in parameters and '_sequence' not in parameters:
                        # Extract just the sequence if previous output is FASTA format
                        clean_sequence = cached_result
                        if '>' in cached_result:
                            # FASTA format - extract just the sequence lines
                            lines = cached_result.split('\n')
                            seq_lines = [line for line in lines if line and not line.startswith('>')]
                            clean_sequence = ''.join(seq_lines).strip()
                        else:
                            # Try to extract just alphanumeric sequence characters
                            import re
                            seq_match = re.search(r'([A-Z]{10,})', cached_result, re.MULTILINE)
                            if seq_match:
                                clean_sequence = seq_match.group(1)
                        
                        parameters['sequence'] = clean_sequence
                        print(f"[Step {idx+1}/{len(steps)}] Using previous result as sequence ({len(clean_sequence)} characters)")
                
                # Execute the tool
                print(f"[Step {idx+1}/{len(steps)}] Running {tool_name}...")
                try:
                    output = self.run_tool(tool_name, **parameters)
                    
                    # Check if output indicates an error
                    if output.startswith('Error') or output.startswith('BLAST search failed'):
                        results.append({
                            'step': idx + 1,
                            'tool': tool_name,
                            'success': False,
                            'error': output,
                            'parameters': parameters
                        })
                        # Stop on first error
                        return False, results
                    
                    # Cache the result for potential use in next step
                    cached_result = output
                    previous_tool = tool_name
                except Exception as e:
                    results.append({
                        'step': idx + 1,
                        'tool': tool_name,
                        'success': False,
                        'error': str(e),
                        'parameters': parameters
                    })
                    return False, results
                
                results.append({
                    'step': idx + 1,
                    'tool': tool_name,
                    'success': True,
                    'output': output,
                    'parameters': parameters
                })
            
            return True, results
        
        except Exception as e:
            results.append({
                'step': len(results) + 1,
                'error': f"Multi-step execution failed: {str(e)}"
            })
            return False, results
    
    def format_multi_step_results(self, results: List[Dict], explanation: str = "") -> str:
        """Format multi-step results for display
        
        Args:
            results: List of result dicts from execute_multi_step
            explanation: Overall explanation of the workflow
        
        Returns:
            Formatted results string
        """
        output = ""
        
        if explanation:
            output += f"📋 Workflow: {explanation}\n"
            output += "=" * 60 + "\n\n"
        
        for result in results:
            step_num = result.get('step', '?')
            tool = result.get('tool', 'unknown')
            
            output += f"Step {step_num}: {tool.upper()}\n"
            output += "-" * 40 + "\n"
            
            if result.get('success'):
                params = result.get('parameters', {})
                if params:
                    output += "Parameters: "
                    param_strs = []
                    for k, v in params.items():
                        if isinstance(v, str) and len(v) > 50:
                            param_strs.append(f"{k}={v[:50]}...")
                        else:
                            param_strs.append(f"{k}={v}")
                    output += ", ".join(param_strs) + "\n"
                
                result_text = result.get('output', '')
                if isinstance(result_text, str) and len(result_text) > 50000:
                    output += f"Output:\n{result_text[:50000]}...\n(Output truncated at 50000 characters)\n\n"
                else:
                    output += f"Output:\n{result_text}\n\n"
            else:
                output += f"❌ Error: {result.get('error', 'Unknown error')}\n\n"
        
        return output


    def bedtools_intersect(self, file_a: str, file_b: str, output_format: str = "bed") -> Dict:
        """Find overlaps between two BED files using BEDTools
        
        Args:
            file_a: Path to first BED file or BED content as string
            file_b: Path to second BED file or BED content as string
            output_format: Output format (bed, count, etc.)
        
        Returns:
            Dict with success status and results
        """
        try:
            # Check if BEDTools is installed
            try:
                subprocess.run(['bedtools', '--version'], capture_output=True, check=True)
            except (FileNotFoundError, subprocess.CalledProcessError):
                return {
                    'success': False,
                    'error': 'BEDTools not installed. Install with: conda install -c bioconda bedtools'
                }
            
            # Create temp files if content is provided as strings
            temp_files = []
            
            if not os.path.exists(file_a):
                with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
                    f.write(file_a)
                    file_a = f.name
                    temp_files.append(file_a)
            
            if not os.path.exists(file_b):
                with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
                    f.write(file_b)
                    file_b = f.name
                    temp_files.append(file_b)
            
            # Run bedtools intersect
            cmd = ['bedtools', 'intersect', '-a', file_a, '-b', file_b]
            
            if output_format == 'count':
                cmd.append('-c')
            elif output_format == 'summary':
                cmd.extend(['-wo'])  # Write original A entry + overlap
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            # Cleanup temp files
            for tf in temp_files:
                try:
                    os.unlink(tf)
                except:
                    pass
            
            if result.returncode == 0:
                output = result.stdout.strip()
                lines = output.split('\n') if output else []
                
                return {
                    'success': True,
                    'output': output,
                    'count': len(lines),
                    'summary': f"Found {len(lines)} overlapping regions"
                }
            else:
                return {
                    'success': False,
                    'error': f"BEDTools error: {result.stderr}"
                }
        
        except Exception as e:
            return {
                'success': False,
                'error': f"BEDTools intersect failed: {str(e)}"
            }


    def blat_search(self, sequence: str, database: str = "hg38", output_format: str = "psl") -> Dict:
        """Search sequence against genome using BLAT via UCSC API
        
        Args:
            sequence: DNA or protein sequence to search
            database: Genome assembly (hg38, mm10, etc.)
            output_format: Output format (psl, json)
        
        Returns:
            Dict with success status and BLAT results
        """
        try:
            # UCSC BLAT API endpoint
            url = "https://genome.ucsc.edu/cgi-bin/hgBlat"
            
            params = {
                'userSeq': sequence,
                'type': 'BLAT',
                'db': database,
                'output': output_format
            }
            
            response = requests.post(url, data=params, timeout=60)
            
            if response.status_code == 200:
                result = response.text
                
                # Parse PSL format to extract top hits with coordinates
                lines = result.split('\n')
                data_lines = [l for l in lines if l and not l.startswith('#') and not l.startswith('psLayout') and '\t' in l]
                
                # Parse top hits (PSL format: matches, mismatches, ..., chrom, start, end, ...)
                parsed_hits = []
                for line in data_lines[:10]:  # Top 10 hits
                    fields = line.split('\t')
                    if len(fields) >= 21:
                        try:
                            matches = int(fields[0])
                            mismatches = int(fields[1])
                            chrom = fields[13]
                            start = int(fields[15])
                            end = int(fields[16])
                            strand = fields[8]
                            
                            # Calculate identity
                            identity = (matches / (matches + mismatches) * 100) if (matches + mismatches) > 0 else 0
                            
                            parsed_hits.append({
                                'chrom': chrom,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'matches': matches,
                                'identity': f"{identity:.1f}%",
                                'ucsc_link': f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={database}&position={chrom}:{start}-{end}"
                            })
                        except (ValueError, IndexError):
                            continue
                
                # Build summary output
                summary_lines = [f"BLAT search found {len(data_lines)} hits in {database}"]
                if parsed_hits:
                    summary_lines.append("\nTop hits:")
                    for i, hit in enumerate(parsed_hits[:5], 1):
                        summary_lines.append(f"{i}. {hit['chrom']}:{hit['start']}-{hit['end']} ({hit['identity']} identity, {hit['matches']} bp)")
                        summary_lines.append(f"   View in UCSC Browser: {hit['ucsc_link']}")
                
                return {
                    'success': True,
                    'output': result,
                    'num_hits': len(data_lines),
                    'parsed_hits': parsed_hits,
                    'summary': "\n".join(summary_lines),
                    'note': "BLAT is more sensitive for near-exact matches than BLAST. Click UCSC Browser links to view genes, conservation, chromHMM, and ChIP-seq tracks."
                }
            else:
                return {
                    'success': False,
                    'error': f"BLAT search failed: HTTP {response.status_code}"
                }
        
        except Exception as e:
            return {
                'success': False,
                'error': f"BLAT search failed: {str(e)}"
            }


    def ucsc_gene_info(self, gene_name: str, database: str = "hg38") -> Dict:
        """Get detailed gene information from UCSC Genome Browser
        
        Args:
            gene_name: Gene symbol or name
            database: Genome assembly (hg38, mm10, etc.)
        
        Returns:
            Dict with gene information including position, exons, strand
        """
        try:
            # UCSC Table Browser API
            url = "https://api.genome.ucsc.edu/getData/track"
            
            params = {
                'genome': database,
                'track': 'ncbiRefSeq',
                'chrom': 'chr1',  # Will be updated based on search
            }
            
            # First, search for the gene
            search_url = f"https://api.genome.ucsc.edu/list/tracks?genome={database}"
            
            response = requests.get(search_url, timeout=30)
            
            if response.status_code == 200:
                # For now, return a structured response indicating capability
                return {
                    'success': True,
                    'gene': gene_name,
                    'database': database,
                    'note': 'UCSC API access configured. Use UCSC Genome Browser for detailed visualization.',
                    'api_endpoint': f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={database}&position={gene_name}"
                }
            else:
                return {
                    'success': False,
                    'error': f"UCSC API error: HTTP {response.status_code}"
                }
        
        except Exception as e:
            return {
                'success': False,
                'error': f"UCSC query failed: {str(e)}"
            }

    def find_neighboring_genes(self, gene_name: str, genome: str = "hg38", track: str = "knownGene") -> str:
        """Find genes neighboring a target gene on the chromosome
        
        Args:
            gene_name: Target gene symbol (e.g., 'CARS1')
            genome: Genome assembly (default: 'hg38')
            track: Gene track to use (default: 'knownGene' for GENCODE)
        
        Returns:
            Formatted string with neighboring genes and distances
        """
        try:
            print(f"Finding neighboring genes for {gene_name}...")
            
            # Step 1: Get target gene location from Ensembl
            ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(ensembl_url, headers=headers, timeout=30)
            if response.status_code != 200:
                return f"Error: Could not find gene {gene_name} in Ensembl"
            
            target_gene = response.json()
            chrom = f"chr{target_gene['seq_region_name']}"
            target_start = target_gene['start']
            target_end = target_gene['end']
            
            print(f"Target gene {gene_name}: {chrom}:{target_start}-{target_end}")
            
            # Step 2: Query UCSC for all genes on this chromosome
            # Get a wider region to find neighbors (500kb each side)
            region_start = max(1, target_start - 500000)  # 500kb upstream
            region_end = target_end + 500000  # 500kb downstream
            
            ucsc_url = f"https://api.genome.ucsc.edu/getData/track"
            params = {
                'genome': genome,
                'track': track,
                'chrom': chrom,
                'start': region_start,
                'end': region_end
            }
            
            print(f"Querying UCSC for genes in {chrom}:{region_start}-{region_end}...")
            response = requests.get(ucsc_url, params=params, timeout=60)
            
            if response.status_code != 200:
                return f"Error: UCSC API returned status {response.status_code}"
            
            # Check if response has content
            if not response.text:
                return f"Error: UCSC API returned empty response for {chrom}:{region_start}-{region_end}"
            
            data = response.json()
            
            # Extract genes from response
            genes = data.get(track, [])
            if isinstance(genes, dict):
                # If returned as dict of chromosomes
                genes = genes.get(chrom, [])
            
            print(f"Found {len(genes)} genes in region")
            
            # Step 3: Find protein-coding genes to the left and right
            left_genes = {}  # Use dict to deduplicate by gene name
            right_genes = {}
            
            for gene in genes:
                gene_start = gene.get('chromStart', gene.get('txStart', 0))
                gene_end = gene.get('chromEnd', gene.get('txEnd', 0))
                # Prefer geneName (gene symbol) over name2 (transcript ID)
                gene_symbol = gene.get('geneName', gene.get('name2', gene.get('name', 'Unknown')))
                gene_type = gene.get('geneType', gene.get('type', 'unknown'))
                
                # Skip non-protein-coding genes
                if gene_type != 'protein_coding':
                    continue
                
                # Skip the target gene itself
                if gene_symbol.upper() == gene_name.upper():
                    continue
                
                # Genes to the left (end before target starts)
                if gene_end < target_start:
                    distance = target_start - gene_end
                    # Keep only the closest transcript for each gene
                    if gene_symbol not in left_genes or distance < left_genes[gene_symbol]['distance']:
                        left_genes[gene_symbol] = {
                            'name': gene_symbol,
                            'start': gene_start,
                            'end': gene_end,
                            'distance': distance
                        }
                
                # Genes to the right (start after target ends)
                elif gene_start > target_end:
                    distance = gene_start - target_end
                    # Keep only the closest transcript for each gene
                    if gene_symbol not in right_genes or distance < right_genes[gene_symbol]['distance']:
                        right_genes[gene_symbol] = {
                            'name': gene_symbol,
                            'start': gene_start,
                            'end': gene_end,
                            'distance': distance
                        }
            
            # Convert to lists and sort by distance
            left_genes = sorted(left_genes.values(), key=lambda x: x['distance'])
            right_genes = sorted(right_genes.values(), key=lambda x: x['distance'])
            
            # Format output
            output = f"Neighboring Protein-Coding Genes for {gene_name}\n"
            output += "=" * 70 + "\n\n"
            output += f"Target Gene: {gene_name}\n"
            output += f"Location: {chrom}:{target_start:,}-{target_end:,}\n"
            output += f"Genome: {genome}, Track: {track}\n\n"
            
            output += "Protein-Coding Genes to the LEFT (upstream):\n"
            output += "-" * 70 + "\n"
            if left_genes:
                for i, gene in enumerate(left_genes[:5], 1):  # Show top 5
                    output += f"{i}. {gene['name']}\n"
                    output += f"   Distance: ~{gene['distance']:,} nucleotides\n"
                    output += f"   Position: {chrom}:{gene['start']:,}-{gene['end']:,}\n\n"
            else:
                output += "No protein-coding genes found within 500kb upstream\n\n"
            
            output += "Protein-Coding Genes to the RIGHT (downstream):\n"
            output += "-" * 70 + "\n"
            if right_genes:
                for i, gene in enumerate(right_genes[:5], 1):  # Show top 5
                    output += f"{i}. {gene['name']}\n"
                    output += f"   Distance: ~{gene['distance']:,} nucleotides\n"
                    output += f"   Position: {chrom}:{gene['start']:,}-{gene['end']:,}\n\n"
            else:
                output += "No protein-coding genes found within 500kb downstream\n\n"
            
            output += f"\nNote: Showing only protein-coding genes (light blue in UCSC Browser).\n"
            output += f"Distance is measured from gene boundaries (not transcription start sites).\n"
            output += f"Search window: ±500kb from target gene.\n"
            
            return output
            
        except Exception as e:
            return f"Error finding neighboring genes: {str(e)}\n{traceback.format_exc()}"


    def gtex_expression(self, gene_name: str, top_n: int = 10) -> Dict:
        """Get gene expression data across tissues from GTEx Portal
        
        Note: GTEx API v2 has limited public access. This function provides
        gene information and direct links to GTEx Portal for visualization.
        
        Args:
            gene_name: Gene symbol (e.g., 'SOCS3', 'TP53')
            top_n: Number of top tissues to return (default: 10)
        
        Returns:
            Dict with success status and GTEx information
        """
        try:
            # Resolve gene name to Ensembl ID using Ensembl API
            ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(ensembl_url, headers=headers, timeout=30)
            
            if response.status_code != 200:
                return {
                    'success': False,
                    'error': f"Could not find Ensembl ID for gene {gene_name}"
                }
            
            ensembl_data = response.json()
            ensembl_id = ensembl_data.get('id')
            gene_description = ensembl_data.get('description', 'No description available')
            
            if not ensembl_id:
                return {
                    'success': False,
                    'error': f"No Ensembl ID found for {gene_name}"
                }
            
            # Build GTEx Portal URL for direct visualization
            gtex_gene_url = f"https://gtexportal.org/home/gene/{gene_name}"
            
            # Format output with gene information and GTEx link
            output = f"GTEx Expression Query for {gene_name}\n"
            output += "=" * 70 + "\n\n"
            output += f"Gene: {gene_name}\n"
            output += f"Ensembl ID: {ensembl_id}\n"
            output += f"Description: {gene_description}\n\n"
            
            output += "GTEx Portal Access:\n"
            output += "-" * 70 + "\n"
            output += f"View expression data at: {gtex_gene_url}\n\n"
            
            output += "The GTEx Portal provides:\n"
            output += "  • Median expression (TPM) across 54 tissues\n"
            output += "  • Interactive bar charts and heatmaps\n"
            output += "  • Expression by age and sex\n"
            output += "  • eQTL data (expression quantitative trait loci)\n"
            output += "  • Isoform-level expression\n\n"
            
            output += "To answer 'Which tissue has highest expression?':\n"
            output += "  1. Click the link above\n"
            output += "  2. Look at the 'Gene Expression' bar chart\n"
            output += "  3. The tallest bar shows the tissue with highest median TPM\n\n"
            
            output += "Note: GTEx v2 API has limited public query access.\n"
            output += "For programmatic access, consider the GTEx Portal bulk download:\n"
            output += "https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression\n"
            
            return {
                'success': True,
                'output': output,
                'ensembl_id': ensembl_id,
                'gene_name': gene_name,
                'gtex_url': gtex_gene_url
            }
        
        except Exception as e:
            return {
                'success': False,
                'error': f"GTEx query failed: {str(e)}"
            }

    def ucsc_table_query(self, genome: str, track: str, table: str = None, 
                        chrom: str = None, start: int = None, end: int = None,
                        filter_field: str = None, filter_value: str = None) -> str:
        """Query UCSC Table Browser using REST API
        
        Args:
            genome: Genome assembly (e.g., 'hg38', 'mm10')
            track: Track name (e.g., 'tRNAs', 'cons100way')
            table: Table name (optional, not used with REST API)
            chrom: Chromosome (e.g., 'chr16', 'chr6')
            start: Start position (0-based)
            end: End position (1-based)
            filter_field: Field name to filter on (e.g., 'score')
            filter_value: Filter value (e.g., '>=65')
        
        Returns:
            str: Query results or count of features
        """
        try:
            import requests
            
            # Map common track names to UCSC REST API track names
            # Based on UCSC documentation for supported track types
            track_map = {
                'trna': 'tRNAs',
                'trnas': 'tRNAs',
                'trnascan': 'tRNAs',
                'trnascan-se': 'tRNAs',
                'trnascan se': 'tRNAs',
                'trna genes': 'tRNAs',
                'trna_genes': 'tRNAs'
            }
            
            # Normalize track name
            track_normalized = track_map.get(track.lower(), track)
            
            # Special handling for score field filters - map to actual field names
            filter_field_map = {
                'score': 'trnaScore',  # The tRNAs track uses 'trnaScore' field
                'trnascore': 'trnaScore'
            }
            
            if filter_field:
                filter_field = filter_field_map.get(filter_field.lower(), filter_field)
            
            # Use UCSC REST API endpoint
            api_url = "https://api.genome.ucsc.edu/getData/track"
            
            # Build parameters according to REST API documentation
            params = {
                'genome': genome,
                'track': track_normalized,
                'maxItemsOutput': '1000000'  # Get all results
            }
            
            # Add chromosome and coordinates if specified
            if chrom:
                params['chrom'] = chrom
            if start is not None and end is not None:
                params['start'] = str(start)
                params['end'] = str(end)
            
            print(f"Querying UCSC REST API: {genome}, track={track_normalized}")
            if chrom:
                if start is not None:
                    print(f"  Region: {chrom}:{start}-{end}")
                else:
                    print(f"  Chromosome: {chrom}")
            
            response = requests.get(api_url, params=params, timeout=60)
            
            if response.status_code != 200:
                return f"Error: UCSC API returned status {response.status_code}\nTried track: {track_normalized}\nURL: {response.url}\nResponse: {response.text[:500]}"
            
            # Parse JSON response
            try:
                data = response.json()
            except:
                return f"Error: Could not parse JSON response from UCSC API\nResponse: {response.text[:500]}"
            
            # Check for errors in response
            if isinstance(data, dict) and 'error' in data:
                return f"Error from UCSC: {data['error']}"
            
            # Extract features from the response
            # The API returns data in different formats depending on track type
            features = []
            if isinstance(data, dict):
                # For tRNAs track, the data is under the track name key
                if track_normalized in data and isinstance(data[track_normalized], list):
                    features = data[track_normalized]
                # Look for other common keys that contain the actual data
                elif 'data' in data and isinstance(data['data'], list):
                    features = data['data']
                elif chrom in data and isinstance(data[chrom], list):
                    features = data[chrom]
                else:
                    # Some responses have the data directly as a single item
                    if 'chrom' in data or 'chromStart' in data:
                        features = [data]
            elif isinstance(data, list):
                features = data
            
            # Apply filters if specified
            if filter_field and filter_value and len(features) > 0:
                filtered = []
                for feature in features:
                    if isinstance(feature, dict):
                        field_val = feature.get(filter_field)
                        if field_val is not None:
                            try:
                                # Parse filter (e.g., '>=65')
                                if filter_value.startswith('>='):
                                    threshold = float(filter_value[2:])
                                    if float(field_val) >= threshold:
                                        filtered.append(feature)
                                elif filter_value.startswith('<='):
                                    threshold = float(filter_value[2:])
                                    if float(field_val) <= threshold:
                                        filtered.append(feature)
                                elif filter_value.startswith('='):
                                    if str(field_val) == filter_value[1:]:
                                        filtered.append(feature)
                            except (ValueError, TypeError):
                                continue
                
                features = filtered
                print(f"  After filtering {filter_field}{filter_value}: {len(features)} features")
            
            # Format output
            output = []
            output.append(f"UCSC REST API Query: {genome}")
            output.append(f"Track: {track_normalized}")
            if chrom:
                if start is not None:
                    region_str = f"{chrom}:{start}-{end}"
                else:
                    region_str = chrom
                output.append(f"Region: {region_str}")
            output.append(f"Total features found: {len(features)}")
            output.append("=" * 60)
            
            # Show first few features
            if len(features) > 0 and isinstance(features[0], dict):
                first_feat = features[0]
                name = first_feat.get('name', first_feat.get('id', 'unnamed'))
                output.append(f"\nFirst feature: {name}")
                
                if len(features) > 5:
                    output.append(f"\nShowing first 5 of {len(features)} features:")
                    for i, feat in enumerate(features[:5], 1):
                        name = feat.get('name', feat.get('id', 'unnamed'))
                        chrom_val = feat.get('chrom', '?')
                        start_pos = feat.get('chromStart', feat.get('start', '?'))
                        end_pos = feat.get('chromEnd', feat.get('end', '?'))
                        output.append(f"  {i}. {name} ({chrom_val}:{start_pos}-{end_pos})")
            
            return "\n".join(output)
            
        except Exception as e:
            import traceback
            return f"UCSC Table Browser query failed: {str(e)}\n{traceback.format_exc()}"

    def pubmed_search(self, query: str, max_results: int = 10, year: int = None) -> str:
        """Search PubMed using NCBI E-utilities
        
        Args:
            query: Search query (e.g., 'gas5[Title]', 'CRISPR cancer')
            max_results: Maximum number of results to return
            year: Limit to specific publication year
        
        Returns:
            str: Formatted search results with PMIDs, titles, and abstracts
        """
        try:
            import requests
            from xml.etree import ElementTree as ET
            
            # Build search query
            search_query = query
            if year:
                search_query += f" AND {year}[pdat]"
            
            print(f"Searching PubMed: {search_query}")
            
            # Step 1: Search for PMIDs
            esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                'db': 'pubmed',
                'term': search_query,
                'retmax': max_results,
                'retmode': 'json'
            }
            
            search_response = requests.get(esearch_url, params=search_params, timeout=30)
            if search_response.status_code != 200:
                return f"Error: PubMed search failed with status {search_response.status_code}"
            
            search_data = search_response.json()
            pmids = search_data.get('esearchresult', {}).get('idlist', [])
            
            if not pmids:
                return f"No PubMed articles found for: {search_query}"
            
            print(f"Found {len(pmids)} articles")
            
            # Step 2: Fetch article details
            efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            fetch_params = {
                'db': 'pubmed',
                'id': ','.join(pmids),
                'retmode': 'xml'
            }
            
            fetch_response = requests.get(efetch_url, params=fetch_params, timeout=30)
            if fetch_response.status_code != 200:
                return f"Error: Could not fetch article details"
            
            # Parse XML
            root = ET.fromstring(fetch_response.content)
            
            # Format output
            output = []
            output.append(f"PubMed Search Results: {search_query}")
            output.append(f"Found {len(pmids)} articles")
            output.append("=" * 60)
            
            for i, article in enumerate(root.findall('.//PubmedArticle'), 1):
                pmid_elem = article.find('.//PMID')
                pmid = pmid_elem.text if pmid_elem is not None else 'Unknown'
                
                title_elem = article.find('.//ArticleTitle')
                title = title_elem.text if title_elem is not None else 'No title'
                
                # Get authors
                authors = []
                for author in article.findall('.//Author')[:3]:  # First 3 authors
                    lastname = author.find('LastName')
                    if lastname is not None:
                        authors.append(lastname.text)
                author_str = ', '.join(authors)
                if len(article.findall('.//Author')) > 3:
                    author_str += ', et al.'
                
                # Get year
                year_elem = article.find('.//PubDate/Year')
                pub_year = year_elem.text if year_elem is not None else 'Unknown'
                
                # Get abstract
                abstract_elem = article.find('.//AbstractText')
                abstract = abstract_elem.text if abstract_elem is not None else 'No abstract available'
                if abstract and len(abstract) > 500:
                    abstract = abstract[:500] + '...'
                
                output.append(f"\n{i}. PMID: {pmid}")
                output.append(f"   Title: {title}")
                output.append(f"   Authors: {author_str}")
                output.append(f"   Year: {pub_year}")
                output.append(f"   Abstract: {abstract}")
                output.append("")
            
            return "\n".join(output)
            
        except Exception as e:
            import traceback
            return f"PubMed search failed: {str(e)}\n{traceback.format_exc()}"

    def track_intersection(self, genome: str, track1: str, track2: str, 
                          chrom: str = None, max_results: int = 10) -> str:
        """Find genomic features that overlap between two tracks
        
        Args:
            genome: Genome assembly (e.g., 'hg38')
            track1: First track name (e.g., 'tRNAs')
            track2: Second track name to intersect with (e.g., 'knownGene')
            chrom: Optional chromosome to limit search
            max_results: Maximum number of overlapping features to show
        
        Returns:
            str: List of overlapping features with details
        """
        try:
            import requests
            
            # Map common track names
            track_map = {
                'trna': 'tRNAs',
                'trnas': 'tRNAs',
                'trnascan': 'tRNAs',
                'gencode': 'knownGene',
                'gencode v48': 'knownGene',
                'all gencode v48': 'knownGene',
                'genes': 'knownGene'
            }
            
            track1_norm = track_map.get(track1.lower(), track1)
            track2_norm = track_map.get(track2.lower(), track2)
            
            api_url = "https://api.genome.ucsc.edu/getData/track"
            
            print(f"Finding overlaps between {track1_norm} and {track2_norm} in {genome}")
            
            # Get data from both tracks
            # Limit to specific chromosome if provided, otherwise get all
            params1 = {
                'genome': genome,
                'track': track1_norm,
                'maxItemsOutput': '1000000'
            }
            if chrom:
                params1['chrom'] = chrom
            
            print(f"Fetching {track1_norm} data...")
            response1 = requests.get(api_url, params=params1, timeout=120)
            if response1.status_code != 200:
                return f"Error fetching {track1_norm}: Status {response1.status_code}\nURL: {response1.url}"
            
            data1 = response1.json()
            features1 = []
            if isinstance(data1, dict) and track1_norm in data1:
                track_data = data1[track1_norm]
                # Handle both list and dict-of-lists formats
                if isinstance(track_data, list):
                    features1 = track_data
                elif isinstance(track_data, dict):
                    # Dict with chromosome keys, flatten all lists
                    for chrom_key, chrom_features in track_data.items():
                        if isinstance(chrom_features, list):
                            features1.extend(chrom_features)
            
            if not features1:
                return f"No features found in {track1_norm}"
            
            print(f"Found {len(features1)} features in {track1_norm}")
            
            # Get track2 data
            params2 = {
                'genome': genome,
                'track': track2_norm,
                'maxItemsOutput': '1000000'
            }
            if chrom:
                params2['chrom'] = chrom
            
            print(f"Fetching {track2_norm} data...")
            response2 = requests.get(api_url, params=params2, timeout=120)
            if response2.status_code != 200:
                return f"Error fetching {track2_norm}: Status {response2.status_code}"
            
            data2 = response2.json()
            features2 = []
            if isinstance(data2, dict) and track2_norm in data2:
                track_data = data2[track2_norm]
                # Handle both list and dict-of-lists formats
                if isinstance(track_data, list):
                    features2 = track_data
                elif isinstance(track_data, dict):
                    # Dict with chromosome keys, flatten all lists
                    for chrom_key, chrom_features in track_data.items():
                        if isinstance(chrom_features, list):
                            features2.extend(chrom_features)
            
            if not features2:
                return f"No features found in {track2_norm}"
            
            print(f"Found {len(features2)} features in {track2_norm}")
            
            # Find overlaps
            print("Computing intersections...")
            overlaps = []
            
            for feat1 in features1:
                # Skip if not a dictionary
                if not isinstance(feat1, dict):
                    continue
                    
                chrom1 = feat1.get('chrom')
                start1 = feat1.get('chromStart', feat1.get('start'))
                end1 = feat1.get('chromEnd', feat1.get('end'))
                name1 = feat1.get('name', 'unnamed')
                
                if start1 is None or end1 is None:
                    continue
                
                # Check for overlaps with track2 features
                for feat2 in features2:
                    # Skip if not a dictionary
                    if not isinstance(feat2, dict):
                        continue
                        
                    chrom2 = feat2.get('chrom')
                    start2 = feat2.get('chromStart', feat2.get('start'))
                    end2 = feat2.get('chromEnd', feat2.get('end'))
                    
                    if start2 is None or end2 is None or chrom1 != chrom2:
                        continue
                    
                    # Check for overlap: features overlap if they share any base
                    if start1 < end2 and start2 < end1:
                        name2 = feat2.get('name', feat2.get('name2', 'unnamed'))
                        overlaps.append({
                            'feature1': name1,
                            'chrom': chrom1,
                            'start1': start1,
                            'end1': end1,
                            'feature2': name2,
                            'start2': start2,
                            'end2': end2
                        })
                        break  # Found overlap for this feat1, move to next
            
            print(f"Found {len(overlaps)} overlapping features")
            
            # Sort overlaps by genomic position (chromosome then start coordinate)
            def sort_key(overlap):
                chrom = overlap['chrom']
                # Extract numeric part of chromosome for proper sorting (chr1, chr2, ..., chr10, chr11, etc.)
                if chrom.startswith('chr'):
                    chrom_num = chrom[3:]
                    # Handle chrX, chrY, chrM specially
                    if chrom_num == 'X':
                        return (23, overlap['start1'])
                    elif chrom_num == 'Y':
                        return (24, overlap['start1'])
                    elif chrom_num == 'M':
                        return (25, overlap['start1'])
                    else:
                        try:
                            return (int(chrom_num), overlap['start1'])
                        except ValueError:
                            return (26, overlap['start1'])  # Other chromosomes
                return (27, overlap['start1'])  # Non-chr format
            
            overlaps.sort(key=sort_key)
            
            # Format output
            output = []
            output.append(f"Track Intersection Analysis: {genome}")
            output.append(f"Track 1: {track1_norm} ({len(features1)} features)")
            output.append(f"Track 2: {track2_norm} ({len(features2)} features)")
            if chrom:
                output.append(f"Region: {chrom}")
            output.append(f"Total overlaps found: {len(overlaps)}")
            output.append("=" * 60)
            
            if overlaps:
                output.append(f"\nFirst overlapping feature: {overlaps[0]['feature1']}")
                output.append(f"  Location: {overlaps[0]['chrom']}:{overlaps[0]['start1']}-{overlaps[0]['end1']}")
                output.append(f"  Overlaps with: {overlaps[0]['feature2']}")
                output.append(f"  Location: {overlaps[0]['chrom']}:{overlaps[0]['start2']}-{overlaps[0]['end2']}")
                
                if len(overlaps) > 1:
                    show_count = min(max_results, len(overlaps))
                    output.append(f"\nShowing first {show_count} of {len(overlaps)} overlaps:")
                    for i, overlap in enumerate(overlaps[:show_count], 1):
                        output.append(f"\n{i}. {overlap['feature1']} ({overlap['chrom']}:{overlap['start1']}-{overlap['end1']})")
                        output.append(f"   ∩ {overlap['feature2']} ({overlap['chrom']}:{overlap['start2']}-{overlap['end2']})")
            else:
                output.append("\nNo overlapping features found.")
            
            return "\n".join(output)
            
        except Exception as e:
            import traceback
            return f"Track intersection failed: {str(e)}\n{traceback.format_exc()}"
    
    # ============================================================
    # BIOPYTHON-POWERED TOOLS (Complement EMBOSS tools)
    # ============================================================
    
    def parse_fasta_biopython(self, fasta_text: str) -> List[Tuple[str, str]]:
        """Parse FASTA format using Biopython (more robust than regex)
        
        Args:
            fasta_text: FASTA formatted text
        
        Returns:
            List of (header, sequence) tuples
        """
        try:
            from io import StringIO
            from Bio import SeqIO
            
            sequences = []
            for record in SeqIO.parse(StringIO(fasta_text), "fasta"):
                sequences.append((record.description, str(record.seq)))
            
            return sequences
        except Exception as e:
            return [(f"Error: {str(e)}", "")]
    
    def multiple_sequence_alignment(self, sequences: List[Tuple[str, str]]) -> str:
        """Perform multiple sequence alignment using Biopython
        
        Args:
            sequences: List of (name, sequence) tuples
        
        Returns:
            str: Alignment visualization
        """
        try:
            from Bio import Align
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            
            if len(sequences) < 2:
                return "Error: Need at least 2 sequences for alignment"
            
            # Create sequence records
            seq_records = [SeqRecord(Seq(seq), id=name[:20]) for name, seq in sequences]
            
            # Simple pairwise alignment for visualization
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            output = ["Multiple Sequence Alignment", "=" * 50, ""]
            
            # Align first two sequences as demo
            alignments = aligner.align(seq_records[0].seq, seq_records[1].seq)
            alignment = alignments[0]
            
            output.append(f"Alignment Score: {alignment.score}")
            output.append(f"\n{seq_records[0].id}:")
            output.append(str(alignment[0]))
            output.append(f"\n{seq_records[1].id}:")
            output.append(str(alignment[1]))
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Multiple alignment failed: {str(e)}"
    
    def get_protein_structure_info(self, pdb_id: str) -> str:
        """Get protein structure information from PDB using Biopython
        
        Args:
            pdb_id: PDB identifier (e.g., '1ABC')
        
        Returns:
            str: Structure information
        """
        try:
            from Bio.PDB import PDBList, PDBParser
            import tempfile
            import shutil
            
            # Create temp directory
            temp_dir = tempfile.mkdtemp()
            
            try:
                # Download PDB file
                pdbl = PDBList()
                pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=temp_dir, file_format='pdb')
                
                # Parse structure
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure(pdb_id, pdb_file)
                
                output = [f"PDB Structure: {pdb_id}", "=" * 50, ""]
                
                # Get models
                models = list(structure.get_models())
                output.append(f"Number of models: {len(models)}")
                
                # Get chains
                for model in models[:1]:  # Just first model
                    chains = list(model.get_chains())
                    output.append(f"\nChains: {len(chains)}")
                    
                    for chain in chains:
                        residues = list(chain.get_residues())
                        output.append(f"  Chain {chain.id}: {len(residues)} residues")
                        
                        # Count atom types
                        atoms = list(chain.get_atoms())
                        output.append(f"    Total atoms: {len(atoms)}")
                
                return "\n".join(output)
                
            finally:
                # Cleanup temp directory
                shutil.rmtree(temp_dir, ignore_errors=True)
                
        except Exception as e:
            return f"PDB structure retrieval failed: {str(e)}\nTry a valid PDB ID like '1ABC' or '3J7Y'"
    
    def calculate_molecular_weight_biopython(self, sequence: str, seq_type: str = 'protein') -> str:
        """Calculate molecular weight using Biopython (alternative to pepstats)
        
        Args:
            sequence: Amino acid or DNA sequence
            seq_type: 'protein' or 'dna'
        
        Returns:
            str: Molecular weight information
        """
        try:
            from Bio.SeqUtils import molecular_weight
            from Bio.Seq import Seq
            
            seq = Seq(sequence.upper())
            
            if seq_type == 'protein':
                mw = molecular_weight(seq, seq_type='protein')
            else:
                mw = molecular_weight(seq, seq_type='DNA')
            
            output = [
                f"Molecular Weight Analysis ({seq_type.upper()})",
                "=" * 50,
                f"Sequence length: {len(sequence)}",
                f"Molecular weight: {mw:.2f} Da",
                f"Molecular weight: {mw/1000:.2f} kDa",
                ""
            ]
            
            if seq_type == 'protein':
                # Additional protein info
                from Bio.SeqUtils.ProtParam import ProteinAnalysis
                protein = ProteinAnalysis(str(seq))
                
                output.append(f"Isoelectric point: {protein.isoelectric_point():.2f}")
                output.append(f"Aromaticity: {protein.aromaticity():.3f}")
                output.append(f"Instability index: {protein.instability_index():.2f}")
                
                # Amino acid composition
                aa_comp = protein.get_amino_acids_percent()
                output.append("\nAmino Acid Composition (%):")
                for aa in sorted(aa_comp.keys()):
                    if aa_comp[aa] > 0:
                        output.append(f"  {aa}: {aa_comp[aa]*100:.1f}%")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Molecular weight calculation failed: {str(e)}"
    
    def phylogenetic_tree(self, sequences: list) -> str:
        """Build phylogenetic tree from sequences using UPGMA clustering
        
        Args:
            sequences: List of (name, sequence) tuples
            
        Returns:
            str: Tree in Newick format with distance matrix
        """
        try:
            from Bio import Phylo
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            from Bio.Align import MultipleSeqAlignment
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            import io
            
            if len(sequences) < 3:
                return "Error: Need at least 3 sequences for phylogenetic tree"
            
            # Create alignment
            seq_records = [SeqRecord(Seq(seq), id=name) for name, seq in sequences]
            alignment = MultipleSeqAlignment(seq_records)
            
            # Calculate distance matrix
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(alignment)
            
            # Construct tree using UPGMA
            constructor = DistanceTreeConstructor(calculator, 'upgma')
            tree = constructor.build_tree(alignment)
            
            # Convert tree to string
            output = io.StringIO()
            Phylo.draw_ascii(tree, file=output)
            tree_ascii = output.getvalue()
            
            # Also get Newick format
            newick_output = io.StringIO()
            Phylo.write(tree, newick_output, 'newick')
            newick = newick_output.getvalue()
            
            result = ["Phylogenetic Tree (UPGMA):\n", tree_ascii]
            result.append(f"\nNewick Format:\n{newick}")
            result.append("\nDistance Matrix:")
            for i, name1 in enumerate(dm.names):
                row = [f"{dm[name1, name2]:.3f}" for name2 in dm.names]
                result.append(f"  {name1}: {' '.join(row)}")
            
            return "\n".join(result)
            
        except Exception as e:
            return f"Phylogenetic tree construction failed: {str(e)}"
    
    def find_motifs(self, sequences: list) -> str:
        """Find common motifs in sequences
        
        Args:
            sequences: List of (name, sequence) tuples
            
        Returns:
            str: Motif analysis results
        """
        try:
            from Bio import motifs
            from Bio.Seq import Seq
            
            if len(sequences) < 2:
                return "Error: Need at least 2 sequences for motif finding"
            
            # Create motif instances
            instances = [Seq(seq) for _, seq in sequences]
            m = motifs.create(instances)
            
            output = ["Motif Analysis:\n"]
            output.append(f"Number of sequences: {len(sequences)}")
            output.append(f"Sequence length: {len(instances[0])}")
            output.append(f"\nConsensus sequence: {m.consensus}")
            output.append(f"Anticonsensus sequence: {m.anticonsensus}")
            
            # Position weight matrix
            output.append("\nPosition Weight Matrix:")
            pwm = m.counts.normalize()
            for base in ['A', 'C', 'G', 'T']:
                if base in pwm:
                    values = [f"{v:.2f}" for v in pwm[base]]
                    output.append(f"  {base}: {' '.join(values[:20])}")  # Show first 20 positions
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Motif finding failed: {str(e)}"
    
    def motif_consensus(self, sequences: list) -> str:
        """Generate consensus motif from sequences
        
        Args:
            sequences: List of (name, sequence) tuples
            
        Returns:
            str: Consensus sequence with IUPAC codes
        """
        try:
            from Bio import motifs
            from Bio.Seq import Seq
            
            instances = [Seq(seq) for _, seq in sequences]
            m = motifs.create(instances)
            
            output = [f"Consensus: {m.consensus}"]
            output.append(f"Degenerate consensus (IUPAC): {m.degenerate_consensus}")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Consensus motif generation failed: {str(e)}"
    
    def primer_analysis(self, sequence: str) -> str:
        """Analyze primer properties (Tm, GC content, secondary structure)
        
        Args:
            sequence: Primer sequence
            
        Returns:
            str: Primer analysis results
        """
        try:
            from Bio.SeqUtils import MeltingTemp as mt
            from Bio.Seq import Seq
            
            seq = Seq(sequence.upper().replace(' ', '').replace('\n', ''))
            
            if len(seq) < 10:
                return "Warning: Primer too short (< 10 bp)"
            if len(seq) > 40:
                return "Warning: Primer too long (> 40 bp)"
            
            # Calculate melting temperature using different methods
            tm_wallace = mt.Tm_Wallace(seq)
            tm_gc = mt.Tm_GC(seq)
            
            # GC content
            gc_count = seq.count('G') + seq.count('C')
            gc_percent = (gc_count / len(seq)) * 100
            
            # Check for runs
            has_runs = any(base*4 in str(seq) for base in 'ATGC')
            
            # GC clamp
            gc_clamp = str(seq)[-5:].count('G') + str(seq)[-5:].count('C')
            
            output = [f"Primer Analysis for: {sequence}\n"]
            output.append(f"Length: {len(seq)} bp")
            output.append(f"GC Content: {gc_percent:.1f}%")
            output.append(f"Tm (Wallace): {tm_wallace:.1f}°C")
            output.append(f"Tm (GC-based): {tm_gc:.1f}°C")
            output.append(f"GC Clamp (last 5 bp): {gc_clamp} G/C bases")
            
            # Quality checks
            output.append("\nQuality Checks:")
            if 40 <= gc_percent <= 60:
                output.append("  ✓ GC content optimal (40-60%)")
            else:
                output.append(f"  ✗ GC content suboptimal ({gc_percent:.1f}%)")
            
            if 50 <= tm_wallace <= 65:
                output.append("  ✓ Tm in optimal range (50-65°C)")
            else:
                output.append(f"  ✗ Tm outside optimal range ({tm_wallace:.1f}°C)")
            
            if gc_clamp >= 2:
                output.append("  ✓ Good GC clamp (≥2)")
            else:
                output.append("  ✗ Weak GC clamp")
            
            if has_runs:
                output.append("  ✗ Contains runs of 4+ identical bases")
            else:
                output.append("  ✓ No problematic runs")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Primer analysis failed: {str(e)}"
    
    def restriction_batch_analysis(self, sequence: str, enzymes=None) -> str:
        """Analyze restriction enzyme sites in batch
        
        Args:
            sequence: DNA sequence
            enzymes: List of enzyme names (None = all commercial enzymes)
            
        Returns:
            str: Restriction analysis results
        """
        try:
            from Bio import Restriction
            from Bio.Seq import Seq
            
            seq = Seq(sequence.upper().replace(' ', '').replace('\n', ''))
            
            if enzymes:
                # Use specified enzymes
                enzyme_list = []
                for enz_name in enzymes:
                    try:
                        enzyme_list.append(getattr(Restriction, enz_name))
                    except AttributeError:
                        continue
                rb = Restriction.RestrictionBatch(enzyme_list)
            else:
                # Use all commercial enzymes
                rb = Restriction.CommOnly
            
            # Analyze sequence
            analysis = rb.search(seq)
            
            # Filter enzymes that cut
            cutters = {enz: sites for enz, sites in analysis.items() if sites}
            
            output = [f"Restriction Analysis ({len(seq)} bp)\n"]
            output.append(f"Enzymes tested: {len(analysis)}")
            output.append(f"Enzymes that cut: {len(cutters)}\n")
            
            if cutters:
                output.append("Cutting Enzymes:")
                for enz, sites in sorted(cutters.items(), key=lambda x: len(x[1]), reverse=True)[:20]:
                    site_list = ', '.join(str(s) for s in sorted(sites)[:10])
                    if len(sites) > 10:
                        site_list += f"... (+{len(sites)-10} more)"
                    output.append(f"  {enz}: {len(sites)} site(s) at {site_list}")
            else:
                output.append("No restriction sites found")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Restriction analysis failed: {str(e)}"
    
    def restriction_map(self, sequence: str) -> str:
        """Generate restriction map showing all sites
        
        Args:
            sequence: DNA sequence
            
        Returns:
            str: Visual restriction map
        """
        try:
            from Bio import Restriction
            from Bio.Seq import Seq
            
            seq = Seq(sequence.upper().replace(' ', '').replace('\n', ''))
            rb = Restriction.CommOnly
            analysis = rb.search(seq)
            
            # Get all cutting positions
            all_sites = []
            for enz, sites in analysis.items():
                if sites:
                    for site in sites:
                        all_sites.append((site, str(enz)))
            
            all_sites.sort()
            
            output = [f"Restriction Map ({len(seq)} bp):\n"]
            
            if all_sites:
                # Show map
                map_scale = min(100, len(seq))
                output.append(f"Position  Enzyme(s)")
                output.append("-" * 40)
                
                for pos, enz in all_sites[:30]:  # Show first 30 sites
                    output.append(f"{pos:>8}  {enz}")
                
                if len(all_sites) > 30:
                    output.append(f"... +{len(all_sites)-30} more sites")
            else:
                output.append("No restriction sites found")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Restriction map generation failed: {str(e)}"
    
    def extract_sequence_features(self, sequence: str, seq_type: str = 'dna') -> str:
        """Extract sequence features and statistics
        
        Args:
            sequence: Input sequence
            seq_type: 'dna', 'rna', or 'protein'
            
        Returns:
            str: Feature analysis
        """
        try:
            from Bio.Seq import Seq
            from Bio.SeqUtils import gc_fraction
            
            seq = Seq(sequence.upper().replace(' ', '').replace('\n', ''))
            
            output = [f"Sequence Features ({seq_type.upper()}):\n"]
            output.append(f"Length: {len(seq)} bp")
            
            if seq_type.lower() in ['dna', 'rna']:
                # Nucleotide features
                gc = gc_fraction(seq) * 100
                output.append(f"GC Content: {gc:.2f}%")
                output.append(f"AT Content: {100-gc:.2f}%")
                
                # Count bases
                if seq_type.lower() == 'dna':
                    output.append(f"\nBase Composition:")
                    output.append(f"  A: {seq.count('A')} ({seq.count('A')/len(seq)*100:.1f}%)")
                    output.append(f"  T: {seq.count('T')} ({seq.count('T')/len(seq)*100:.1f}%)")
                    output.append(f"  G: {seq.count('G')} ({seq.count('G')/len(seq)*100:.1f}%)")
                    output.append(f"  C: {seq.count('C')} ({seq.count('C')/len(seq)*100:.1f}%)")
                
                # Find repeats
                output.append(f"\nDinucleotide Repeats:")
                for di in ['AA', 'TT', 'GG', 'CC', 'AT', 'TA', 'GC', 'CG']:
                    count = str(seq).count(di)
                    if count > 5:
                        output.append(f"  {di}: {count}")
            
            elif seq_type.lower() == 'protein':
                # Protein features
                from Bio.SeqUtils.ProtParam import ProteinAnalysis
                pa = ProteinAnalysis(str(seq))
                
                output.append(f"Molecular Weight: {pa.molecular_weight():.2f} Da")
                output.append(f"Isoelectric Point: {pa.isoelectric_point():.2f}")
                output.append(f"Aromaticity: {pa.aromaticity():.3f}")
                
                # Secondary structure prediction
                sec_struc = pa.secondary_structure_fraction()
                output.append(f"\nSecondary Structure Prediction:")
                output.append(f"  Helix: {sec_struc[0]*100:.1f}%")
                output.append(f"  Turn: {sec_struc[1]*100:.1f}%")
                output.append(f"  Sheet: {sec_struc[2]*100:.1f}%")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Feature extraction failed: {str(e)}"
    
    def transcribe_sequence(self, sequence: str) -> str:
        """Transcribe DNA to RNA
        
        Args:
            sequence: DNA sequence
            
        Returns:
            str: RNA sequence
        """
        try:
            from Bio.Seq import Seq
            
            seq = Seq(sequence.upper().replace(' ', '').replace('\n', ''))
            rna = seq.transcribe()
            
            return f"DNA: {sequence}\nRNA: {rna}"
            
        except Exception as e:
            return f"Transcription failed: {str(e)}"
    
    def back_transcribe(self, sequence: str) -> str:
        """Back-transcribe RNA to DNA
        
        Args:
            sequence: RNA sequence
            
        Returns:
            str: DNA sequence
        """
        try:
            from Bio.Seq import Seq
            
            seq = Seq(sequence.upper().replace(' ', '').replace('\n', ''))
            dna = seq.back_transcribe()
            
            return f"RNA: {sequence}\nDNA: {dna}"
            
        except Exception as e:
            return f"Back-transcription failed: {str(e)}"
    
    def codon_optimize(self, sequence: str, organism: str = 'human') -> str:
        """Optimize codon usage for expression (simplified)
        
        Args:
            sequence: Protein or DNA sequence
            organism: Target organism ('human', 'ecoli', 'yeast')
            
        Returns:
            str: Codon optimization analysis
        """
        try:
            from Bio.Seq import Seq
            from Bio.Data import CodonTable
            
            seq = Seq(sequence.upper().replace(' ', '').replace('\n', ''))
            
            # If DNA, translate first
            if all(c in 'ATGC' for c in str(seq)):
                protein = seq.translate()
                is_dna = True
            else:
                protein = seq
                is_dna = False
            
            # Get codon table
            table = CodonTable.unambiguous_dna_by_name["Standard"]
            
            output = [f"Codon Optimization Analysis (Target: {organism}):\n"]
            
            if is_dna:
                output.append(f"Original DNA: {sequence[:60]}...")
                output.append(f"Protein: {protein[:30]}...")
                
                # Analyze codon usage
                codons = [str(seq)[i:i+3] for i in range(0, len(seq)-2, 3)]
                output.append(f"\nCodon count: {len(codons)}")
                
                # Find rare codons (simplified - actual optimization needs codon usage tables)
                rare_codons = ['CTA', 'CGA', 'CGG', 'ATA', 'TTA']
                rare_count = sum(1 for c in codons if c in rare_codons)
                output.append(f"Potential rare codons: {rare_count}")
                
                if rare_count > 0:
                    output.append("\nNote: Consider optimizing rare codons for better expression")
            else:
                output.append(f"Protein: {protein}")
                output.append("\nNote: Provide DNA sequence for detailed codon optimization")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Codon optimization failed: {str(e)}"
    
    def pairwise_comparison(self, seq1: str, seq2: str) -> str:
        """Detailed pairwise sequence comparison with statistics
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            str: Comparison statistics
        """
        try:
            from Bio import pairwise2
            from Bio.Seq import Seq
            
            s1 = Seq(seq1.upper().replace(' ', '').replace('\n', ''))
            s2 = Seq(seq2.upper().replace(' ', '').replace('\n', ''))
            
            # Perform alignment
            alignments = pairwise2.align.globalxx(s1, s2, one_alignment_only=True)
            
            if not alignments:
                return "No alignment found"
            
            alignment = alignments[0]
            
            # Calculate identity
            matches = sum(1 for a, b in zip(alignment.seqA, alignment.seqB) if a == b and a != '-')
            identity = (matches / max(len(s1), len(s2))) * 100
            
            output = ["Pairwise Sequence Comparison:\n"]
            output.append(f"Sequence 1 length: {len(s1)}")
            output.append(f"Sequence 2 length: {len(s2)}")
            output.append(f"Alignment score: {alignment.score}")
            output.append(f"Identity: {identity:.2f}%")
            output.append(f"Matches: {matches}")
            
            # Show alignment (first 60 chars)
            output.append(f"\nAlignment preview:")
            output.append(f"Seq1: {alignment.seqA[:60]}")
            output.append(f"Seq2: {alignment.seqB[:60]}")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Pairwise comparison failed: {str(e)}"
    
    def predict_secondary_structure(self, sequence: str) -> str:
        """Predict protein secondary structure
        
        Args:
            sequence: Protein sequence
            
        Returns:
            str: Secondary structure prediction
        """
        try:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            
            seq = sequence.upper().replace(' ', '').replace('\n', '')
            pa = ProteinAnalysis(seq)
            
            # Get secondary structure fractions
            helix, turn, sheet = pa.secondary_structure_fraction()
            
            output = ["Secondary Structure Prediction:\n"]
            output.append(f"Sequence length: {len(seq)} aa")
            output.append(f"\nEstimated fractions:")
            output.append(f"  α-Helix: {helix*100:.1f}%")
            output.append(f"  β-Turn: {turn*100:.1f}%")
            output.append(f"  β-Sheet: {sheet*100:.1f}%")
            output.append(f"  Random coil: {(1-helix-turn-sheet)*100:.1f}%")
            
            output.append(f"\nNote: This is a simple statistical prediction.")
            output.append("For accurate predictions, use tools like DSSP or PSIPRED.")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Secondary structure prediction failed: {str(e)}"
    
    def hydrophobicity_plot(self, sequence: str) -> str:
        """Calculate hydrophobicity profile using Kyte-Doolittle scale
        
        Args:
            sequence: Protein sequence
            
        Returns:
            str: Hydrophobicity analysis
        """
        try:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            
            seq = sequence.upper().replace(' ', '').replace('\n', '')
            pa = ProteinAnalysis(seq)
            
            # Get hydrophobicity
            hydro = pa.protein_scale(window=9)
            
            if not hydro:
                return "Error: Sequence too short for hydrophobicity analysis (need >9 aa)"
            
            output = ["Hydrophobicity Analysis (Kyte-Doolittle):\n"]
            output.append(f"Sequence length: {len(seq)} aa")
            output.append(f"Window size: 9")
            output.append(f"\nAverage hydrophobicity: {sum(hydro)/len(hydro):.3f}")
            output.append(f"Maximum: {max(hydro):.3f}")
            output.append(f"Minimum: {min(hydro):.3f}")
            
            # Identify hydrophobic regions (positive values)
            hydrophobic_regions = [i for i, h in enumerate(hydro) if h > 1.0]
            if hydrophobic_regions:
                output.append(f"\nHydrophobic regions (score > 1.0): {len(hydrophobic_regions)} positions")
            
            # Show profile (sample every 10th position for brevity)
            output.append(f"\nHydrophobicity profile (sample):")
            for i in range(0, len(hydro), 10):
                output.append(f"  Position {i+5}: {hydro[i]:.2f}")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Hydrophobicity analysis failed: {str(e)}"
    
    def sequence_entropy(self, sequence: str) -> str:
        """Calculate Shannon entropy of sequence
        
        Args:
            sequence: DNA, RNA, or protein sequence
            
        Returns:
            str: Entropy analysis
        """
        try:
            import math
            from collections import Counter
            
            seq = sequence.upper().replace(' ', '').replace('\n', '')
            
            # Calculate frequency
            counts = Counter(seq)
            total = len(seq)
            
            # Calculate Shannon entropy
            entropy = 0
            for count in counts.values():
                p = count / total
                entropy -= p * math.log2(p)
            
            # Maximum possible entropy
            alphabet_size = len(counts)
            max_entropy = math.log2(alphabet_size)
            
            # Normalized entropy
            norm_entropy = entropy / max_entropy if max_entropy > 0 else 0
            
            output = ["Sequence Entropy Analysis:\n"]
            output.append(f"Sequence length: {len(seq)}")
            output.append(f"Alphabet size: {alphabet_size}")
            output.append(f"Shannon entropy: {entropy:.3f} bits")
            output.append(f"Maximum entropy: {max_entropy:.3f} bits")
            output.append(f"Normalized entropy: {norm_entropy:.3f}")
            
            output.append(f"\nSymbol frequencies:")
            for symbol, count in sorted(counts.items()):
                freq = count / total
                output.append(f"  {symbol}: {count} ({freq*100:.1f}%)")
            
            if norm_entropy > 0.9:
                output.append("\n✓ High complexity sequence (diverse)")
            elif norm_entropy < 0.5:
                output.append("\n⚠ Low complexity sequence (repetitive)")
            else:
                output.append("\n• Moderate complexity")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Entropy calculation failed: {str(e)}"
    
    def sequence_complexity(self, sequence: str) -> str:
        """Analyze sequence complexity using multiple metrics
        
        Args:
            sequence: Input sequence
            
        Returns:
            str: Complexity analysis
        """
        try:
            from collections import Counter
            import math
            
            seq = sequence.upper().replace(' ', '').replace('\n', '')
            
            # Calculate various complexity metrics
            
            # 1. Base composition evenness
            counts = Counter(seq)
            total = len(seq)
            
            # 2. Longest homopolymer run
            max_run = 1
            current_run = 1
            for i in range(1, len(seq)):
                if seq[i] == seq[i-1]:
                    current_run += 1
                    max_run = max(max_run, current_run)
                else:
                    current_run = 1
            
            # 3. Dinucleotide diversity
            dinucs = [seq[i:i+2] for i in range(len(seq)-1)]
            dinuc_diversity = len(set(dinucs))
            
            # 4. Shannon entropy
            entropy = 0
            for count in counts.values():
                p = count / total
                entropy -= p * math.log2(p)
            
            output = ["Sequence Complexity Analysis:\n"]
            output.append(f"Length: {len(seq)}")
            output.append(f"Unique symbols: {len(counts)}")
            output.append(f"Shannon entropy: {entropy:.3f} bits")
            output.append(f"Longest run: {max_run} ('{seq[0]}')")
            output.append(f"Dinucleotide diversity: {dinuc_diversity}")
            
            # Complexity assessment
            output.append("\nComplexity Assessment:")
            if max_run > 10:
                output.append(f"  ⚠ Long homopolymer run ({max_run} bp)")
            if entropy < 1.5:
                output.append("  ⚠ Low entropy (repetitive)")
            if dinuc_diversity < 10:
                output.append("  ⚠ Low dinucleotide diversity")
            
            if max_run <= 10 and entropy >= 1.5 and dinuc_diversity >= 10:
                output.append("  ✓ Good sequence complexity")
            
            return "\n".join(output)
            
        except Exception as e:
            return f"Complexity analysis failed: {str(e)}"


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
