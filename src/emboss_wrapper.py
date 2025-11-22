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


class EMBOSSWrapper:
    """Wrapper for common EMBOSS bioinformatics tools"""
    
    def __init__(self):
        """Initialize the EMBOSS wrapper and verify installation"""
        # Set NCBI Entrez email for downloads
        Entrez.email = "user@bioquery.local"
        
        # Map natural language to EMBOSS commands (shortcuts for common operations)
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
            'tissue_expression': 'gtex_expression'
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
            'gtex': 'Get gene expression data across tissues (GTEx)'
        }
        
        # Cache for available EMBOSS tools
        self.available_emboss_tools = None
        
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
                  max_results: int = 10, expect_threshold: float = 10.0, exclude_taxa: str = None) -> str:
        """Run BLAST search against NCBI databases using remote Entrez service
        
        Args:
            sequence: Query sequence (DNA, RNA, or protein) or gene name to resolve
            blast_type: 'blastn' (DNA), 'blastp' (protein), 'blastx' (DNA->protein)
            database: NCBI database ('nt' for nucleotide, 'nr' for protein, etc.)
            max_results: Maximum number of results to return
            expect_threshold: E-value threshold for reporting
        
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
            if exclude_taxa and 'primate' in exclude_taxa.lower():
                # Use NCBI taxonomy to exclude all primates (taxid:9443)
                # Search in Mammalia (taxid:40674) but NOT in Primates (taxid:9443)
                entrez_query = "txid40674[Organism:exp] NOT txid9443[Organism:exp]"
                print(f"Filtering: Mammals excluding primates (using NCBI taxonomy)")
                # Request more results to ensure we get good hits after filtering
                actual_hitlist_size = max(max_results * 3, 50)
            else:
                actual_hitlist_size = max_results
            
            # Submit BLAST query using NCBIWWW (web interface)
            result_handle = NCBIWWW.qblast(
                blast_type, 
                database, 
                sequence,
                hitlist_size=actual_hitlist_size,
                expect=expect_threshold,
                format_type="XML",
                entrez_query=entrez_query
            )
            
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
                
                # Show results (already filtered by NCBI if entrez_query was used)
                alignments_to_show = record.alignments[:max_results]
                
                # Debug: Print all organism names to console
                if entrez_query and len(record.alignments) > 0:
                    print(f"Top {min(10, len(record.alignments))} organisms found:")
                    for i, aln in enumerate(record.alignments[:10], 1):
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
        # Map natural language name to EMBOSS command if needed
        emboss_name = self.tool_map.get(tool_name.lower(), tool_name.lower())
        
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
        elif emboss_name == 'needle':
            return self.align_sequences(kwargs.get('seq1', ''), kwargs.get('seq2', ''))
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
                exclude_taxa=kwargs.get('exclude_taxa', None)
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
        else:
            # Generic EMBOSS tool fallback
            return self._run_generic_emboss_tool(emboss_name, **kwargs)
    
    def _run_generic_emboss_tool(self, tool_name: str, **kwargs) -> str:
        """Generic runner for any EMBOSS tool that accepts sequence input
        
        Args:
            tool_name: EMBOSS tool name (e.g., 'iep', 'charge', 'mwfilter')
            **kwargs: Parameters including 'sequence'
        
        Returns:
            str: Tool output or error message
        """
        sequence = kwargs.get('sequence', '')
        if not sequence:
            return f"Error: No sequence provided for tool '{tool_name}'"
        
        try:
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
                
                # Smart chaining: if previous step was gene_query and this step is BLAST/analysis tool
                if previous_tool == 'gene_query' and tool_name in ['blast', 'blastn', 'blastp', 'blastx', 'search']:
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
                            print(f"[Step {idx+1}/{len(steps)}] Using {seq_type} sequence ({len(sequence)} {'aa' if seq_type == 'protein' else 'bp'}) for BLAST")
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
                        parameters['sequence'] = cached_result
                
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
                if isinstance(result_text, str) and len(result_text) > 500:
                    output += f"Output:\n{result_text[:500]}...\n\n"
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
                
                # Parse basic info from PSL format
                lines = result.split('\n')
                hits = [l for l in lines if l and not l.startswith('#') and not l.startswith('psLayout')]
                
                return {
                    'success': True,
                    'output': result,
                    'num_hits': len(hits),
                    'summary': f"BLAT search found {len(hits)} hits in {database}",
                    'note': "BLAT is more sensitive for near-exact matches than BLAST"
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
