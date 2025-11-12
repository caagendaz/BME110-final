#!/usr/bin/env python
"""Test genome query NLP parsing"""

from nlp_handler import NLPHandler

def test_genome_query():
    nlp = NLPHandler()
    
    # Test a genome query
    query = 'Get the sequence from human genome hg38, chromosome 1, from position 1000000 to 1001000'
    print(f'Query: {query}\n')
    success, result = nlp.parse_user_query(query)
    if success:
        print(f'✓ Tool: {result.get("tool")}')
        print(f'✓ Parameters: {result.get("parameters")}')
        print(f'✓ Explanation: {result.get("explanation")}')
    else:
        print(f'✗ Error: {result.get("error")}')

if __name__ == "__main__":
    test_genome_query()
