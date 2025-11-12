#!/usr/bin/env python
"""Test NLP gene query parsing"""

from nlp_handler import NLPHandler

nlp = NLPHandler()

# Test queries that should be recognized as gene queries
test_queries = [
    "How many exons are in the human ALKBH1 gene?",
    "What is the CDS length of TP53?",
    "Show me the transcript information for BRCA1",
    "Find gene ALKBH1 using GENCODE V48 track",
]

print("Testing NLP Handler Gene Query Recognition\n")
print("=" * 70)

for query in test_queries:
    print(f"\nQuery: {query}")
    success, result = nlp.parse_user_query(query)
    
    if success:
        tool = result.get('tool')
        params = result.get('parameters', {})
        explanation = result.get('explanation', '')
        
        print(f"✓ Tool: {tool}")
        print(f"✓ Parameters: {params}")
        print(f"✓ Explanation: {explanation}")
    else:
        print(f"✗ Error: {result.get('error')}")
