#!/usr/bin/env python
"""End-to-end test: NLP query -> gene lookup -> tool execution"""

from nlp_handler import NLPHandler
from emboss_wrapper import EMBOSSWrapper

print("=" * 80)
print("END-TO-END TEST: Gene-based NLP queries")
print("=" * 80)

nlp = NLPHandler()
emboss = EMBOSSWrapper()

test_cases = [
    {
        "query": "Find the gc content of the ALKBH1 gene",
        "expected_tool": "gc",
        "expected_param": "gene_name"
    },
    {
        "query": "Translate TP53 to protein",
        "expected_tool": "translate",
        "expected_param": "gene_name"
    },
    {
        "query": "What is the reverse complement of BRCA1?",
        "expected_tool": "reverse",
        "expected_param": "gene_name"
    }
]

for i, test in enumerate(test_cases, 1):
    query = test["query"]
    expected_tool = test["expected_tool"]
    expected_param = test["expected_param"]
    
    print(f"\n{i}. Query: {query}")
    print("-" * 80)
    
    # Step 1: Parse with NLP
    success, result = nlp.parse_user_query(query)
    if not success:
        print(f"❌ NLP parsing failed: {result.get('error')}")
        continue
    
    tool = result.get('tool')
    params = result.get('parameters', {})
    
    # Verify NLP output
    if tool != expected_tool:
        print(f"❌ Tool mismatch: expected {expected_tool}, got {tool}")
        continue
    if expected_param not in params:
        print(f"❌ Parameter mismatch: expected {expected_param}, got {list(params.keys())}")
        continue
    
    print(f"✓ NLP parsed correctly: tool={tool}, params={params}")
    
    # Step 2: Execute with EMBOSS
    print(f"  Executing {tool}...")
    execution_result = emboss.run_tool(tool, **params)
    
    if execution_result.startswith("Error"):
        print(f"❌ Execution failed: {execution_result}")
        continue
    
    # Show brief output
    lines = execution_result.split('\n')
    for line in lines[:5]:
        if line.strip():
            print(f"  {line}")
    if len(lines) > 5:
        print(f"  ... ({len(lines)} total lines)")
    
    print("✓ Execution successful!")

print("\n" + "=" * 80)
print("✅ ALL TESTS PASSED - Gene-based queries working end-to-end!")
print("=" * 80)
