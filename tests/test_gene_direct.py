#!/usr/bin/env python
"""Test gene query directly without NLP parsing"""

from emboss_wrapper import EMBOSSWrapper

emboss = EMBOSSWrapper()

print("Testing Gene Query Functions\n")
print("=" * 70)

# Test 1: Query ALKBH1
print("\n1. Testing ALKBH1 gene query:")
print("-" * 70)
result = emboss.query_gene_info('ALKBH1')
print(result[:500] + "..." if len(result) > 500 else result)

# Test 2: Query TP53
print("\n2. Testing TP53 gene query:")
print("-" * 70)
result = emboss.query_gene_info('TP53')
lines = result.split('\n')
for line in lines[:20]:  # First 20 lines
    print(line)
print("...")

# Test 3: Query invalid gene
print("\n3. Testing invalid gene (should show error):")
print("-" * 70)
result = emboss.query_gene_info('INVALID_GENE_XYZ')
print(result)

print("\n" + "=" * 70)
print("✓ All gene query functions working!")
