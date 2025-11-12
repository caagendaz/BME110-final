#!/usr/bin/env python
"""Test gene query"""

from emboss_wrapper import EMBOSSWrapper

emboss = EMBOSSWrapper()

print("Testing gene query for ALKBH1...\n")
result = emboss.query_gene_info('ALKBH1', genome='hg38', track='gencode')
print(result)
