#!/usr/bin/env python

import sys
import math
from Bio import SeqIO

N = 0
L = 0

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    L = L + len(seq_record)
    N = N + 1

print(min(18, math.floor(math.log2(max(L / N, 100)))))
