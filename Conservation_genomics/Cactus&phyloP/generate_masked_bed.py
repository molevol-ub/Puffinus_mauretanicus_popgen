#!/usr/bin/env python

import sys

chrom = ""
pos = -1
start = -1
in_masked_region = False
infile = sys.argv[1]
outfile = sys.argv[2]

lowcase=["a","c","t","n","g"]

outwirte=open(outfile, "w")

with open(sys.argv[1], "r") as fh:
    for line in fh:
        if line.startswith(">"):
            if in_masked_region:  # last masked region from previous chrom
                outwrite.write(f"{chrom}\t{start}\t{pos}")
                start = -1  # not needed actually
                in_masked_region = False
            pos = 0
            chrom = line.split(" ")[0].replace(">", "").strip()
        else:
            for c in line.strip():
                if not in_masked_region and c in lowcase:
                    in_masked_region = True
                    start = pos
                elif in_masked_region and c not in lowcase:
                    in_masked_region = False
                    outwrite.write(f"{chrom}\t{start}\t{pos}")

                pos += 1

if in_masked_region:  # last masked region in last chrom
    outwrite.write(f"{chrom}\t{start}\t{pos}")
