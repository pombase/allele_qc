"""
Split the peptide and gene sequences including introns and UTRs from peptide.fa and
cds+introns+utrs.fa into individual files named `systematic gene name.fa`, and contain
either the DNA or peptide sequence in fasta format.

E.g. for the first gene:

File data/peptide_sequences/SPAC1002.01.fa contains:

```
>SPAC1002.01.1:pep mrx11|mitochondrial expression network (MIOREX) component Mrx11
MLPPTIRISGLAKTLHIPSRSPLQALKGSFILLNKRKFHYSPFILQEKVQSSNHTIRSDT
KLWKRLLKITGKQAHQFKDKPFSHIFAFLFLHELSAILPLPIFFFIFHSLDWTPTGLPGE
YLQKGSHVAASIFAKLGYNLPLEKVSKTLLDGAAAYAVVKVC*
```

File data/nucleotide_sequences/SPAC1002.01.fa


"""

import re

with open('peptide.fa') as ins: