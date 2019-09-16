#!/usr/bin/env python3

from Bio import SeqIO

# orthogroup_file = 'data/MultipleSequenceAlignments/OG0000016.fa'
# proteins_file = 'test.fa'

orthogroup_file = snakemake.input[0]
proteins_file = snakemake.output[0]

wasp_proteins = []
wasp_prefixes = ['Vger', 'Vpen', 'Vvul']

orthogroup_seqs = list(SeqIO.parse(orthogroup_file, 'fasta'))
for og_seq in orthogroup_seqs:
    if og_seq.id[0:4] in wasp_prefixes:
        og_seq.seq = og_seq.seq.ungap('-')
        wasp_proteins.append(og_seq)

SeqIO.write(wasp_proteins, proteins_file, 'fasta')
