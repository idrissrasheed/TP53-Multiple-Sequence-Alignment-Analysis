import os
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq


input_file = "TP53_Gene.fasta"
records = SeqIO.parse(input_file, "fasta")
records = list(records) # make a copy

maxlen = max(len(record.seq) for record in records)

# Pad sequences to have same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '.')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# Write to temporary file and print alignment
output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')
alignment = AlignIO.read(output_file, "fasta")
print(alignment)


