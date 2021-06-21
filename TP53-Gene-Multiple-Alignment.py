import os
from Bio import AlignIO
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline

# Provide an email address
email = input("Enter email address ")

Entrez.email = email

# Enter data from Genbank
with Entrez.einfo(db="nucleotide") as handle:
    record = Entrez.read(handle)
for field in record["DbInfo"]["FieldList"]:
    print("%(Name)s, %(FullName)s, %(Description)s" % field)

# Gather all accession numbers for TP53 homo sapiens 
with Entrez.esearch(db="nucleotide", 
                    term="TP53 – tumor protein p53 Homo sapiens (human)", 
                    idtype="acc", 
                    retmax="10000") as handle:
    results = Entrez.read(handle)
accs = results["IdList"] #save the list of accession numbers

#Use efetch to download records and save to directory
filename = "TP53_Gene.gbk"
if not os.path.isfile(filename):
    # Downloading...
    with Entrez.efetch(db="nucleotide", 
                       id=accs, 
                       rettype="gb", 
                       retmode="text") as net_handle: 
        with open(filename, "w") as out_handle:
            out_handle.write(net_handle.read())
    print("Saved")
    
#Read in unfiltered data
unfiltered = SeqIO.parse("TP53_Gene.gbk", "genbank")
#Drop data without (close to) full length sequences
full_length_records = []
for record in unfiltered:
    if len(record.seq) > 29000:
        full_length_records.append(record)
#Write filtered data to file
SeqIO.write(full_length_records, "TP53_Gene.fasta", "fasta")
#Align sequences with MUSCLE (using parameters to make the alignment
#process as fast as possible)
muscle_cline = MuscleCommandline(input="TP53_Gene.fasta", 
                                 out="TP53_Gene_aligned.fasta", 
                                 diags = True, 
                                 maxiters = 1, 
                                 log="align_log.txt")
#Print out filenames 
print(muscle_cline)

align = AlignIO.read("TP53_Gene.fasta")