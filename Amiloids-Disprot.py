from pyfaidx import Fasta

genes = Fasta('/home/lora/Amilods/Disprot/Amiloids-Disprot/DisProt_2022_06.fasta')
seq_1000 = dict()  # list for all sequences which are shorter than 1000 residues
i = 0 # counter for number of peptide in fasta file
while i < len(genes.keys()):
    if genes[i][:].end < 1000:                       # deleting all sequences longer then 1000 amino acids
        seq_1000[genes[i][:].name] = genes[i][:]
    i += 1
for key in seq_1000:
    print(key, seq_1000[key])
