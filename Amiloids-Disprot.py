from pyfaidx import Fasta
genes = Fasta('/home/lora/Amilods/Disprot/Amiloids-Disprot/DisProt_2022_06.fasta')
seq_1000 = []
for key in genes:
    if len(key) < 1000:
        seq_1000.append(key)
for key in seq_1000:
    print(key)