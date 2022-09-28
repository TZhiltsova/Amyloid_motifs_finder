from pyfaidx import Fasta

genes = Fasta('/home/lora/Amilods/Disprot/Amiloids-Disprot/DisProt_2022_06.fasta')
waltz_db = open('/home/lora/Amilods/Disprot/Amiloids-Disprot/WALTZ_DB_amiloid_seq')
waltz_seq = []    # list for seq from waltz db
seq_1000 = dict()  # dict for all sequences which are shorter than 1000 residues
filter_seq = dict()   # dict for seq that have amiloid-forming region
i = 0 # counter for number of peptide in fasta file
for elem in waltz_db:
    waltz_seq.append(elem)     # putting elements of file into list
while i < len(genes.keys()):    # deleting all sequences longer then 1000 amino acids
    if genes[i][:].end < 1000:
        seq_1000[genes[i][:].name] = genes[i][:]
    i += 1
for key in seq_1000:
    for seq in waltz_seq:
        for a in range(len(seq_1000[key])-1):       # frame for finding of amiloid-forming sequences in disprot peptides
            if seq_1000[key][a:len(seq)] == seq:
                filter_seq[key] = seq_1000[key]
for key in filter_seq:
    print(key, filter_seq[key])