from pyfaidx import Fasta


waltz_seq = ()
waltz_dict = {}    #dictionry for sequences from WaltzDB
with open('WALTZ_DB_amiloid_seq') as waltz_db:
    waltz_seq = waltz_db.readlines()  # reading seq from the file
    for i, count in enumerate(waltz_seq):
        if count.strip().isalpha() == False:  # deleting empty lines
            break
        else:
            waltz_dict[i] = count
seq_1000 = dict()  # dict for all sequences which are shorter than 1000 residues
filter_seq = dict()  # dict for seq that have amiloid-forming region
with Fasta('DisProt_2022_06.fasta') as genes:
    for i in range((len(genes.keys()))):
        if genes[i][:].end < 1000:
            seq_1000[genes[i][:].name] = genes[i][:]
for key, seq in seq_1000.items():
    for val in waltz_dict.values():
        if val in seq:
            filter_seq[key] = seq
for key, val in filter_seq.items():
    print(key, val)
