from pyfaidx import Fasta


waltz_seq = []
waltz_dict = {}    #dictionry for sequences from WaltzDB
with open('WALTZ_DB_amiloid_seq') as waltz_db:
    for line in waltz_db:
        line = line.strip()
        waltz_seq.append(line) # reading seq from the file
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
        if val in str(seq):
            #if str(seq) not in filter_seq.values():
            filter_seq[key] = seq

discript = []
with Fasta('DisProt_2022_06.fasta') as genes:
    for record in genes:
        line = record.long_name
        discript.append(line[:])

amiloid_seq = {}
for key, val in filter_seq.items():
    for disc in discript:
        if key in disc:
            amiloid_seq[disc] = val


with open('Amiloid+Disprot_full_corr', 'w') as A_D:
    for key, val in amiloid_seq.items():
        A_D.write(key + '\n')
        A_D.write(str(val) + '\n')
        A_D.write('\n')
