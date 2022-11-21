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

short_amiloid_list_for_head = {}
for key, seq in seq_1000.items():
    short_amiloid_list = []
    for val in waltz_dict.values():
        if val in str(seq):
            #if str(seq) not in filter_seq.values():
            filter_seq[key] = seq
            short_amiloid_list.append(val)
            short_amiloid_list_for_head[key] = short_amiloid_list

discriptor = []
with Fasta('DisProt_2022_06.fasta') as genes:
    for record in genes:
        line = record.long_name
        discriptor.append(line[:])

amiloid_seq_no_amil_disc = {}
for key, val in filter_seq.items():
    for disc in discriptor:
        if key in disc:
            amiloid_seq_no_amil_disc[disc] = val

amil_seq = {}
for key, val in amiloid_seq_no_amil_disc.items():
    head_for_fasta = key
    for key_amil, val_amil in short_amiloid_list_for_head.items():
        if key_amil in key:
            for i in val_amil:
                if ' amilod_seq=' not in head_for_fasta:
                    head_for_fasta += ' amilod_seq=' + i
                else:
                    head_for_fasta += ', ' + i
    amil_seq[head_for_fasta] = val

amil_seq_no_repeats = {}
for key, val in amil_seq.items():
    amil_seq_no_repeats[str(val)] = key


with open('Amiloid+Disprot_full_corr', 'w') as A_D:
    for key, val in amil_seq_no_repeats.items():
        A_D.write(str(val) + '\n')
        A_D.write(key + '\n')
