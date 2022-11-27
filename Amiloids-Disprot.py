from pyfaidx import Fasta
import csv

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

tsv_list = {}
count = 0
with open('DisProt release_2022_06.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    for rows in file:
        while '' in rows:
            rows.remove('')
        for i in range(0, len(rows)):
            tsv_list[rows[len(rows)-1]] = rows[0]

new_tsv = {}
short_amiloid_list_for_head = {}
for key, seq in seq_1000.items():
    short_amiloid_list = []
    for val in waltz_dict.values():
        if val in str(seq):
            #if str(seq) not in filter_seq.values():
            filter_seq[key] = seq
            short_amiloid_list.append(val)
            short_amiloid_list_for_head[key] = short_amiloid_list
            for pep, uni in tsv_list.items():
                if str(seq) in pep:
                    new_tsv[pep] = uni

delet_repeats = []
for key, seq in filter_seq.items():
    delet_repeats.append([key, str(seq)])

delet_repeats2 = []
for couples in delet_repeats:
    count = 0
    count2 = 0
    for coup in delet_repeats2:
        if couples[1] in coup[1]:
            count += 1
        elif coup[1] in couples[1]:
            count2 += 1
    if count == 0:
        delet_repeats2.append(couples)


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
fasta_uniprot = {}
for key, val in amiloid_seq_no_amil_disc.items():
    head_for_fasta_new = key
    head_for_fasta = key
    amils_for_uniprot = ''
    uni_for_uniprot = ''
    pep_for_uniprot = ''
    for pep, uni in new_tsv.items():
        if str(val) in pep:
            head_for_fasta_new = '>UniProt|' + uni + ' ' + head_for_fasta.replace('>', '')
            uni_for_uniprot = '>UniProt|' + uni
            pep_for_uniprot = pep
    for key_amil, val_amil in short_amiloid_list_for_head.items():
        if key_amil in key:
            for i in val_amil:
                if ' amilod_seq=' not in head_for_fasta_new:
                    head_for_fasta_new += ' amilod_seq=' + i
                else:
                    head_for_fasta_new += ', ' + i
                if ' amilod_seq=' not in uni_for_uniprot:
                    uni_for_uniprot += ' amilod_seq=' + i
                else:
                    uni_for_uniprot += ', ' + i
    fasta_uniprot[pep_for_uniprot] = uni_for_uniprot
    amil_seq[head_for_fasta_new] = val

'''
fasta_uniprot_no_repeats = {}
for pep, uni in fasta_uniprot.items():
    fasta_uniprot_no_repeats[uni] = pep
for uni, pep in fasta_uniprot_no_repeats.items():
    print(uni, pep, sep='\n')

amil_seq_no_repeats = {}
for key, val in amil_seq.items():
    amil_seq_no_repeats[str(val)] = key
'''

table = []
with open('DisProt release_2022_06.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    count = 0
    for rows in file:
        if count == 0:
            table.append(rows)
            count += 1
        for pep, uni in new_tsv.items():
            if uni in rows and pep in rows:
                table.append(rows)

with open('table_amiloids.tsv', 'w') as tab:
    for elem in table:
        for words in elem:
            tab.write(words + ', ')
        tab.write('\n')

with open('Amiloid+Disprot_full_corr', 'w') as A_D:
    for key, val in amil_seq.items():
        A_D.write(str(val) + '\n')
        A_D.write(key + '\n')