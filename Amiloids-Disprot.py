from pyfaidx import Fasta
import csv

# reading of the file with amiloids motifs
waltz_seq = []       # list for short motifs of amiloid
waltz_dict = {}      # dictionry for sequences from WaltzDB
with open('WALTZ_DB_amiloid_seq') as waltz_db:
    for line in waltz_db:
        line = line.strip()
        waltz_seq.append(line)
    for number, amiloid_motif in enumerate(waltz_seq):
        if amiloid_motif.strip().isalpha() == False:    # deleting empty lines
            break
        else:
            waltz_dict[number] = amiloid_motif

# filter all sequences which are shorter than 1000 residues
seq_1000 = {}  # dict for all sequences which are shorter than 1000 residues
with Fasta('DisProt_2022_06.fasta') as genes:
    for number_of_peptide in range((len(genes.keys()))):
        if genes[number_of_peptide][:].end < 1000:
            seq_1000[genes[number_of_peptide][:].name] = genes[number_of_peptide][:]

# reading tsv file, which contains UniProt names for all peptides from DisProt. They will be used in final fasta
tsv_list = {}  # dict for pairs "UniProt_name - Sequence"
with open('DisProt release_2022_06.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    for rows in file:
        while '' in rows:     # deleting empty elements from tsv
            rows.remove('')
        tsv_list[rows[len(rows)-1]] = rows[0]

# finding peptides with amiloid motifs
new_tsv = {}   # dict for remembering UniProt names of filtered peptides
amiloid_motifs_list_for_header = {}  # remembering amiloid motifs for each peptide, which have amiloid-forming region
filter_seq = {}  # dict for seq that have amiloid-forming regions
for key, seq in seq_1000.items():
    amiloid_motifs_list = []   # updated list for all amiloid motifs in one peptide
    for val in waltz_dict.values():
        if val in str(seq):  #
            filter_seq[key] = seq
            amiloid_motifs_list.append(val)
            for pep, uni in tsv_list.items():
                if str(seq) in pep:
                    new_tsv[pep] = uni
        amiloid_motifs_list_for_header[key] = amiloid_motifs_list

'''
#   algorithm for deleting repeats in fasta
filter_seq_no_repeats = {}
for key, seq in filter_seq.items():
    r = 0
    new_key = str(key).replace('disprot|', '')
    finish = str(new_key).find('r')
    for dp_id, am_a in filter_seq_no_repeats.items():
        if new_key[:finish] == dp_id:
            r += 1
            if str(am_a) == str(seq) or len(str(am_a)) > len(str(seq)):
                continue
            elif len(str(am_a)) < len(str(seq)):
                filter_seq_no_repeats['disprot|' + new_key[:finish]] = seq
    if r == 0:
        filter_seq_no_repeats['disprot|' + new_key[:finish]] = seq

'''
discriptor = []
with Fasta('DisProt_2022_06.fasta') as genes:
    for record in genes:
        line = record.long_name
        '''
        finding_pos = str(line)[str(line).find('pos'):]
        finding_pos = finding_pos.replace(finding_pos[:finding_pos.find(' ')], '')
        finding_pos_2 = str(line)[:str(line).find('pos')]
        final_disc = finding_pos_2 + finding_pos
        '''
        discriptor.append(str(line))

amiloid_seq_no_amil_disc = {}
for key, val in filter_seq.items():
    for disc in discriptor:
        if key in disc:
            amiloid_seq_no_amil_disc[disc] = val

amil_seq = {}
fasta_uniprot = {}
for key, val in amiloid_seq_no_amil_disc.items():
    head_for_fasta = str(key)
    for pep, uni in new_tsv.items():
        print(uni)
        if str(val) in pep:
            head_for_fasta_new = '>UniProt|' + uni + ' ' + head_for_fasta.replace('>', '')
    for key_amil, val_amil in amiloid_motifs_list_for_header.items():
        if key_amil in key:
            for i in val_amil:
                if ' amilod_seq=' not in head_for_fasta_new:
                    head_for_fasta_new += ' amilod_seq=' + i
                else:
                    head_for_fasta_new += '|' + i
    amil_seq[head_for_fasta_new] = val

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

table_short = []
t = 1
for line in table:
    uniprot = line.index('acc')
    name = line.index('name')
    organism = line.index('organism')
    disprot_id = line.index('disprot_id')
    break

for elem in table:
    count_short_tab = 0
    table_short_words = []
    for lines in table_short:
        if elem[uniprot] in lines:
            count_short_tab += 1
    if count_short_tab == 0:
        table_short_words.append(elem[uniprot])
        table_short_words.append(elem[name])
        table_short_words.append(elem[organism])
        table_short_words.append(elem[disprot_id])
        table_short.append(table_short_words)

with open('table_amiloids.tsv', 'w') as tab:
    for elem in table:
        for words in elem:
            tab.write(words + ', ')
        tab.write('\n')

with open('Amiloid+Disprot_full_corr', 'w') as A_D:
    for key, val in amil_seq.items():
        A_D.write(str(key) + '\n')
        A_D.write(str(val) + '\n')

with open('UniProt+DisProt_amiloids.tsv', 'w') as tab:
    for elem in table_short:
        for words in elem:
            tab.write(words + ', ')
        tab.write('\n')
