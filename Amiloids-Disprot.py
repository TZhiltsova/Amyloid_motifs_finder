from pyfaidx import Fasta
import csv

def waltz_reader(waltz_pass: str) -> dict:
    '''
        :parameter waltz_pass: pass to the Waltz_DB file
        :return: Dictionary with list of all amyloid motifs from Waltz_DB file
    '''
    waltz_seq = []
    waltz_dict = {}
    with open(waltz_pass) as waltz_db:
        waltz_seq = [line.strip() for line in waltz_db]
        for number, amyloid_motif in enumerate(waltz_seq):
            if amyloid_motif.strip().isalpha() == False:    # deleting empty lines
                break
            else:
                waltz_dict[number] = amyloid_motif
    return waltz_dict
def fasta_reader(fasta_pass: str) -> dict:
    '''
    :param fasta_pass: pass to the studying fasta file
    :return: two parameters: dictionary with all sequences from fasta as values and DisProt code as keys
                            dictionary with full description of peptides as key instead of short code
    '''
    seq_1000 = {}
    discriptor = []
    with Fasta(fasta_pass) as genes:
        for record in genes:
            line = record.long_name
            discriptor.append(str(line))
        for number_of_peptide in range((len(genes.keys()))):
            if genes[number_of_peptide][:].end < 1000:
                seq_1000[genes[number_of_peptide][:].name] = genes[number_of_peptide][:]
    return seq_1000, discriptor
def tsv_reader(tsv_pass: str) -> dict:
    '''
    :param tsv_pass: takes tsv file where sequences presented in according to there UniProt code
    :return: dictionary with all sequences from fasta as values and UniProt code as keys
    '''
    tsv_list = {}
    with open(tsv_pass) as tsv:
        file = csv.reader(tsv, delimiter='\t')
        for list_of_elem in file:
            while '' in list_of_elem:     # deleting empty elements from tsv
                list_of_elem.remove('')
            tsv_list[list_of_elem[len(list_of_elem)-1]] = list_of_elem[0]
    return tsv_list
def am_finder(waltz_dict, seq_1000, tsv_list) -> dict:
    '''
    :param waltz_dict: Dictionary with list of all amyloid motifs from Waltz_DB file
    :param seq_1000: dictionary with all sequences from fasta as values and DisProt code as keys
    :param tsv_list: dictionary with all sequences from fasta as values and UniProt code as keys
    :return: three parameters: dictionary with peptide in which sequences amyloid motifs where found where DisProt code is a key
                            dictionary with peptide in which sequences amyloid motifs where found where UniProt code is a key
                            dictionary with amyloid motifs in according to peptide in which sequences it was found
    '''
    new_tsv = {}
    amiloid_motifs_list_for_header = {}
    filter_seq = {}
    for key, sequence in seq_1000.items():
        amiloid_motifs_list = []
        for val in waltz_dict.values():
            if val in str(sequence):  #
                filter_seq[key] = sequence
                pos = str(sequence).find(val) + 1
                val = str(pos) + val
                amiloid_motifs_list.append(val)
                for peptide_seq, uniprot_number in tsv_list.items():
                    if str(sequence) in peptide_seq:
                        new_tsv[peptide_seq] = uniprot_number
            amiloid_motifs_list_for_header[key] = amiloid_motifs_list
    return filter_seq, new_tsv, amiloid_motifs_list_for_header

def creating_fasta(filter_seq, discriptor, new_tsv):
    '''
    :param filter_seq: dictionary with peptide in which sequences amyloid motifs where found where DisProt code is a key
    :param discriotor: dictionary with full description of peptides as key instead of short code
    :param new_tsv: dictionary with peptide in which sequences amyloid motifs where found where UniProt code is a key
    :return: fasta file with all peptides in which sequences amyloid motifs where found (with repeats)
    '''
    amiloid_seq_no_amil_disc = {}
    for key, val in filter_seq.items():
        for disc in discriptor:
            if key in disc:
                amiloid_seq_no_amil_disc[disc] = val

    amil_seq = {}
    fasta_uniprot = {}
    for key, val in amiloid_seq_no_amil_disc.items():
        head_for_fasta = str(key)
        for peptide_seq, uniprot_number in new_tsv.items():
            if str(val) in peptide_seq:
                head_for_fasta_new = head_for_fasta + ' ' + 'UniProt|' + uniprot_number
        for key_amil, val_amil in amiloid_motifs_list_for_header.items():
            if key_amil in key:
                for i in val_amil:
                    if ' amyloid_seq=' not in head_for_fasta_new:
                        head_for_fasta_new += ' amyloid_seq=' + i
                    else:
                        head_for_fasta_new += '|' + i
        amil_seq[head_for_fasta_new] = val

    with open('Amyloid+Disprot_full_corr.fasta', 'w') as amyl_disprot:
        for key, val in amil_seq.items():
            amyl_disprot.write(str(key) + '\n')
            amyl_disprot.write(str(val) + '\n')
    return

def creating_tsv(tsv_pass, new_tsv):
    '''
    :param tsv_pass: takes tsv file where sequences presented in according to there UniProt code
    :param new_tsv: dictionary with peptide in which sequences amyloid motifs where found where UniProt code is a key
    :return: tsv file with peptides in which sequences amyloid motifs where found (without repeats)
    '''
    table = []
    with open(tsv_pass) as tsv:
        file = csv.reader(tsv, delimiter='\t')
        count = 0
        for rows in file:
            if count == 0:
                table.append(rows)
                count += 1
            for peptide_seq, uniprot_number in new_tsv.items():
                if uniprot_number in rows and peptide_seq in rows:
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

    with open('UniProt+DisProt_amyloids.tsv', 'w') as tab:
        for elem in table_short:
            for words in elem:
                tab.write(words + '\t')
            tab.write('\n')
    return

# Waltz_pass = input(str('print'+'\n'))
waltz_pass = 'WALTZ_DB_amiloid_seq'
waltz_dict = waltz_reader(waltz_pass)
fasta_pass = 'DisProt_2022_06.fasta'
seq_1000, discriptor = fasta_reader(fasta_pass)
tsv_pass = 'DisProt release_2022_06.tsv'
tsv_list = tsv_reader(tsv_pass)
filter_seq, new_tsv, amiloid_motifs_list_for_header = am_finder(waltz_dict, seq_1000, tsv_list)
creating_fasta(filter_seq, discriptor, new_tsv)
creating_tsv(tsv_pass, new_tsv)