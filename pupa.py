import csv

row = []
with open('DisProt release_2022_06.tsv') as tsv:
    file = csv.reader(tsv, delimiter='\t')
    for rows in file:
        while '' in rows:
            rows.remove('')
        row.append(rows)
for line in row:
    print(len(line))
