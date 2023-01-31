# Algorithm for searching of short amyloid motifs from Waltz-DB in DisProt database

---

Material for master's degree work in ITMO University.

Student: Zhiltsova T.

Supervisor: Kajava A., Pyankov I.

## Abstract

---

The main purpose of this work in general is to expand amount of known amyloid forming peptides. 
This may help to improve  prediction quality of amyloid predictors. 
On this step our tusk was to develop algorithm of searching of short amyloid motifs obtained 
from Waltz-DB in peptide sequences from DisProt. In DisProt presented sequences which are partly or fully unstructured, 
so this let us suppose that motifs in this regions will aggregate with higher chances, then same motifs in structured 
parts of peptide. Also, we're studied only peptides which are shorter than 1000 residues to prevent impact of other 
part of peptide on amyloid region.   

## Results

---

We found 14 peptides in DisProt database that contain amyloid motifs from Waltz-DB. Half of them have already been studied
on amyloid activity while half of them has a potential for further experiments on their biological properties. Our results
are presented in "results.zip". 

## User Guid

---

Download [Waltz-DB amyloid hexapeptides dataset](http://waltzdb.switchlab.org/sequences) (switch "Classification" 
to "amyloid") in csv format and [DisProt dataset](https://disprot.org/download) in fasta and tsv format
and place it in one directory with all supplement files from this repo. 

### Launch
All algorithm is presented in Amiloids-Disprot.py. You can launch it after making the previous step. No any special 
preparation are needed.

Output fasta file will be named "Amyloid+Disprot_full_corr" and will contain all sequences with amyloid motifs (they will
be shown in descriptor). Sheet with list and UniProt/DisProt code of found peptides will be named 
"UniProt+DisProt_amyloids". 

All requirements are shown in "requirements.txt".

## References

---

1. Jacinte Beerten, J. WALTZ-DB: a benchmark database of amyloidogenic hexapeptides / J. Beerten, J. Van Durme, R. Gallardo, E. Capriotti, L. Serpell, F. Rousseau, J. Schymkowitz // Bioinformatics. – 2015. – Vol. 31. – p. 1698 – 1700.
2. Louros, N. WALTZ-DB 2.0: an updated database containing structural information of experimentally determined amyloid-forming peptides / N. Louros, K. Konstantoulea, M. De Vleeschouwer, M. Ramakers, J. Schymkowitz, F. Rousseau // Nucleic Acids Research. – 2020. – Vol. 48. – p. D389 – D393.
3. Sickmeier, M. DisProt: the Database of Disordered Proteins / M. Sickmeier, J.A. Hamilton, T. LeGall, V. Vacic, M.S. Cortese, A. Tantos, B. Szabo, P. Tompa, J. Chen, V.N. Uversky, Z. Obradovic, A.K. Dunker // Nucleic Acids Research. – 2007. – Vol. 35. – p. D786 – D793.
4. Hatos, A. DisProt: intrinsic protein disorder annotation in 2020 / A. Hatos, B. Hajdu-Soltész, A.M. Monzon, N. Palopoli, L. Álvarez, B. Aykac-Fas, C. Bassot, G.I. Benítez, M. Bevilacqua // Nucleic Acids Research. – 2020. – Vol. 48. – p. D269 – D276.
5. Shirley, M.D. Efficient "pythonic" access to FASTA files using pyfaidx / M.D. Shirley​, Z. Ma, B.S. Pedersen, S.J. Wheelan // PeerJ Preprints. – 2015.
6. H.M. Berman, H.M. The Protein Data Bank / H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne // Nucleic Acids Research. – 2000. – Vol. 28. – p. 235 – 242.
