# Algorithm for searching of short amiloid motifs from Waltz-DB in DisProt database

Material for master's degree work in ITMO University.

Student: Zhiltsova T.

Supervisor: Kajava A., Pyankov I.

## Abstract

The main purpose of this work in general is to expand amount of known amyloid forming peptides. 
This may help to improve  prediction quality of amyloid predictors. 
On this step our tusk was to develop algorithm of searching of short amyloid motifs obtained 
from Waltz-DB in peptide sequences from DisProt. In DisProt presented sequences which are partly or fully unstructured, 
so this let us suppose that motifs in this regions will aggregate with higher chances, then same motifs in structured 
parts of peptide. Also, we're studied only peptides which are shorter than 1000 residues to prevent impact of other 
part of peptide on amyloid region.   

## User Guid

Download [Waltz-DB amyloid hexapeptides dataset](http://waltzdb.switchlab.org/sequences) (switch "Classification" 
to "amyloid") in csv format and [DisProt dataset](https://disprot.org/download) and place it in one directory with all 
supplement files from this repo. 

### Launch
