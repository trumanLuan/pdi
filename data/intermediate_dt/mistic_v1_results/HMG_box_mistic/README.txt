MISTIC ReadMe File:

You have downloaded results obtained using MISTIC Server.
Below is a description for each file (depending on the type of data submitted for the job, some files may not be available).
For more information visit http://mistic.leloir.org.ar/help.php


- ***.fasta:
This is your original Multiple Sequence Alignment (MSA) as it was received by the server. 
If you uploaded your MSA, it will have the same name as the file. 
If it was uploaded via using copy&paste on the text box, it will be named as "msa_tmpdir***.fasta". 
If you used a Pfam Id, it will be named as the Id. For example, for family PF00001, the MSA will be PF00001.fasta


- msa_gapstrip.fasta:
If you selected a reference sequence, the resulting MSA will be stored in this fasta file.
Sequences that are composed of more than 50% of gaps will be removed (this is, 50% of the length of the reference sequence)


- MI_data:
Raw Mutual Information file
Comment lines start with a '#'
Information about parameters used, time and date will be at the header
Sequences removed during the trimming of the reference sequence will be annotated.(see 'Reference Sequence')
Example:
# # Remove ***fasta id*** 17 16 

The first number is the number of gaps found for that sequence and the second one is the number of gaps allowed.
Then follows clustering percent identity used (Cluster ID) and number of clusters found (Ncluster).

Z-score Mutual Information values per pair of positions in the MSA are found after the comment lines.
resnum1 (column 1): Position 1 in the Mutual Information (MI) pair (numbered from 1 to N)
resname1 (column 2)*: Aminoacid of the reference sequence (one-letter-code) in position 1 of the MI pair
resnum2 (column 3): Position 2 in the MI pair
resname1 (column 4)*: Aminoacid of the reference sequence (one-letter-code) in position 2 of the MI pair
MIvalue (column 5): Z-score MI value for that pair of positions. 
If value -999.990000 is found, it means that no MI was calculated for that pair of positions.
This may be due to too many gaps for one of the positions. (See 'Gap Removal')


- pdb_MI_data:
Raw Mutual Information file aligned to the PDB sequence (if a pdb was submitted)
MI seq: sequence for the 'reference sequence'.
PDB (PDBID) (CHAIN) nres
PDBseq: sequence of the PDB
QAL MI_entry: region aligned of the 'reference sequence'
DAL (PDBID):  region aligned of the PDB sequence
pdbnum1 (column 1): Position 1 in the Mutual Information (MI) pair, using the numbering of the PDB.
resname1 (column 2): Aminoacid of the reference sequence (one-letter-code) in position 1 of the MI pair
pdbnum2 (column 3): Position 1 in the Mutual Information (MI) pair, using the numbering of the PDB
resname2 (column 4): Position 2 in the MI pair, using PDB numbering
MIvalue (column 5): Z-score MI value for that pair of positions. 
Distance (column 6): Minimum distance between heavy atoms of the two residues in the PDB file

If '-9' is found in pdbnum1 or pdbnum2 it means there is a gap in the alignment between the 'reference sequence' and the PDB sequence. Also, MI values and distance
values for the pair are set to -999.990000 and -99.999 respectively.


- conservation_data:
Conservation per MSA column using Kullback-Leibler.
column 1: Column of the MSA
column 2*: Aminoacid in the reference sequence (one-letter-code)
column 3: Kullback-Leibler divergence value
column 5-24: Aminoacid relative frequencies for that column of the MSA



-pdb_conservation_data:
Same as before but using PDB numbering
Positions with '-9' mean that there is a gap in the alignment between 'reference sequence' and PDB sequence


- CMI_data:
Per position accumulation of MI (PDB numbering is used when possible)

- pMI_data:
Per position accumulation of Cumulative MI.
resnum (column 1): Column of the MSA
resname (column 2): Aminoacid of the reference sequence (one-letter-code)
pdbnum (column 3): Corresponding column of the MSA using PDB numbering
CMI (column 4): Cumulative MI value
Interactions (column 5): Number of MI interactions with values >6.5
pMI (column 6): Proximity MI value. 
Contacts (column 5): Number of contacts in the structure at <5 Angstroms



(*) If no reference sequence was selected, characters for the first sequence in the alignment will be shown, including gaps.


How to cite us:
"Mistic: Mutual Information Server to Infer Coevolution". FL Simonetti, E Teppa, A Chernomoretz, M Nielsen, C Marino Buslje


