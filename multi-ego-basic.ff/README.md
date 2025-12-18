This is the force field including the bonded multi-eGO parameters.

AMINOACIDS
----------
These are the guiding lines to parametrize the backbone dihedral angles:
0. steric interactions are indentified following  

1. Ideal geometries for secondary structure are use to guide the centers of the minima:

* alpha helix : -60, -45
* 3_10 helix: -49, -26 
* strand: –135, 135
* left-alpha: 45,60

2. Three classes of alpha/beta population:
Jiang, et al. Biopolymers 1998
Hayward, et al. J Struct Biol 2021

3. 3J coupling are matched using Karplus relation with parameters from:
#https://imserc.northwestern.edu/guide/eNMR/proteins/J/HNHA.html

J coupling experimental table from https://doi.org/10.1073/pnas.0510420103:
ALA  6.06 
ARG  6.85 
ASN 7.45  
ASP  6.93 
CYS 7.31  
GLN  7.14 
GLU  6.63 
GLY 5.85  
HIS  7.89 
ILE  7.33 
LEU  6.88 
LYS  6.83 
MET  7.02 
PHE  7.18 
SER  7.02 
THR  7.35 
TRP  6.91 
TYR  7.13 
VAL  7.3  

4. helical populations are fine tuned with Ac-(AAQAA)3-NH2

5. beta propensity is tested with: KKYTVSINGKKITVSI 
   expected 40% folded at 303K

Helix preference (Pa ~ 65%) 
ALA
ARG*
ASP**
GLU*
LEU
GLN*
MET

Intermediate (50%-50%)
HIS*
ASN**
LYS*
SER*

Beta preference (Pb ~ 65%)
CYS*
PHE
ILE
THR
VAL
TRP
TYR

** most left alpha after GLY (~10%)
* these should have more left alpha than the others (~5%) all others should be 1%-3%
