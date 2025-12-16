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

3. 3J coupling are matched 
Add table:


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
