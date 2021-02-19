import pandas as pd
from protein_configuration import fibril_chain_length, fibril_residue_offset, fibril_atom_number, fibril_atom_offset, protein
import MDAnalysis as mda

from gromologist import Top, Pdb


# This script includes all the functions used to read the input files.


# The following three functions are used to read the atoms information from the topology

def read_pep_atoms():
    directory = 'input_%s/pep_atoms' % (protein)
    pep_atoms = pd.read_csv(directory, sep = '\\s+', header = None)
    # header=0 because it counts ; as a column. Therefore, I don't care about the header and I'll recreate one from
    # scratch.
    pep_atoms["charge"] = ""
    # addition of a charge empty column
    pep_atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", "charge"]
    # I've created the header with the correct column names
    # The following function is necessary to associate the resID with the residue number.
    # Changing the atom type column -> patoms column 5 + 3
    pep_atoms["type"] = pep_atoms["atom"].apply(str) + '_' + pep_atoms["resnr"].apply(str)
    # Afterwards the function make_pep_atp_dict will create a dictionary based on this column
    
    #print(pep_atoms.to_string())
    
    return pep_atoms


def read_fib_atoms():
    directory = 'input_%s/fib_atoms' % (protein)
    fib_atoms = pd.read_csv(directory, sep = '\\s+', header = None)
    # header=0 because it counts ; as a column.
    # Therefore, I don't care about the header and I'll recreate one from scratch.
    fib_atoms["charge"] = ""
    # addition of a charge empty column
    fib_atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", "charge"]
    
    # THIS IS IMPORTANT TO CHANGE WHEN USING A DIFFERENT FIBRIL
    fib_atoms['resnr'] = fib_atoms['resnr'] % fibril_chain_length #64 for B2m and 11 for TTR #flag
    fib_atoms['resnr'] = fib_atoms['resnr'].replace(0, fibril_chain_length) #flag

    # In B2m model, there are not N or C terminal in the model, therefore it is necessary to change the atomid and resid
    # Unfortunately SMOG always renumber everything and so i renumber using python
    # The first residue of the fibril is 23 and the first atom is 179
    # The atomID is not renumbered because those values will be replaced using the fib_atoms with the atomtype.

    #fib_atoms['; nr'] = fib_atoms['; nr']+178 #flag
    #fib_atoms['cgnr'] = fib_atoms['cgnr']+178 #flag
    fib_atoms['resnr'] = fib_atoms['resnr'] + fibril_residue_offset #flag

    # Likewise, the same procedure will be applied also in dihedrals and pairs
    # STARE ATTENTO ALLA RINUMERAZIONE! -> DALLA FIBRILLA NON PARTONO DA 1 MA DA 23

    # I've created the header with the correct column names
    fib_atoms["type"] = fib_atoms["atom"].apply(str) + '_' + fib_atoms["resnr"].apply(str)
    return fib_atoms

    
def read_gro_atoms(): # GRETA TO KEEP
    # Reading the atoms section from gromos topology
    # Requires a manual clean to delete all the comment lines
    directory = 'input_%s/pep_gro_atoms' % (protein)
    gro_atoms = pd.read_csv(directory, sep = "\\s+", header = None)
    gro_atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", 'charge', 'mass']
    gro_atoms["atom_nmr"] = gro_atoms["atom"].apply(str) + '_' + gro_atoms["resnr"].apply(str)
    gro_atoms['res_atom'] = gro_atoms['residue'] + '_' + gro_atoms['atom']
    return gro_atoms

# Those two functions are used to read the dihedrals from SMOG topology

def read_pep_dihedrals():
    # Reading the peptide dihedrals
    directory = 'input_%s/pep_dihedrals' % (protein)
    pep_dihedrals = pd.read_csv(directory, sep = "\\s+", header = None)
    pep_dihedrals.columns = [";ai", "aj", "ak", "al", "func", "phi0", "Kd", "mult"]
    pep_dihedrals['mult'] = pep_dihedrals['mult'].fillna(value = '')
    # dihedrals = dihedrals.replace(np.nan, '', regex=True)
    # NaN are replaced with and empty line.
    return pep_dihedrals


def read_fib_dihedrals():
    # Reading the fib_atomstide dihedrals
    directory = 'input_%s/fib_dihedrals' % (protein)
    fib_dihedrals = pd.read_csv(directory, sep = "\\s+", header = None)
    fib_dihedrals.columns = [";ai", "aj", "ak", "al", "func", "phi0", "Kd", "mult"]
    
    # NaN are replaced with and empty line.
    fib_dihedrals['mult'] = fib_dihedrals['mult'].fillna(value = '')

    # AtomID renumber to match the native structure because I use atomnumber
    fib_dihedrals[';ai'] = fib_dihedrals[';ai']-fibril_atom_number+fibril_atom_offset #flag
    fib_dihedrals['aj'] = fib_dihedrals['aj']-fibril_atom_number+fibril_atom_offset #flag
    fib_dihedrals['ak'] = fib_dihedrals['ak']-fibril_atom_number+fibril_atom_offset #flag
    fib_dihedrals['al'] = fib_dihedrals['al']-fibril_atom_number+fibril_atom_offset #flag
    return fib_dihedrals

def read_pep_pairs():
    # Reading the peptide pairs
    directory = 'input_%s/pep_pairs' % (protein)
    pep_pairs = pd.read_csv(directory, sep = "\\s+", header = None)
    pep_pairs.columns = [";ai", "aj", "type", "A", "B"]
    return pep_pairs

def read_fib_pairs():
    # Reading the fib_atomstide pairs
    directory = 'input_%s/fib_pairs' %(protein)
    fib_pairs = pd.read_csv(directory, sep = "\\s+", header = None)
    fib_pairs.columns = [";ai", "aj", "type", "A", "B"]
    # The atomID is not renumbered because those values will be replaced using the fib_atoms with the atomtype.
    #fib_pairs[';ai'] = fib_pairs[';ai']+178 #flag
    #fib_pairs['aj'] = fib_pairs['aj']+178 #flag

    return fib_pairs


##############
# GRETA
##############


def read_pdbs():
    
    native_directory = 'GRETA/native_%s/native.pdb' %(protein)
    fibril_directory = 'GRETA/fibril_%s/conf.pdb' %(protein)
    native_pdb = mda.Universe(native_directory, guess_bonds = True)
    fibril_pdb = mda.Universe(fibril_directory, guess_bonds = True)

    return native_pdb, fibril_pdb

def read_top():
    native_directory = 'GRETA/native_%s/topol.top' %(protein)
    #native_directory = 'gromologist/examples/01_pentapeptide/topol.top'
    native_pdb = 'GRETA/native_%s/native.pdb' %(protein)
    native_topology = Top(native_directory, gmx_dir='/home/emanuele/MAGROS/GRETA', pdb=native_pdb)
    
    return native_topology

def read_gro_bonds():
    native_directory = 'GRETA/native_%s/gro_bonds' %(protein)
    native_bonds = pd.read_csv(native_directory, sep = "\\s+", header = None)
    native_bonds.columns = ["ai", "aj", "func", "def"]

    return native_bonds

def read_gro_angles():
    native_directory = 'GRETA/native_%s/gro_angles' %(protein)
    native_angles = pd.read_csv(native_directory, sep = "\\s+", header = None)
    native_angles.columns = ["ai", "aj", "ak", "func", "def"]

    return native_angles

def read_gro_dihedrals():
    native_directory = 'GRETA/native_%s/gro_dihedrals' %(protein)
    native_dihedrals = pd.read_csv(native_directory, sep = "\\s+", header = None)
    native_dihedrals.columns = ["ai", "aj", "ak", "al", "func", "def"]

    return native_dihedrals

def read_gro_impropers():
    native_directory = 'GRETA/native_%s/gro_impropers' %(protein)
    native_dihedrals = pd.read_csv(native_directory, sep = "\\s+", header = None)
    native_dihedrals.columns = ["ai", "aj", "ak", "al", "func", "def"]

    return native_dihedrals