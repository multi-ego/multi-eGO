import pandas as pd
import MDAnalysis as mda
from gromologist import Top
import warnings

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

def read_pdbs(protein, greta_to_keep):
    if greta_to_keep == 'native':
        native_directory = 'inputs/native_%s/native.pdb' %(protein)
        native_pdb = mda.Universe(native_directory, guess_bonds = True)

    elif greta_to_keep == 'fibril':
        fibril_directory = 'inputs/fibril_%s/fibril.pdb' %(protein)
        fibril_pdb = mda.Universe(fibril_directory, guess_bonds = True)

    elif greta_to_keep == 'all':
        native_directory = 'inputs/native_%s/native.pdb' %(protein)
        fibril_directory = 'inputs/fibril_%s/fibril.pdb' %(protein)
        native_pdb = mda.Universe(native_directory, guess_bonds = True)
        fibril_pdb = mda.Universe(fibril_directory, guess_bonds = True)

    else:
        print("I dont' understand --build-from=",greta_to_keep)
        exit()

    return native_pdb, fibril_pdb

def read_top(protein):
    native_directory = 'inputs/native_%s/topol.top' %(protein)
    native_pdb = 'inputs/native_%s/native.pdb' %(protein)
    native_topology = Top(native_directory, gmx_dir='/home/emanuele/MAGROS', pdb=native_pdb)
    
    return native_topology

def read_native_pairs(protein, distance_residue, distance_cutoff):
    native_directory = 'inputs/native_%s/monomer_pairs_md_ex%s_co%s.txt' %(protein, distance_residue, distance_cutoff)
    native_pairs = pd.read_csv(native_directory, sep = '\\s+', header = None)
    native_pairs.columns = ['ai', 'aj', 'counts', 'ratio', 'distance']

    return native_pairs
