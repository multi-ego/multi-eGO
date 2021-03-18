import pandas as pd
from protein_configuration import protein
import MDAnalysis as mda
from gromologist import Top

def read_pdbs():
    
    native_directory = 'native_%s/native.pdb' %(protein)
    fibril_directory = 'fibril_%s/conf.pdb' %(protein)
    native_pdb = mda.Universe(native_directory, guess_bonds = True)
    fibril_pdb = mda.Universe(fibril_directory, guess_bonds = True)

    return native_pdb, fibril_pdb

def read_top():
    native_directory = 'native_%s/topol.top' %(protein)
    #native_directory = 'gromologist/examples/01_pentapeptide/topol.top'
    native_pdb = 'native_%s/native.pdb' %(protein)
    native_topology = Top(native_directory, gmx_dir='/home/emanuele/MAGROS', pdb=native_pdb)
    
    return native_topology