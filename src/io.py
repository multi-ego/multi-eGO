import pandas as pd

def read_molecular_contacts(path):
    '''
    This function add column headers and add a column based whether the contacts are intramolecular or intermolecular
    '''

    print('\t-', f"Reading {path}")
    contact_matrix = pd.read_csv(path, header=None, sep='\s+')
    contact_matrix.columns = ['molecule_number_ai', 'ai', 'molecule_number_aj', 'aj', 'distance_m', 'distance_NMR', 'distance_14', 'probability', 'cutoff', 'unimod']
    contact_matrix.columns = ['molecule_number_ai', 'ai', 'molecule_number_aj', 'aj', 'distance_m', 'distance', 'distance_14', 'probability', 'cutoff', 'unimod']
    contact_matrix.drop('unimod', axis=1, inplace=True)
    contact_matrix['molecule_number_ai'] = contact_matrix['molecule_number_ai'].astype(str)
    contact_matrix['ai'] = contact_matrix['ai'].astype(str)
    contact_matrix['molecule_number_aj'] = contact_matrix['molecule_number_aj'].astype(str)
    contact_matrix['aj'] = contact_matrix['aj'].astype(str)

    return contact_matrix
