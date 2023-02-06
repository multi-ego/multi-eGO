import chunk
import pandas as pd
import MDAnalysis as mda
import warnings
from os import listdir
from os.path import isfile, join


warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')
warnings.filterwarnings("ignore", category=DeprecationWarning) 
pd.options.mode.chained_assignment = None  # default='warn'


def find_files(ensemble, parameters, is_md=True):

    file_paths = {}
    directory = f'inputs/{parameters["protein"]}/{ensemble}'
    try:
        file_list = [f for f in listdir(directory) if isfile(join(directory, f))]
    except:
        raise Exception(f'Missing {directory}, check the name in input')
    try:
        file_paths[f'{ensemble}_topology'] = f'{directory}/{[f for f in file_list if ".top" in f][0]}'
    except:
        raise Exception(f'Missing topology in {directory}')
    
    list_ndx = [f for f in file_list if ".ndx" in f]
    if not list_ndx and parameters['egos'] != 'rc':
        raise Exception('\n\n\nIt is required a .ndx file to continue. \nTry egos=rc first.')        
    ensemble_contacts = {}
    try:
        ensemble_contacts['intra'] = f'{directory}/{[f for f in list_ndx if "intramat.ndx" in f][0]}'
    except:
        pass
    try:
        ensemble_contacts['inter'] = f'{directory}/{[f for f in list_ndx if "intermat.ndx" in f][0]}'
    except:
        pass
    try:
        ensemble_contacts['rc_intra'] = f'{directory}/{[f for f in list_ndx if "intra_random_coil_contacts.ndx" in f][0]}'
    except:
        pass
    try:
        ensemble_contacts['rc_inter'] = f'{directory}/{[f for f in list_ndx if "inter_random_coil_contacts.ndx" in f][0]}'
    except:
        pass
    
    file_paths[f'{ensemble}_contacts'] = ensemble_contacts

    return file_paths


def read_ensemble_mdmat_contacs(contact_map_files, idx_sbtype_dict):
    '''
    This function aggregates the two mdmat information #TODO finish
    '''
    print(f'\t- Reading {contact_map_files}')
    atomic_mat = pd.DataFrame()

    for contact_type, directory in contact_map_files.items():
        temp_mat_df = read_mdmat_dataframe(contact_type, directory, idx_sbtype_dict)
        atomic_mat = pd.concat([atomic_mat, temp_mat_df], axis=0)
    
    return atomic_mat


def read_mdmat_dataframe(contact_type, directory, idx_sbtype_dict):
    temp_mat_df = pd.read_csv(directory, header=None, sep = '\s+')
    temp_mat_df.columns = ['ai', 'aj', 'distance', 'distance_NMR', 'probability']
    temp_mat_df.drop(columns=['distance'], inplace=True)
    temp_mat_df.columns = ['ai', 'aj', 'distance', 'probability']

    if 'intra' in contact_type:
        temp_mat_df['same_chain'] = 'Yes'
    elif 'inter' in contact_type:
        temp_mat_df['same_chain'] = 'No'
    temp_mat_df = temp_mat_df.replace({'ai':idx_sbtype_dict})
    temp_mat_df = temp_mat_df.replace({'aj':idx_sbtype_dict})

    temp_mat_df = temp_mat_df[~temp_mat_df['ai'].astype(str).str.startswith('H')]
    temp_mat_df = temp_mat_df[~temp_mat_df['aj'].astype(str).str.startswith('H')]
    temp_mat_df['distance'].loc[(temp_mat_df['probability'] < (0.000001))&(temp_mat_df['distance'] == 0.)] = 0.550000 
    temp_mat_df['probability'].loc[temp_mat_df['probability'] < (0.000001)] = 0.000001 

    # Questo è per tenerlo il più simile a prima, nel MD non rinomino le chain ma nel rc si.
    if 'rc' in contact_type:
        temp_mat_df[['type_ai', 'residue_ai']] = temp_mat_df.ai.str.split("_", expand = True)
        temp_mat_df[['type_aj', 'residue_aj']] = temp_mat_df.aj.str.split("_", expand = True)
        temp_mat_df['residue_ai'] = temp_mat_df['residue_ai'].astype(int)
        temp_mat_df['residue_aj'] = temp_mat_df['residue_aj'].astype(int)
        temp_mat_df.drop(columns=['type_ai', 'type_aj'], inplace=True)
        #temp_mat_df['probability'].loc[temp_mat_df['probability'] < (0.000001)] = 0.000001 
        
        new_colnames = []
        for colname in temp_mat_df.columns:
            new_colnames.append(f'rc_{colname}')
        temp_mat_df.columns = new_colnames

    return temp_mat_df
