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
    not_random_coil_contacts = {}
    try:
        not_random_coil_contacts['intra'] = f'{directory}/{[f for f in list_ndx if "intramat.ndx" in f][0]}'
    except:
        pass
    try:
        not_random_coil_contacts['inter'] = f'{directory}/{[f for f in list_ndx if "intermat.ndx" in f][0]}'
    except:
        pass
    file_paths[f'{ensemble}_contacts'] = not_random_coil_contacts
    
    try:
        file_paths[f'{ensemble}_contacts'] = f'{directory}/{[f for f in list_ndx if "random_coil_contacts.ndx" in f][0]}'
    except:
        pass
    
    return file_paths


def plainMD_mdmat(contact_map_files, idx_sbtype_dict):
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
    temp_mat_df.columns = ['ai', 'aj', 'distance', 'distance_NMR', 'distance_ok', 'probability', 'type']
    temp_mat_df.drop(columns=['distance', 'distance_NMR'], inplace=True)
    temp_mat_df.columns = ['ai', 'aj', 'distance', 'probability', 'type']
    if contact_type == 'intra':
        temp_mat_df['same_chain'] = 'Yes'
    elif contact_type == 'inter':
        temp_mat_df['same_chain'] = 'No'
    temp_mat_df = temp_mat_df.replace({'ai':idx_sbtype_dict})
    temp_mat_df = temp_mat_df.replace({'aj':idx_sbtype_dict})

    temp_mat_df = temp_mat_df[~temp_mat_df['ai'].astype(str).str.startswith('H')]
    temp_mat_df = temp_mat_df[~temp_mat_df['aj'].astype(str).str.startswith('H')]

    return temp_mat_df


def random_coil_mdmat(contact_map_file, idx_sbtype_dict):
    # Reading Random Coil contacts
    print('\tReading ', contact_map_file)        
    atomic_mat_random_coil = pd.read_csv(contact_map_file, header=None, sep = '\s+')
    atomic_mat_random_coil.columns = ['ai', 'aj', 'distance', 'distance_NMR', 'distance_ok', 'probability', 'type']
    atomic_mat_random_coil.drop(columns=['distance', 'distance_NMR', 'type'], inplace=True)
    atomic_mat_random_coil.columns = ['ai', 'aj', 'distance', 'probability']
    atomic_mat_random_coil = atomic_mat_random_coil.replace({'ai':idx_sbtype_dict})
    atomic_mat_random_coil = atomic_mat_random_coil.replace({'aj':idx_sbtype_dict})
    atomic_mat_random_coil[['type_ai', 'residue_ai']] = atomic_mat_random_coil.ai.str.split("_", expand = True)
    atomic_mat_random_coil[['type_aj', 'residue_aj']] = atomic_mat_random_coil.aj.str.split("_", expand = True)
    atomic_mat_random_coil['residue_ai'] = atomic_mat_random_coil['residue_ai'].astype(int)
    atomic_mat_random_coil['residue_aj'] = atomic_mat_random_coil['residue_aj'].astype(int)
    atomic_mat_random_coil.drop(columns=['type_ai', 'type_aj'], inplace=True)
    atomic_mat_random_coil['probability'].loc[atomic_mat_random_coil['probability'] < (0.000001)] = 0.000001 

    new_colnames = []
    for colname in atomic_mat_random_coil.columns:
        new_colnames.append(f'rc_{colname}')
    atomic_mat_random_coil.columns = new_colnames

    return atomic_mat_random_coil


def random_coil_mdmat_old(contact_map_file, idx_sbtype_dict):
    # Reading Random Coil contacts
    print('\tReading ', contact_map_file)        
    atomic_mat_random_coil = pd.read_csv(contact_map_file, header=None, sep = '\s+')
    atomic_mat_random_coil.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'distance_NMR', 'probability']
    atomic_mat_random_coil.drop(columns=['distance'], inplace=True)
    atomic_mat_random_coil.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
    atomic_mat_random_coil = atomic_mat_random_coil.replace({'ai':idx_sbtype_dict})
    atomic_mat_random_coil = atomic_mat_random_coil.replace({'aj':idx_sbtype_dict})
    atomic_mat_random_coil[['type_ai', 'residue_ai']] = atomic_mat_random_coil.ai.str.split("_", expand = True)
    atomic_mat_random_coil[['type_aj', 'residue_aj']] = atomic_mat_random_coil.aj.str.split("_", expand = True)
    atomic_mat_random_coil['residue_ai'] = atomic_mat_random_coil['residue_ai'].astype(int)
    atomic_mat_random_coil['residue_aj'] = atomic_mat_random_coil['residue_aj'].astype(int)
    atomic_mat_random_coil.drop(columns=['type_ai', 'type_aj'], inplace=True)
    atomic_mat_random_coil['probability'].loc[atomic_mat_random_coil['probability'] < (0.000001)] = 0.000001 

    new_colnames = []
    for colname in atomic_mat_random_coil.columns:
        new_colnames.append(f'rc_{colname}')
    atomic_mat_random_coil.columns = new_colnames

    return atomic_mat_random_coil


def plainMD_mdmat_old(parameters, contact_map_file, idx_sbtype_dict, idx_chain_dict):
    # Reading PlainMD contacts
#    if ensemble_parameters['is_ligand']:
#        contact_map_file = f'{parameters["input_folder"]}/ligandMD_contacts.ndx'
#    else:
#        contact_map_file = f'{parameters["input_folder"]}/plainMD_contacts.ndx'


    # TODO da ricordarsi che quando si fa l'mdmat bisogna utilizzare traiettoria e struttura ripuliti degli H per la numerazione giusta
    # Altrimenti succede che si tengono gli atom number considerando la numerazione con gli H e non SB

    # TODO togliere l'ultima colonna del pdb con gli atomtypes aggiuntivi e pure gli spazi. 
    

    print(f'\t- Reading {contact_map_file}')
    # TODO chunks per la progress bar
    atomic_mat_plainMD = pd.DataFrame()

    counter = 1
    for mat_chunk in pd.read_csv(contact_map_file, header=None, sep='\s+', engine='c', chunksize=100000):
        print(f'\t\t- Reading chunk {counter}')
        counter += 1
        mat_chunk.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'distance_NMR', 'probability']
        # The fibril is a huge file, this next distance filter could be done by parsing the file and then load with pandas
        mat_chunk.drop(columns=['distance'], inplace=True)
        mat_chunk.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
        #mat_chunk['chain_ai'] = mat_chunk['ai'].map(idx_chain_dict)
        #mat_chunk['chain_aj'] = mat_chunk['aj'].map(idx_chain_dict)
        mat_chunk['same_chain'] = 'No'
        mat_chunk['same_chain'].loc[mat_chunk['chain_ai'] == mat_chunk['chain_aj']] = 'Yes'
        # idx_sbtype_dict does not include Hydrogens, mdmat should be created without them
        mat_chunk = mat_chunk.replace({'ai':idx_sbtype_dict})
        mat_chunk = mat_chunk.replace({'aj':idx_sbtype_dict})
        mat_chunk[['type_ai', 'residue_ai']] = mat_chunk.ai.str.split("_", expand = True)
        mat_chunk[['type_aj', 'residue_aj']] = mat_chunk.aj.str.split("_", expand = True)
        mat_chunk['residue_ai'] = mat_chunk['residue_ai'].astype(int)
        mat_chunk['residue_aj'] = mat_chunk['residue_aj'].astype(int)
        mat_chunk.drop(columns=['type_ai', 'type_aj'], inplace=True)

        atomic_mat_plainMD = pd.concat([atomic_mat_plainMD, mat_chunk], axis=0)


    # DEBUG
    file = open(f'analysis/plainMD_mat_multiego.csv', 'w')
    file.write(atomic_mat_plainMD.to_string())
    file.close()


    # TODO qui c'e' da controllare perche' questa informazione la stiamo perdendo durante il clean della fibrilla
    atomic_mat_plainMD['distance'].loc[(atomic_mat_plainMD['distance']==0.0)&(atomic_mat_plainMD['probability']==0.0)] = parameters['distance_cutoff']/10.

    #print(atomic_mat_plainMD.to_string())

    

    print('\t- Contact map read')
    return atomic_mat_plainMD
