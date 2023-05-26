import pandas as pd
import time
import glob
import os

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

def write_nonbonded(topology_dataframe, lj_potential, parameters, output_folder):
    # pd.set_option('display.colheader_justify', 'right')
    header = make_header(vars(parameters))
    file = open(f'{output_folder}/ffnonbonded.itp', 'w')
    file.write(header)
    file.write('\n[ atomtypes ]\n')
    atomtypes = topology_dataframe[['sb_type', 'atomic_number', 'mass', 'charge', 'ptype', 'c6', 'c12']].copy()
    atomtypes['c6'] = atomtypes['c6'].map(lambda x: '{:.6e}'.format(x))
    atomtypes['c12'] = atomtypes['c12'].map(lambda x: '{:.6e}'.format(x))
    file.write(dataframe_to_write(atomtypes))

    if not lj_potential.empty:
        file.write("\n\n[ nonbond_params ]\n")
        lj_potential['c6'] = lj_potential['c6'].map(lambda x: '{:.6e}'.format(x))
        lj_potential['c12'] = lj_potential['c12'].map(lambda x: '{:.6e}'.format(x))
        lj_potential.insert(5, ';', ';')
        lj_potential.drop(columns= ['molecule_name_ai', 'molecule_name_aj'], inplace=True)
        file.write(dataframe_to_write(lj_potential))

    file.close()

def write_model(meGO_ensemble, meGO_LJ_potential, meGO_LJ_14, parameters, output_dir):
    print('- Writing Multi-eGO model')
    write_topology(meGO_ensemble['reference_topology_dataframe'], meGO_ensemble['molecule_type_dict'], meGO_ensemble['meGO_bonded_interactions'], meGO_LJ_14, parameters, output_dir)
    write_nonbonded(meGO_ensemble['reference_topology_dataframe'], meGO_LJ_potential, parameters, output_dir)

    print('\n- The model is baked with the following parameters:\n')
    for argument, value in vars(parameters).items():
        if type(value) is list:
            print('\t- {:<20} = {:<20}'.format(argument, ", ".join(value)))
        elif type(value) is not str:
            print('\t- {:<20} = {:<20}'.format(argument, str(value)))
        else:
            print('\t- {:<20} = {:<20}'.format(argument, value))
    if parameters.egos != 'rc':
        print(f'''
        - LJ parameterization completed with a total of {len(meGO_LJ_potential)} contacts.
        - The average epsilon is {meGO_LJ_potential['epsilon'].loc[meGO_LJ_potential['epsilon']>0.].mean():{5}.{3}}
        - The maximum epsilon is {meGO_LJ_potential['epsilon'].max():{5}.{3}}
        - The maximum sigma is {meGO_LJ_potential['sigma'].loc[meGO_LJ_potential['epsilon']>0.].max():{5}.{3}}, suggested cut-off at {2.5*meGO_LJ_potential['sigma'].loc[meGO_LJ_potential['epsilon']>0.].max():{4}.{3}} nm
        ''')
    print(
        f'\nAnd it can be found in the following folder:\n{output_dir}')
    print('\nNessuno è più basito, nessuno è più sorpreso. Ognuno di voi ha capito tutto.\nCarlo is happy!\t\^o^/\n')

def dataframe_to_write(df):
    if df.empty:
        # TODO insert and improve the following warning
        print('A topology parameter is empty. Check the reference topology.')
        return '; The following parameters where not parametrized on multi-eGO.\n; If this is not expected, check the reference topology.'
    else:
        df.rename(columns={df.columns[0]: f"; {df.columns[0]}"}, inplace= True)
        return df.to_string(index=False)

def make_header(parameters):
    now = time.strftime("%d-%m-%Y %H:%M", time.localtime())

    header = f'''
; Multi-eGO force field provided by Emanuele Scalone and Carlo Camilloni at Camilloni Lab
; Created on the {now}
; With the following parameters:
'''
    for parameter, value in parameters.items():
        if type(value) is list:
            header += ';\t- {:<15} = {:<20}\n'.format(parameter, ", ".join(value))
        elif not value:
            value = ''
            header += ';\t- {:<15} = {:<20}\n'.format(parameter, ", ".join(value))
        else:
            header += ';\t- {:<15} = {:<20}\n'.format(parameter, value)
    return header

def write_topology(topology_dataframe, molecule_type_dict, bonded_interactions_dict, lj_14, parameters, output_folder):
    molecule_footer = []
    header = make_header(vars(parameters))
    file = open(f'{output_folder}/topol_GRETA.top', 'w')
    header += '''
; Include forcefield parameters
#include "multi-ego-basic.ff/forcefield.itp"
'''

    file.write(header)
    for molecule, bonded_interactions in bonded_interactions_dict.items():
        # TODO here I defined an empty exclusion. In the case we are not reading a protein topology, the exclusion part is not read and needed to be added in someway.
        # Hence, an empty exclusion gives error. Here I define an empty variable so it does not gets stuck
        exclusions = pd.DataFrame(columns=['ai', 'aj'])
        # TODO here only proteins have custom pairs and exclusions. Nucleic acids and others will use the one in topol.top used as reference
        # if molecule_type_dict[molecule] == 'protein':
        pairs = lj_14[molecule]
        pairs.insert(5, ';', ';')
        pairs['c6'] = pairs["c6"].map(lambda x: '{:.6e}'.format(x))
        pairs['c12'] = pairs["c12"].map(lambda x: '{:.6e}'.format(x))
        bonded_interactions_dict[molecule]['pairs'] = pairs
        exclusions = pairs[['ai', 'aj']].copy()

        molecule_footer.append(molecule)
        molecule_header = f'''
[ moleculetype ]
; Name\tnrexcl
{molecule}\t\t\t3

'''

        file.write(molecule_header)
        file.write('[ atoms ]\n')
        atom_selection_dataframe = topology_dataframe.loc[topology_dataframe['molecule_name'] == molecule][[
            'number', 'sb_type', 'resnum', 'resname', 'name', 'cgnr']].copy()
        file.write(f'{dataframe_to_write(atom_selection_dataframe)}\n\n')
        # Here are written bonds, angles, dihedrals and impropers
        for bonded_type, interactions in bonded_interactions.items():
            if interactions.empty:
                continue
            else:
                if bonded_type == 'impropers':
                    file.write(f'[ dihedrals ]\n')
                else:
                    file.write(f'[ {bonded_type} ]\n')
                file.write(dataframe_to_write(interactions))
                file.write('\n\n')
        # Here are written pairs and exclusions
        # file.write(f'[ pairs ]\n')
        # file.write(dataframe_to_write(pairs))
        file.write(f'[ exclusions ]\n')
        file.write(dataframe_to_write(exclusions))

    footer = f'''

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

[ system ]
{parameters.system}

[ molecules ]
; Compound #mols
'''

    file.write(footer)
    for molecule in molecule_footer:
        file.write(f'{molecule}\t\t\t1\n')

    file.close()

def get_name(parameters):
    if parameters.egos =="rc":
        name = f'{parameters.system}_{parameters.egos}'
    else:
        name = f'{parameters.system}_{parameters.egos}_e{parameters.epsilon}_{parameters.inter_epsilon}'
    return name

def create_output_directories(parameters):
    if parameters.egos == 'rc':
        name = f'{parameters.system}_{parameters.egos}'
    else: name = f'{parameters.system}_{parameters.egos}_e{parameters.epsilon}_{parameters.inter_epsilon}'
    output_folder = f'outputs/{name}'
    
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    if os.path.isfile(f'{output_folder}/ffnonbonded.itp'): os.remove(f'{output_folder}/ffnonbonded.itp')
    if os.path.isfile(f'{output_folder}/topol_GRETA.top'):
        os.remove(f'{output_folder}/topol_GRETA.top')

    return output_folder

def check_files_existence(ego, protein, md_ensembles):
    for ensemble in md_ensembles:
        ensemble = f'inputs/{protein}/{ensemble}'
        if not os.path.exists(ensemble):
            raise FileNotFoundError(f"Folder {ensemble}/ does not exist.")
        else:
            top_files = glob.glob(f'{ensemble}/*.top')
            if not top_files:
                raise FileNotFoundError(f"No .top files found in {ensemble}/")
            ndx_files = glob.glob(f'{ensemble}/*.ndx')
            if not ndx_files and not ego =="rc":
                raise FileNotFoundError(f"No .ndx files found in {ensemble}/")
