from time import localtime, strftime
import pandas as pd


pd.set_option('display.colheader_justify', 'right')

def get_name(parameters):
    if parameters.egos=="rc":
        name = f'{parameters.protein}_{parameters.egos}'
    else:
        name = f'{parameters.protein}_{parameters.egos}_e{parameters.epsilon}_{parameters.inter_epsilon}'

    return name


def dataframe_to_write(df):
    if df.empty:
        # TODO insert and improve the following warning
        print('A topology parameter is empty. Check the reference topology.')
        return '; The following parameters where not parametrized on multi-eGO.\n; If this is not expected, check the reference topology.'
    else:
        df.rename(columns={df.columns[0]: f"; {df.columns[0]}"}, inplace = True)
        return df.to_string(index=False)


def make_header(parameters):
    now = strftime("%d-%m-%Y %H:%M", localtime())
    
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
        if molecule_type_dict[molecule] == 'protein':
            pairs = lj_14[molecule]
            pairs.insert(5, ';', ';')
            pairs['c6'] = pairs["c6"].map(lambda x:'{:.6e}'.format(x))
            pairs['c12'] = pairs["c12"].map(lambda x:'{:.6e}'.format(x))
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
        atom_selection_dataframe = topology_dataframe.loc[topology_dataframe['molecule_name'] == molecule][['number', 'sb_type', 'resnum', 'resname', 'name', 'cgnr']].copy()
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
        #file.write(f'[ pairs ]\n')
        #file.write(dataframe_to_write(pairs))
        file.write(f'[ exclusions ]\n')
        file.write(dataframe_to_write(exclusions))

    footer = f'''

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

[ system ]
{parameters.protein}

[ molecules ]
; Compound #mols
'''

    file.write(footer)
    for molecule in molecule_footer:
        file.write(f'{molecule}\t\t\t1\n')
    
    file.close()


def write_nonbonded(topology_dataframe, lj_potential, parameters, output_folder):
    #pd.set_option('display.colheader_justify', 'right')
    header = make_header(vars(parameters))
    file = open(f'{output_folder}/ffnonbonded.itp', 'w')
    file.write(header)
    file.write('\n[ atomtypes ]\n')
    atomtypes = topology_dataframe[['sb_type', 'atomic_number', 'mass', 'charge', 'ptype', 'c6', 'c12']].copy()
    atomtypes['c6'] = atomtypes['c6'].map(lambda x:'{:.6e}'.format(x))
    atomtypes['c12'] = atomtypes['c12'].map(lambda x:'{:.6e}'.format(x))
    file.write(dataframe_to_write(atomtypes))
    if not lj_potential.empty:
        file.write("\n\n[ nonbond_params ]\n")
        lj_potential['c6'] = lj_potential['c6'].map(lambda x:'{:.6e}'.format(x))
        lj_potential['c12'] = lj_potential['c12'].map(lambda x:'{:.6e}'.format(x))
        lj_potential.insert(5, ';', ';')
        lj_potential.drop(columns = ['rc_distance', 'molecule_name_ai', 'molecule_name_aj', 'rc_molecule_name_ai', 'rc_molecule_name_aj', 'rc_ai', 'rc_aj', 'rc_same_chain', 'rc_source'], inplace=True)
        file.write(dataframe_to_write(lj_potential))

