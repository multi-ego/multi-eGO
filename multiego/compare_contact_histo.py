#! python

import os
import pandas as pd
import numpy as np
import parmed as pmd
import multiprocessing
import argparse
import itertools
import time
import sys
import warnings

import resources.type_definitions as type_definitions

d = { type_definitions.gromos_atp.name[i] : type_definitions.gromos_atp.c12[i] for i in range(len(type_definitions.gromos_atp.name))}

def run_intra_(arguments):
    '''
    Preforms the main routine of the histogram analysis to obtain the intra- and intermat files.
    Is used in combination with multiprocessing to speed up the calculations.

    Parameters
    ----------
    arguments : dict
        Contains all the command-line parsed arguments

    Returns
    -------
    out_path : str
        Path to the temporary file which contains a partial pd.DataFrame with the analyzed data 
    '''
    (args, protein_ref_indices_i, protein_ref_indices_j, original_size_j, c12_cutoff, mi, mj, frac_target_list) = arguments
    process = multiprocessing.current_process()
    columns = ['mi', 'ai', 'mj', 'aj', 'dist', 'c12dist', 'hdist', 'p', 'cutoff']
    df = pd.DataFrame(columns=columns)
    # We do not consider old histograms
    frac_target_list = [ x for x in frac_target_list if x[0]!="#" and x[-1]!="#" ]
    for i, ref_f in enumerate(frac_target_list):
        results_df = pd.DataFrame()
        ai = ref_f.split('.')[-2].split('_')[-1]
        
        if True:
            all_ai = [ ai for _ in range(1, original_size_j+1) ]
            range_list = [ str(x) for x in range(1, original_size_j+1) ]

            results_df['ai'] = np.array(all_ai).astype(int)
            results_df['aj'] = np.array(range_list).astype(int)

            results_df['mi'] = mi
            results_df['mj'] = mj
            results_df['dist'] = 0
            results_df['c12dist'] = 0
            results_df['hdist'] = 0
            results_df['p'] = 0
            results_df['cutoff'] = 0
            
        if np.isin(int(ai), protein_ref_indices_i):
            cut_i = np.where(protein_ref_indices_i == int(ai))[0][0]

            # column mapping
            ref_f = f'{args.histo}/{ref_f}'
            ref_df = pd.read_csv(ref_f, header=None, sep='\s+', usecols=[ 0, *protein_ref_indices_j ])
            ref_df_columns = ['distance', *[ str(x) for x in protein_ref_indices_j ]]
            ref_df.columns = ref_df_columns
            ref_df.set_index('distance', inplace=True)
            ref_df.loc[len(ref_df)] = c12_cutoff[cut_i]

            # calculate data
            dist = ref_df.apply(lambda x: weighted_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            c12dist = ref_df.apply(lambda x: c12_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            hdist = ref_df.apply(lambda x: weighted_havg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            p = ref_df.apply(lambda x: calculate_probability(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values

            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'dist'] = dist
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'c12dist'] = c12dist
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'hdist'] = hdist
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'p'] = p
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'cutoff'] = c12_cutoff[cut_i]
    
        df = pd.concat([df, results_df])
        df = df.sort_values(by = ['p', 'c12dist', 'dist'], ascending=True)

    df.fillna(0)
    out_path = f'mat_{process.pid}_t{time.time()}.part'
    df.to_csv(out_path, index=False)

    return out_path

def run_inter_(arguments):
    '''
    Preforms the main routine of the histogram analysis to obtain the intra- and intermat files.
    Is used in combination with multiprocessing to speed up the calculations.

    Parameters
    ----------
    arguments : dict
        Contains all the command-line parsed arguments

    Returns
    -------
    out_path : str
        Path to the temporary file which contains a partial pd.DataFrame with the analyzed data 
    '''
    (args, protein_ref_indices_i, protein_ref_indices_j, original_size_j, c12_cutoff, mi, mj, frac_target_list) = arguments
    # frac_target_list, frac_cumulative_probabilities = [ x[0] for x in frac_target_list ], [ x[1] for x in frac_target_list ]
    process = multiprocessing.current_process()
    columns = ['mi', 'ai', 'mj', 'aj', 'dist', 'c12dist', 'hdist', 'p', 'cutoff']
    df = pd.DataFrame(columns=columns)
    # We do not consider old histograms
    frac_target_list = [ x for x in frac_target_list if x[0]!="#" and x[-1]!="#" ]
    for i, ref_f in enumerate(frac_target_list):
        results_df = pd.DataFrame()
        ai = ref_f.split('.')[-2].split('_')[-1]
        
        if True:
            all_ai = [ ai for _ in range(1, original_size_j+1) ]
            range_list = [ str(x) for x in range(1, original_size_j+1) ]

            results_df['ai'] = np.array(all_ai).astype(int)
            results_df['aj'] = np.array(range_list).astype(int)

            results_df['mi'] = mi
            results_df['mj'] = mj
            results_df['dist'] = 0
            results_df['c12dist'] = 0
            results_df['hdist'] = 0
            results_df['p'] = 0
            results_df['cutoff'] = 0
            
        if np.isin(int(ai), protein_ref_indices_i):
            cut_i = np.where(protein_ref_indices_i == int(ai))[0][0]

            # column mapping
            ref_f = f'{args.histo}/{ref_f}'
            ref_df = pd.read_csv(ref_f, header=None, sep='\s+', usecols=[ 0, *protein_ref_indices_j ])
            ref_df_columns = ['distance', *[ str(x) for x in protein_ref_indices_j ]]
            ref_df.columns = ref_df_columns
            ref_df.set_index('distance', inplace=True)
            ref_df.loc[len(ref_df)] = c12_cutoff[cut_i]

            # repeat for cumulative
            c_ref_f = ref_f.replace('inter_mol_', 'inter_mol_c_')
            c_ref_df = pd.read_csv(c_ref_f, header=None, sep='\s+', usecols=[ 0, *protein_ref_indices_j ])
            c_ref_df_columns = ['distance', *[ str(x) for x in protein_ref_indices_j ]]
            c_ref_df.columns = c_ref_df_columns
            c_ref_df.set_index('distance', inplace=True)
            c_ref_df.loc[len(c_ref_df)] = c12_cutoff[cut_i]

            # calculate data
            dist = ref_df.apply(lambda x: weighted_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            c12dist = ref_df.apply(lambda x: c12_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            hdist = ref_df.apply(lambda x: weighted_havg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            p = c_ref_df.apply(lambda x: get_cumulative_probability(c_ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values

            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'dist'] = dist
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'c12dist'] = c12dist
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'hdist'] = hdist
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'p'] = p
            results_df.loc[results_df['aj'].isin(protein_ref_indices_j), 'cutoff'] = c12_cutoff[cut_i]
    
        df = pd.concat([df, results_df])

        df = df.sort_values(by = ['p', 'c12dist', 'dist'], ascending=True)

    df.fillna(0)
    out_path = f'mat_{process.pid}_t{time.time()}.part'
    df.to_csv(out_path, index=False)

    return out_path

def read_topologies(mego_top, target_top):
    '''
    Reads the input topologies using parmed. Ignores warnings to prevent printing
    of GromacsWarnings regarding 1-4 interactions commonly seen when using
    parmed in combination with multi-eGO topologies.

    Parameters
    ----------
    mego_top : str
        Path to the multi-eGO topology obtained from gmx pdb2gmx with multi-ego-basic force fields
    target_top : str
        Path to the toplogy of the system on which the analysis is to be performed
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        topology_mego = pmd.load_file(mego_top)
        topology_ref = pmd.load_file(target_top)

    n_mol=len(list(topology_mego.molecules.keys()))
    mol_names=list(topology_mego.molecules.keys())
    mol_list=np.arange(1,n_mol+1,1)

    return topology_mego, topology_ref, n_mol, mol_names, mol_list

def map_if_exists(atom_name):
    '''
    Maps an atom name to a multi-eGO atom name if possible

    Parameters
    ----------
    atom_name : str
        The atom name with which to attempt the mapping
    
    Return
    ------
    atom_name : str
        Mapped atom name. Equal to the input if mapping was not possible
    '''
    if atom_name in type_definitions.from_ff_to_multiego.keys(): return type_definitions.from_ff_to_multiego[atom_name]
    else: return atom_name

def hallfunction(values, weights):
    '''
    A substitute to the allfunction without considering the cutoff.
    Does not change the data except for the removal of the cutoff from the array.

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Return
    ------
    norm : float
        The normalization constant
    v : np.array
        The unchanged x values of the histogram
    w : np.array
        The unchanged weights of the histogram
    '''
    v = values[:-1]
    w = weights[:-1]
    norm = np.sum(w)
    return norm, v, w

def allfunction(values, weights):
    '''
    TODO rename pls

    Preprocesses arrays (histograms) to allow for proper analysis. Last values are removed from the arrays
    and should correspond to the respective cutoff for the histogram. The histograms are truncated
    according to the cutoff.

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Returns
    -------
    cutoff : float
        The cutoff which is deduced by reading the last value of the weights array
    i : int
        The index at which the cutoff is greter or equal than the values array
    norm : float
        The new normalization constant after truncation
    v : np.array
        The truncated x values of the histogram according to the cutoff
    w : np.array
        The truncated weights of the histogram according to the cutoff
    '''
    v = values[:-1]
    cutoff = weights[len(weights)-1]
    w = weights[:-1]
    i = np.where(v <= cutoff)
    if not np.any(i): return 0,0,0,0,0  # check if empty
    i = i[0]
    w = w[i]
    v = v[i]
    norm = np.sum(w)
    i = i[-1]
    return cutoff, i, norm, v, w

def zero_callback(values, weights):
    '''
    A zero callback doing nothing but returning the normalization constant, values and weights.
    The first two return values so the function can be switched with allfunction.

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Returns
    -------
    None : None
        None
    None : None
        None
    np.sum(weights) : float
        The normalization constant
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights  
    '''
    return None, None, np.sum(weights), values, weights

def get_cumulative_probability(values, weights, callback=allfunction):
    cutoff, i, norm, v, w = callback(values, weights)
    return weights[i]

def weighted_havg(values, weights, callback=hallfunction):
    norm, v, w = callback(values, weights)
    if norm == 0.: return 0
    return np.sum(v * w) / norm

def weighted_avg(values, weights, callback=allfunction):
    '''
    Calculates the weighted average as 1 / n * \sum_i^n v[i] * w[i]

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights
    callback : function
        Preprocesses the data before going in to the analysis

    Returns
    -------
    The weighted average
    '''
    cutoff, i, norm, v, w = callback(values, weights)
    if norm == 0.: return 0
    return np.sum(v * w) / norm


def c12_avg(values, weights, callback=allfunction):
    '''
    Calculates the c12 averaging of a histogram as 1 / ( (\sum_i^n w[i] * (1 / x[i])^12 ) / norm )^(1/12)

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights
    callback : function
        Preprocesses the data before going in to the analysis

    Returns
    -------
    The c12 average
    '''
    cutoff, i, norm, v, w = callback(values, weights)
    if np.sum(w)==0: return 0
    r = np.where(w > 0.)
    v = v[r[0][0]:v.size]
    w = w[r[0][0]:w.size]
    
    #exp_aver = (1./0.1)/np.log(np.sum(w*np.exp(1./v/0.1))/norm)
    res = np.maximum(cutoff/4.5, 0.1)
    exp_aver = (1./res)/np.log(np.sum(w*np.exp(1./v/res))/norm)

    return exp_aver


def warning_cutoff_histo(cutoff, max_adaptive_cutoff):
    '''
    Prints warning if the histogram cutoff is smaller as the maximum adaptive cutoff.

    Parameters
    ----------
    cutoff : float
        The cutoff of the histogram calculations. Parsed from the command-line in the standard programm.
    max_adaptive_cutoff : float
        The maximum adaptive cutoff calculated from the LJ c12 parameters.
    '''    
    print(f"""
    #############################
    
    -------------------
    WARNING
    -------------------

    Found an adaptive cutoff greater then the cutoff used to generate the histogram:
    histogram cutoff = {cutoff}
    maximum adaptive cutoff = {max_adaptive_cutoff}

    Be careful!. This could create errors.
    If this is not wanted, please recalculate the histograms setting the cutoff to at least cutoff={max_adaptive_cutoff}
    
    #############################
    """)


def calculate_intra_probabilities(args):
    '''
    Starts the main routine for calculating the intramat by:
     - reading the topologies
     - calculating the cutoffs
     - and caclulating the probabilities
    The operation is finalized by writing out a csv with the name pattern intramat<_name>_{mol_i}_{mol_j}.ndx
    
    Parameters
    ----------
    args : dict
        The command-line parsed parameters
    '''
    topology_mego, topology_ref, N_molecules, molecules_name, mol_list = read_topologies(args.mego_top, args.target_top)

    print(f"""
    Topology contains {N_molecules} molecules species. Namely {molecules_name}. 
    Calculating intramat for all species
    """)
    for i in range(N_molecules):
        print(f"\n Calculating intramat for molecule {mol_list[i]}: {molecules_name[i]}")
        df = pd.DataFrame()
        topology_df = pd.DataFrame()

        prefix=f"intra_mol_{mol_list[i]}_{mol_list[i]}"
        target_list = [ x for x in os.listdir(args.histo) if prefix in x ]

        protein_mego = topology_mego.molecules[list(topology_mego.molecules.keys())[i]][0]
        protein_ref = topology_ref.molecules[list(topology_ref.molecules.keys())[i]][0]
        original_size = len(protein_ref.atoms)
        protein_ref_indices = np.array([ i+1 for i in range(len(protein_ref.atoms)) if protein_ref[i].element_name != 'H' ])
        protein_ref = [ a for a in protein_ref.atoms if a.element_name != 'H' ]

        sorter = [ str(x.residue.number) + map_if_exists(x.name) for x in protein_ref ]
        sorter_mego = [ str(x.residue.number) + x.name for x in protein_mego ]

        topology_df['ref_ai'] = protein_ref_indices
        topology_df['ref_type'] = [ a.name for a in protein_ref ]
        topology_df['sorter'] = sorter
        topology_df['ref_ri'] = topology_df['sorter'].str.replace('[a-zA-Z]+[0-9]*', '', regex=True).astype(int)
        topology_df.sort_values(by='sorter', inplace=True)
        topology_df['mego_type'] = [ a[0].type for a in sorted(zip(protein_mego, sorter_mego), key=lambda x: x[1]) ]
        topology_df['mego_name'] = [ a[0].name for a in sorted(zip(protein_mego, sorter_mego), key=lambda x: x[1]) ]
        # need to sort back otherwise c12_cutoff are all wrong
        topology_df.sort_values(by='ref_ai', inplace=True)
        topology_df['c12'] = topology_df['mego_type'].map(d)

        #define all cutoff
        c12_cutoff = CUTOFF_FACTOR * np.power(np.sqrt(topology_df['c12'].values * topology_df['c12'].values[:,np.newaxis]),1./12.)
        if np.any(c12_cutoff>args.cutoff): warning_cutoff_histo(args.cutoff, np.max(c12_cutoff) )
        #c12_cutoff = c12_cutoff*0+0.75

        ########################
        # PARALLEL PROCESS START
        ########################

        chunks = np.array_split(target_list, args.proc)
        pool = multiprocessing.Pool(args.proc)
        results = pool.map(run_intra_, [ (args, protein_ref_indices, protein_ref_indices, original_size, c12_cutoff, mol_list[i], mol_list[i], x) for x in chunks ])
        pool.close()
        pool.join()

        ########################
        # PARALLEL PROCESS END
        ########################

        # concatenate and remove partial dataframes
        for name in results:
            part_df = pd.read_csv(name)
            df = pd.concat([df, part_df])
        [ os.remove(name) for name in results ]

        df = df.astype({
             'mi': 'int32', 
             'mj': 'int32', 
             'ai': 'int32', 
             'aj': 'int32', 
            })

        df = df.sort_values(by = ['mi', 'mj', 'ai', 'aj'])
        df.drop_duplicates(subset=['mi', 'ai', 'mj', 'aj'], inplace=True)

        df['mi'] = df['mi'].map('{:}'.format)
        df['mj'] = df['mj'].map('{:}'.format)
        df['ai'] = df['ai'].map('{:}'.format)
        df['aj'] = df['aj'].map('{:}'.format)
        df['dist'] = df['dist'].map('{:,.6f}'.format)
        df['c12dist'] = df['c12dist'].map('{:,.6f}'.format)
        df['hdist'] = df['hdist'].map('{:,.6f}'.format)
        df['p'] = df['p'].map('{:,.6e}'.format)
        df['cutoff'] = df['cutoff'].map('{:,.6f}'.format)

        df.index = range(len(df.index))
        out_name = args.out_name+'_' if args.out_name else ''
        
        output_file=f'{args.out}/intramat_{out_name}{mol_list[i]}_{mol_list[i]}.ndx'
        print(f"Saving output for molecule {mol_list[i]} in {output_file}")
        
        df.to_csv(output_file, index=False, sep=' ', header=False)

def calculate_inter_probabilities(args):
    '''
    Starts the main routine for calculating the intermat by:
     - reading the topologies
     - figuring out all the interacting molecules
     - calculating the cutoffs
     - and caclulating the probabilities
    The operation is finalized by writing out a csv with the name pattern intermat<_name>_{mol_i}_{mol_j}.ndx
    
    Parameters
    ----------
    args : dict
        The command-line parsed parameters
    '''
    topology_mego, topology_ref, N_species, molecules_name, mol_list = read_topologies(args.mego_top, args.target_top)
    pairs=list(itertools.combinations_with_replacement(mol_list,2))
    
    chain_list=[]
    chains = [x for x in topology_mego.molecules]

    for i in chains: 
        chain_list.append((i, len(topology_mego.molecules[i][0].atoms), len(topology_mego.split()[list(topology_mego.molecules.keys()).index(i)][1])))

    #number of molecules per species
    N_mols=[]
    for chain in chain_list:
        N_mols.append(chain[2])
    N_mols=np.array(N_mols)

    print(f"""
    Topology contains {N_species} molecules species. Namely {molecules_name}. 
    Calculating intermat for all species\n\n
    """)
    for pair in pairs:

        df = pd.DataFrame()

        topology_df_i = pd.DataFrame()
        topology_df_j = pd.DataFrame()
 
        mol_i=pair[0]
        mol_j=pair[1]        
        
        print(f"\nCalculating intermat between molecule {mol_i} and {mol_j}: {molecules_name[mol_i-1]} and {molecules_name[mol_j-1]}")
        prefix = f"inter_mol_{mol_i}_{mol_j}"
        # prefix_cum = f'inter_mol_c_{mol_i}_{mol_j}'
        target_list = [ x for x in os.listdir(args.histo) if prefix in x ]
        # target_list_cum = [ x for x in os.listdir(args.histo) if prefix_cum in x ]
        # target_list_norm = sorted(target_list_norm)
        # target_list_cum = sorted(target_list_cum)
        # target_list = list(zip(target_list_norm, target_list_cum))
        # for n, c in target_list:
        #     n = n.replace('inter_mol_','')
        #     c = c.replace('inter_mol_c_','')
        #     assert n == c, f'inter_mol {n} and inter_mol_d {c} are not the same'

        protein_mego_i = topology_mego.molecules[list(topology_mego.molecules.keys())[mol_i-1]][0]
        protein_mego_j = topology_mego.molecules[list(topology_mego.molecules.keys())[mol_j-1]][0]

        protein_ref_i = topology_ref.molecules[list(topology_ref.molecules.keys())[mol_i-1]][0]
        protein_ref_j = topology_ref.molecules[list(topology_ref.molecules.keys())[mol_j-1]][0]

        original_size_i = len(protein_ref_i.atoms)
        original_size_j = len(protein_ref_j.atoms)

        if mol_i==mol_j:
            if N_mols[mol_i-1]==0:
                print(f"Skipping intermolecular calculation between {mol_i} and {mol_j} cause the number of molecules of this species is only {N_mols[mol_i-1]}")
                columns=['mi', 'ai', 'mj', 'aj', 'dist', 'c12dist' , 'hdist' , 'p' , 'cutoff']
                matrix_index = pd.MultiIndex.from_product([ range(1,original_size_i+1) , range(1, original_size_j+1)], names=['ai', 'aj'])
                indeces_ai=np.array(list(matrix_index)).T[0]
                indeces_aj=np.array(list(matrix_index)).T[1]
                df=pd.DataFrame(columns=columns)
                df['mi'] = [ mol_i for x in range(1, original_size_i*original_size_j+1) ]
                df['mj'] = [ mol_j for x in range(1, original_size_i*original_size_j+1) ]
                df['ai'] = indeces_ai
                df['aj'] = indeces_aj
                df['dist']    = 0.
                df['c12dist'] = 0.
                df['hdist']   = 0.
                df['p']       = 0.
                df['cutoff']  = 0.
                df['mi'] = df['mi'].map('{:}'.format)
                df['mj'] = df['mj'].map('{:}'.format)
                df['ai'] = df['ai'].map('{:}'.format)
                df['aj'] = df['aj'].map('{:}'.format)
                df['dist'] = df['dist'].map('{:,.6f}'.format)
                df['c12dist'] = df['c12dist'].map('{:,.6f}'.format)
                df['hdist'] = df['hdist'].map('{:,.6f}'.format)
                df['p'] = df['p'].map('{:,.6e}'.format)
                df['cutoff'] = df['cutoff'].map('{:,.6f}'.format)

                df.index = range(len(df.index))
                out_name = args.out_name+'_' if args.out_name else ''
                output_file=f'{args.out}/intermat_{out_name}{mol_i}_{mol_j}.ndx'
        
                df.to_csv(output_file, index=False, sep=' ', header=False)
                continue

        protein_ref_indices_i = np.array([ i+1 for i in range(len(protein_ref_i.atoms)) if protein_ref_i[i].element_name != 'H' ])
        protein_ref_indices_j = np.array([ i+1 for i in range(len(protein_ref_j.atoms)) if protein_ref_j[i].element_name != 'H' ])

        protein_ref_i = [ a for a in protein_ref_i.atoms if a.element_name != 'H' ]
        protein_ref_j = [ a for a in protein_ref_j.atoms if a.element_name != 'H' ]

        sorter_i = [ str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_i ]
        sorter_mego_i = [ str(x.residue.number) + x.name for x in protein_mego_i ]

        sorter_j = [ str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_j ]
        sorter_mego_j = [ str(x.residue.number) + x.name for x in protein_mego_j ]

        #preparing topology of molecule i
        topology_df_i['ref_ai'] = protein_ref_indices_i
        topology_df_i['ref_type'] = [ a.name for a in protein_ref_i ]
        topology_df_i['sorter'] = sorter_i
        topology_df_i['ref_ri'] = topology_df_i['sorter'].str.replace('[a-zA-Z]+[0-9]*', '', regex=True).astype(int)
        topology_df_i.sort_values(by='sorter', inplace=True)
        topology_df_i['mego_type'] = [ a[0].type for a in sorted(zip(protein_mego_i, sorter_mego_i), key=lambda x: x[1]) ]
        topology_df_i['mego_name'] = [ a[0].name for a in sorted(zip(protein_mego_i, sorter_mego_i), key=lambda x: x[1]) ]
        # need to sort back otherwise c12_cutoff are all wrong
        topology_df_i.sort_values(by='ref_ai', inplace=True)
        topology_df_i['c12'] = topology_df_i['mego_type'].map(d)

        #preparing topology of molecule j
        topology_df_j['ref_ai'] = protein_ref_indices_j
        topology_df_j['ref_type'] = [ a.name for a in protein_ref_j ]
        topology_df_j['sorter'] = sorter_j
        topology_df_j['ref_ri'] = topology_df_j['sorter'].str.replace('[a-zA-Z]+[0-9]*', '', regex=True).astype(int)
        topology_df_j.sort_values(by='sorter', inplace=True)
        topology_df_j['mego_type'] = [ a[0].type for a in sorted(zip(protein_mego_j, sorter_mego_j), key=lambda x: x[1]) ]
        topology_df_j['mego_name'] = [ a[0].name for a in sorted(zip(protein_mego_j, sorter_mego_j), key=lambda x: x[1]) ]
        # need to sort back otherwise c12_cutoff are all wrong
        topology_df_j.sort_values(by='ref_ai', inplace=True)
        topology_df_j['c12'] = topology_df_j['mego_type'].map(d)

        #define all cutoff
        c12_cutoff = CUTOFF_FACTOR * np.power(np.sqrt(topology_df_j['c12'].values * topology_df_i['c12'].values[:,np.newaxis]),1./12.)
        if np.any(c12_cutoff>args.cutoff): warning_cutoff_histo(args.cutoff, np.max(c12_cutoff) )
        #c12_cutoff = c12_cutoff*0 + 0.75

        ########################
        # PARALLEL PROCESS START
        ########################

        chunks = np.array_split(target_list, args.proc)
        pool = multiprocessing.Pool(args.proc)
        results = pool.map(run_inter_, [ (args, protein_ref_indices_i, protein_ref_indices_j, original_size_j, c12_cutoff, mol_i, mol_j, x) for x in chunks ])
        pool.close()
        pool.join()

        ########################
        # PARALLEL PROCESS END
        ########################

        # concatenate and remove partial dataframes
        for i, name in enumerate(results):
            part_df = pd.read_csv(name)
            df = pd.concat([df, part_df])
            os.remove(name)

        df = df.astype({
             'mi': 'int32', 
             'mj': 'int32', 
             'ai': 'int32', 
             'aj': 'int32'
            })

        df = df.sort_values(by = ['mi', 'mj', 'ai', 'aj'])
        df.drop_duplicates(subset=['mi', 'ai', 'mj', 'aj'], inplace=True)
        
        df['mi'] = df['mi'].map('{:}'.format)
        df['mj'] = df['mj'].map('{:}'.format)
        df['ai'] = df['ai'].map('{:}'.format)
        df['aj'] = df['aj'].map('{:}'.format)
        df['dist'] = df['dist'].map('{:,.6f}'.format)
        df['c12dist'] = df['c12dist'].map('{:,.6f}'.format)
        df['hdist'] = df['hdist'].map('{:,.6f}'.format)
        df['p'] = df['p'].map('{:,.6e}'.format)
        df['cutoff'] = df['cutoff'].map('{:,.6f}'.format)

        df.index = range(len(df.index))
        out_name = args.out_name+'_' if args.out_name else ''
        output_file=f'{args.out}/intermat_{out_name}{mol_i}_{mol_j}.ndx'
        
        df.to_csv(output_file, index=False, sep=' ', header=False)

def calculate_probability(values, weights, callback=allfunction):
    '''
    Calculates a plain probability accoring to \sum_x x * dx

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights
    callback : function
        Preprocesses the data before going in to the analysis

    Returns
    -------
    The probability of the histogram
    '''
    dx = values[1] - values[0]
    cutoff, i, norm, v, w = callback(values, weights)
    return np.minimum( np.sum(w * dx), 1 )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--histo', type=str, required=True, help='Path to the directory containing the histograms. The histogram files should contain the prefix "intra_" for intra molecular contact descriptions and "inter_" for  inter molecular.')
    parser.add_argument('--target_top', required=True, help='Path to the topology file of the system on which the histograms were calculated on')
    parser.add_argument('--mego_top', required=True, help='''Path to the standard multi-eGO topology of the system generated by pdb2gmx''')
    parser.add_argument('--inter', action='store_true', help='Sets the caculation to be adapted to intermolecular calculations of histograms')
    parser.add_argument('--out', default='./', help='''Sets the output path''')
    parser.add_argument('--out_name', help='''Sets the output name of files to be added to the default one: intermat_<out_name>_mi_mj.ndx or intramat_<out_name>_mi_mj.ndx''')
    parser.add_argument('--proc', default=1, type=int, help='Sets the number of processes to perform the calculation')
    parser.add_argument('--cutoff', required=True, type=float, help='To be set to the max cutoff used for the accumulation of the histograms')
    args = parser.parse_args()
    
    #check if output file exists
    if not os.path.exists(args.out):
        print(f"The path '{args.out}' does not exist.")
        sys.exit()

    N_BINS = args.cutoff / ( 0.01 / 4 )
    DX = args.cutoff / N_BINS
    PREFIX = 'inter_mol_' if args.inter else 'intra_mol_'
    CUTOFF_FACTOR=1.45
    print(f"""
    Starting with cutoff = {args.cutoff},
                  n_bins = {N_BINS},
                  dx     = {DX}
                  on {args.proc} threads
    """)

    if not args.inter:
        calculate_intra_probabilities(args)
    else:
        calculate_inter_probabilities(args)
