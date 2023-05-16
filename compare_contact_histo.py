import os
import sys
import pandas as pd
import numpy as np
import parmed as pmd
import multiprocessing

################################################################
import warnings                                                #
warnings.filterwarnings("ignore")                              #
################################################################

################################################################
# TODO                  maybe                                  #
################################################################
# - set datatypes to double or rather float64                  #
################################################################

REF_PATH = sys.argv[1]      # histograms to analyze
REF_TOP = sys.argv[2]       # 
MEGO_TOP = sys.argv[3]
OUT = sys.argv[4]
CUTOFF = 0.75
N_BINS = CUTOFF / ( 0.01 / 4 )
DX = CUTOFF / N_BINS
N_THREADS = 1
EPSILON = 0.5

print(f"""
Starting with cutoff = {CUTOFF},
              n_bins = {N_BINS},
              dx     = {DX}
              on {N_THREADS} threads
""")

from_ff_to_multiego = {
    'OC1' : 'O1',
    'OC2' : 'O2',
    'OT1' : 'O1',
    'OT2' : 'O2',
    'C13' :'CN1',
    'C14' :'CN2',
    'C15' :'CN3',
    'N'   :'N',
    'C12' :'CA',
    'C11' :'CB',
    'O12' :'OA',
    'P'   :'P',
    'O13' :'OB',
    'O14' :'OC',
    'O11' :'OD',
    'C1'  :'CC',
    'C2'  :'CD',
    'O21' :'OE',
    'C21' :'C1A',
    'O22' :'OF',
    'C22' :'C1B',
    'C23' :'C1C',
    'C24' :'C1D',
    'C25' :'C1E',
    'C26' :'C1F',
    'C27' :'C1G',
    'C28' :'C1H',
    'C29' :'C1I',
    'C210':'C1J',
    'C211':'C1K',
    'C212':'C1L',
    'C213':'C1M',
    'C214':'C1N',
    'C215':'C1O',
    'C216':'C1P',
    'C217':'C1Q',
    'C218':'C1R',
    'C3'  :'CE',
    'O31' :'OG',
    'C31' :'C2A',
    'O32' :'OH',
    'C32' :'C2B',
    'C33' :'C2C',
    'C34' :'C2D',
    'C35' :'C2E',
    'C36' :'C2F',
    'C37' :'C2G',
    'C38' :'C2H',
    'C39' :'C2I',
    'C310':'C2J',
    'C311':'C2K',
    'C312':'C2L',
    'C313':'C2M',
    'C314':'C2N',
    'C315':'C2O',
    'C316':'C2P',
}

gromos_atp = pd.DataFrame(
    {'name': ['O', 'OA', 'N', 'C', 'CH1', 
            'CH2', 'CH3', 'CH2r', 'NT', 'S',
            'NR', 'OM', 'NE', 'NL', 'NZ',
            'CH3p', 'P', 'OE', 'CR1'],            
     'c12': [2.631580e-07, 5.018430e-07, 8.752940e-07, 2.598570e-06, 6.555574e-05, # CH1
             1.543890e-05, 8.595562e-06, 1.193966e-05, 2.596154e-06, 2.724050e-06, 
             1.506347e-06, 1.724403e-07, 8.752940e-07, 8.752940e-07, 8.752940e-07,
             8.736473e-06, 3.893600e-06, 3.558824e-07, 6.29856e-06]
     }
)

d = { gromos_atp.name[i] : gromos_atp.c12[i] for i in range(len(gromos_atp.name))}

def run_(frac_target_list):
    process = multiprocessing.current_process()
    columns = ['mi', 'ai', 'mj', 'aj', 'dist', 'c12dist', 'p', 'sigma', 'cutoff', 'is_gauss']

    df = pd.DataFrame(columns=columns)
    for i, ref_f in enumerate(frac_target_list):
        results_df = pd.DataFrame()
        ai = ref_f.split('.')[-2].split('_')[-1]
        if not np.isin(ai, protein_ref_indices):
            all_ai = [ ai for _ in range(1, original_size+1) ]
            range_list = [ str(x) for x in range(1, original_size+1) ]
            results_df['ai'] = [ *all_ai, *range_list ]
            results_df['aj'] = [ *range_list, *all_ai ]
            results_df['mi'] = 1
            results_df['mj'] = 1
            results_df['sigma'] = -1
            results_df['dist'] = 0
            results_df['c12dist'] = 0
            results_df['p'] = 0
            results_df['cutoff'] = 0
            results_df['is_gauss'] = 0
        else:
            cut_i = np.where(protein_ref_indices == int(ai))[0][0]
         
            # column mapping
            ref_f = f'{REF_PATH}/{ref_f}'
            ref_df = pd.read_csv(ref_f, header=None, sep='\s+', usecols=[ 0, *protein_ref_indices ])
            ref_df_columns = ['distance', *[ str(x) for x in protein_ref_indices ]]
            ref_df.columns = ref_df_columns
            ref_df.set_index('distance', inplace=True)
            ref_df.loc[len(ref_df)] = c12_cutoff[cut_i]

            # calculate data
            dist = ref_df.apply(lambda x: weighted_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            c12dist = ref_df.apply(lambda x: c12_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            sigma = ref_df.apply(lambda x: calc_sigma(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            p = ref_df.apply(lambda x: calculate_probability(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values

            results_df['mi'] = [ 1 for x in range(1, len(ref_df.columns)+1) ]
            results_df['mj'] = [ 1 for x in range(1, len(ref_df.columns)+1) ]
            results_df['ai'] = [ ai for _ in range(1, len(ref_df.columns)+1) ]
            results_df['aj'] = [ str(x) for x in protein_ref_indices ]
            results_df['sigma'] = sigma 
            results_df['dist'] = dist
            results_df['c12dist'] = c12dist
            results_df['p'] = p
            results_df['cutoff'] = c12_cutoff[cut_i]
            results_df['is_gauss'] = ref_df.apply(lambda x: single_gaussian_check(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
    
        df = pd.concat([df, results_df])
        df = df.sort_values(by = ['p', 'c12dist', 'dist'], ascending=True)

    df.fillna(0)
    out_path = f'mat_{process.pid}.part'
    df.to_csv(out_path, index=False)

    return out_path

def map_if_exists(x):
    if x in from_ff_to_multiego.keys(): return from_ff_to_multiego[x]
    else: return x

def allfunction(values, weights):
    v = values[:-1]
    cutoff = weights[len(weights)-1]
    w = weights[:-1]
    i = np.where(v <= cutoff)
    if not np.any(i): return 0,0,0,0,0  # check if empty
    i = i[0]
    w = w[i]
    v = v[i]
    norm = np.sum(w)
    return cutoff, i, norm, v, w

def zero_callback(values, weights):
    return None, None, np.sum(weights), values, weights

def weighted_avg(values, weights, callback=allfunction):
    cutoff, i, norm, v, w = callback(values, weights)
    if norm == 0.: return 0
    return np.sum(v * w) / norm

def single_gaussian_check(values, weights, callback=allfunction):
    cutoff, i, norm, values, weights = callback(values, weights)
    values, weights = remove_monotonic(values, weights)
    a = weights[:-2]
    b = weights[2:]
    slope = (b - a) / (2. * DX)
    danger_sign=0
    danger_trend=0
    increasing=1
    for i in range(slope.size-1):
        if slope[i+1] != 0. or slope[i] != 0.:
            if increasing and slope[i] > 15. and slope[i+1] < slope[i]:
                increasing=0
                danger_trend+=1
            if not increasing and slope[i+1] > slope[i]:
                increasing=1
                danger_trend+=1
            if danger_trend > 2: return 0
            if not danger_sign and danger_trend > 1: return 0
            if slope[i] * slope[i+1] < 0.:
                if danger_sign: return 0
                danger_sign = 1

    return 1

def c12_avg(values, weights, callback=allfunction):
    single_gaussian = single_gaussian_check(values, weights)
    cutoff, i, norm, v, w = callback(values, weights)
    if norm == 0.: return 0

    # calculate the sigma
    sigma = calc_sigma(v, w, zero_callback)
    
    r = np.where(w > 0.)
    if r[0].size == 0: return 0
    maxd = v[r[0][-1]]
    mind = v[r[0][0]]
    sigma_cut = (maxd-mind)/4.
    i_start = r[0][0]
    
    # apply sigma rules
    if sigma > sigma_cut or not single_gaussian:
        i_stop = int(w.size - (w.size - i_start) / 2)
        w = w[i_start:i_stop] 
        v = v[i_start:i_stop]

    norm = np.sum(w)    
    # if single_gaussian: return np.sum(weights * values) / norm
    return np.power( 1. / ( np.sum(w*np.power(1./v, 12.)) / norm ), 1. / 12.)

def calc_sigma(values, weights, callback=allfunction):
    cutoff, i, norm, v, w = callback(values, weights)
    wo = w
    v, w = remove_monotonic(v, w)
    norm = np.sum(w)
    if norm == 0.: return -1.
    avg = np.sum(w * v)
    avg2 = np.sum(w * v * v)
    return np.sqrt( avg2 / norm - avg * avg / norm / norm )

def calculate_probability(values, weights, callback=allfunction):
    cutoff, i, norm, v, w = callback(values, weights)
    return np.minimum( np.sum(w * DX), 1 )

df = pd.DataFrame()
target_list = [ x for x in os.listdir(REF_PATH) if 'intra_mol_' in x ]

topology_df = pd.DataFrame()
topology_mego = pmd.load_file(MEGO_TOP)
topology_ref = pmd.load_file(REF_TOP)
protein_mego = topology_mego.molecules[list(topology_mego.molecules.keys())[0]][0]
protein_ref = topology_ref.molecules[list(topology_ref.molecules.keys())[0]][0]
original_size = len(protein_ref.atoms)
protein_ref_indices = np.array([ i+1 for i in range(len(protein_ref.atoms)) if protein_ref[i].element_name != 'H' ])
protein_ref = [ a for a in protein_ref.atoms if a.element_name != 'H' ]

sorter = [ str(x.residue.number) + map_if_exists( x.name) for x in protein_ref ]
sorter_mego = [ str(x.residue.number) + x.name for x in protein_mego ]

topology_df['ref_ai'] = protein_ref_indices
topology_df['ref_type'] = [ a.name for a in protein_ref ]
topology_df['sorter'] = sorter
topology_df.sort_values(by='sorter', inplace=True)
topology_df['mego_type'] = [ a[0].type for a in sorted(zip(protein_mego, sorter_mego), key=lambda x: x[1]) ]
#need to sort back otherwise c12_cutoff are all wrong
topology_df.sort_values(by='ref_ai', inplace=True)
topology_df['c12'] = topology_df['mego_type'].map(d)

c12_cutoff = 1.45 * np.power(np.sqrt(topology_df['c12'].values * topology_df['c12'].values[:,np.newaxis]),1./12.)
#c12_cutoff = 0. * np.power(np.sqrt(topology_df['c12'].values * topology_df['c12'].values[:,np.newaxis]),1./12.) + 0.55

if __name__ == '__main__':
    ########################
    # PARALLEL PROCESS START
    ########################

    chunks = np.array_split(target_list, N_THREADS)
    pool = multiprocessing.Pool(N_THREADS)
    results = pool.map(run_, chunks)
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

    print(df)

    df = df.astype({
         'mi': 'int32', 
         'mj': 'int32', 
         'ai': 'int32', 
         'aj': 'int32', 
         'is_gauss': 'int32'
        })

    df = df.sort_values(by = ['mi', 'mj', 'ai', 'aj'])

    df['mi'] = df['mi'].map('{:}'.format)
    df['mj'] = df['mj'].map('{:}'.format)
    df['ai'] = df['ai'].map('{:}'.format)
    df['aj'] = df['aj'].map('{:}'.format)
    df['is_gauss'] = df['is_gauss'].map('{:}'.format)
    df['dist'] = df['dist'].map('{:,.6f}'.format)
    df['sigma'] = df['sigma'].map('{:,.6f}'.format)
    df['c12dist'] = df['c12dist'].map('{:,.6f}'.format)
    df['p'] = df['p'].map('{:,.6f}'.format)

    df.drop_duplicates(subset=['mi', 'ai', 'mj', 'aj'], inplace=True)

    df.index = range(len(df.index))
    df.to_csv(OUT, index=False, sep=' ', header=False)
