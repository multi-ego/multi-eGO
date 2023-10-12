import numpy as np 
import sys
import argparse
import os
import parmed as pmd
import warnings

def find_atom(top, res_num):
    atom_idx = 0
    for i in range(res_num-1):
        atom_idx += len(top.residues[i].atoms)

    return atom_idx


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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TODO!')
    parser.add_argument('--type', choices=['split', 'group'], required=True, help='Type of operation.\n split: splits the md intramat into blocks associated to the domains.\n group: group the rc and the inter domain rc into blocks associated to the domains.')
    parser.add_argument('--md_intra', type=str,required=False, help='intramat to work on')
    parser.add_argument('--rc_intra', type=str, help='random coil intramat')
    parser.add_argument('--dom_rc_intra', type=str, help='inter domain random coil intramat')
    parser.add_argument('--target_top', type=str)
    parser.add_argument('--mego_top', type=str)
    parser.add_argument('--dom_res', nargs='+', type=int, default=[], help='list of residue indeces associated to the starting residues of new domains. Example: 15,44,..')
    parser.add_argument('--out', type=str, default="./", help='path for ouput')

    args = parser.parse_args()

    # checking the options provided in the commandline
    if args.type != 'split' and args.type is None:
        print('--type=choose either split or group. ')
        sys.exit()

    if args.type == 'split' and (args.md_intra is None or args.target_top is None or args.mego_top is None):
        print('--type=split requires 3 inputs: --md_intra PATH_TO_intramat_md --target_top PATH_to_target_topology --mego_top PATH_to_mego_topology')
        sys.exit()

    if args.type != 'group' and args.type is None:
        print('--type=choose either split or group. ')
        sys.exit()

    if args.type == 'group' and (args.rc_intra is None or args.dom_rc_intra is None or args.target_top is None):
        print('--type=group requires 3 inputs: --rc_intra PATH_TO_intramat_rc --dom_rc_intra PATH_TO_intramat_inter_domain_rc --target_top PATH_to_target_topology')
        sys.exit()

    if args.type == 'group' and args.mego_top == None:
        args.mego_top = args.target_top

    if args.out:
        if not os.path.isdir(args.out):
            print(f"{args.out} does not exists. Insert an existing directory")
            exit()
        else:
            if args.out[-1]!="/":
                args.out=args.out + "/"

#read topology

topology_mego, topology_ref, N_species, molecules_name, mol_list = read_topologies(args.mego_top, args.target_top)

molecule_name = molecules_name[0]

if molecule_name not in list(topology_ref.molecules.keys()):
    print(f'''
        ERROR: 
        Molecule "{molecule_name}" found in mego topology is not found in {args.target_top}. 
        Molecules in {args.target_top} are: {list(topology_ref.molecules.keys())}
        Either you used the wrong topology or you need to use the same name in mego and ref topologies for the system''')
    exit()

top_prot = topology_ref.molecules[molecule_name][0]

if N_species > 1: 
    print('maximum 1 molecule specie')
    exit()

#define atom_num and res_num
n_atoms = len(top_prot.atoms)
n_res   = len(top_prot.residues)

print(f'\n Total number of residues {n_res} and total number of atoms {n_atoms} \n')

if args.type=="split":

    print("\n Generating intramat for inter-domain rancomd coil")
    print("")

    intramat = args.input
    
    #read intramat
    intra_md = np.loadtxt(intramat, unpack=True)
    dim = int(np.sqrt(len(intra_md[0])))

    #define domain mask
    domain_mask = np.full((dim, dim), False)

    #consistency check
    if(dim!=n_atoms):
        print(f'ERROR: number of atoms in intramat ({dim}) does not correspond to that of topology ({n_atoms})')
        exit()
    
    res_idx = args.dom_res
    #find atom indeces associated to starting domain residues
    atom_idxs = []
    for ri in res_idx:
        atom_idxs.append(find_atom(topology_ref, ri))

    full_blocks = [0, *atom_idxs, dim+1]
    full_blocks_res = [1, *res_idx, n_res]

    for i, _ in enumerate([0, *atom_idxs]):

        print(f"Dividing {intramat} at residues {full_blocks_res[i]} - {full_blocks_res[i+1]} and atoms {full_blocks[i] + 1} - {full_blocks[i+1]}")
        map = np.array([ True if x >= full_blocks[i] and x < full_blocks[i+1] else False for x in range(dim)])  
        map = map * map[:,np.newaxis]
        domain_mask = np.logical_or(domain_mask, map)

    domain_mask_linear = domain_mask.reshape(dim**2)

    intra_md[4] = np.where(domain_mask_linear, intra_md[4], 0.)
    intra_md[5] = np.where(domain_mask_linear, intra_md[5], 0.)

    if '/' in intramat:
        intramat = intramat.split('/')[-1]

    np.savetxt(f'{args.out}split_{"-".join(np.array(args.dom_res, dtype=str))}_{intramat}',intra_md.T, delimiter=" ", fmt = ['%i', '%i', '%i', '%i', '%2.6f', '%.6e', '%2.6f' ])
    print(f"Finished splitting")


if args.type=="group":
    
    print("Group intramat_rc with intramat inter-domain_rc")
    print("")

    intra1 = args.input
    intra2 = args.input2

    #read intramats
    intra_rc = np.loadtxt(intra1, unpack=True)
    intra_domain_rc = np.loadtxt(intra2, unpack=True)

    #first consistency check
    if intra_rc.shape!=intra_domain_rc.shape:
        print("intramats of input 1 and 2 must have the same dimensions (they should be of the same system!)")
        exit()

    dim = int(np.sqrt(len(intra_rc[0])))

    #define domain mask
    domain_mask = np.full((dim, dim), False)

    #second consistency check
    if(dim!=n_atoms):
        print(f'ERROR: number of atoms in intramat ({dim}) does not correspond to that of topology ({n_atoms})')
        exit()
        
    res_idx = args.dom_res
    
    #find atom indeces associated to starting domain residues
    atom_idxs = []
    for ri in res_idx:
        atom_idxs.append(find_atom(topology_ref, ri))

    full_blocks = [0, *atom_idxs, dim+1]
    full_blocks_res = [1, *res_idx, n_res]
    for i, _ in enumerate([0, *atom_idxs]):

        print(f"Group {intra1} and {intra2} at residues {full_blocks_res[i]} - {full_blocks_res[i+1]} and atoms {full_blocks[i] + 1} - {full_blocks[i+1]} ")

        map = np.array([ True if x >= full_blocks[i] and x < full_blocks[i+1] else False for x in range(dim)])
        map = map * map[:,np.newaxis]
        domain_mask = np.logical_or(domain_mask, map)

    domain_mask_linear = domain_mask.reshape(dim**2)
    intra_rc[4] = np.where(domain_mask_linear, intra_rc[4], intra_domain_rc[4])
    intra_rc[5] = np.where(domain_mask_linear, intra_rc[5], intra_domain_rc[5])
    intra_rc[6] = np.where(domain_mask_linear, intra_rc[6], intra_domain_rc[6])

    if '/' in intra1:
        intra1 = intra1.split('/')[-1]

    np.savetxt(f'{args.out}group_{"-".join(np.array(args.dom_res, dtype=str))}_{intra1}',intra_rc.T, delimiter=" ", fmt = ['%i', '%i', '%i', '%i', '%2.6f', '%.6e', '%2.6f' ])
    print(f"Finished group")

    