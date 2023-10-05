import numpy as np 
import sys
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TODO!')
    parser.add_argument('--type', choices=['split', 'group'], required=True, help='Type of EGO.\n rc: creates a force-field for random coil simulations.\n production: creates a force-field combining random coil simulations and training simulations.')
    parser.add_argument('--input', type=str,required=True, help='intramat to work on')
    parser.add_argument('--input2', type=str, help='second intramat')

    parser.add_argument('--out', type=str, default="./", help='path for ouput')
    parser.add_argument('--iblocks', nargs='+', type=int, default=[], help='list of atom indeces associated to the the ending point of domains (excluded last) 15,44,..')

    args = parser.parse_args()

    # checking the options provided in the commandline
    if args.type != 'split' and args.type is None:
        print('--type=choose either split or group. ')
        sys.exit()

    if args.type != 'group' and args.type is None:
        print('--type=choose either split or group. ')
        sys.exit()

    if args.type == 'group' and args.input2 is None:
        print('--type=group requires 2 inputs: --input PATH_TO_intramat_rc --input2 PATH_TO_intramat_inter_domain_rc')
        sys.exit()

    if args.out:
        if not os.path.isdir(args.out):
            print(f"{args.out} does not exists. Insert an existing directory")
            exit()
        else:
            if args.out[-1]!="/":
                args.out=args.out + "/"


if args.type=="split":

    print("Generating intramat for inter-domain rancomd coil")
    print("")

    intramat = args.input
    
    #read intramat
    intra_md = np.loadtxt(intramat, unpack=True)

    dim = int(np.sqrt(len(intra_md[0])))
    domain_mask = np.full((dim, dim), False)
    full_blocks = [0, *args.iblocks, dim+1]

    for i, _ in enumerate([0, *args.iblocks]):

        print(f"Dividing {intramat} at atoms {full_blocks[i]} and {full_blocks[i+1]}")
        map = np.array([ True if x >= full_blocks[i] and x < full_blocks[i+1] else False for x in range(dim)])
        map = map * map[:,np.newaxis]
        domain_mask = np.logical_or(domain_mask, map)

    domain_mask_linear = domain_mask.reshape(dim**2)

    intra_md[4] = np.where(domain_mask_linear, intra_md[4], 0.)
    intra_md[5] = np.where(domain_mask_linear, intra_md[5], 0.)

    if '/' in intramat:
        intramat = intramat.split('/')[-1]

    np.savetxt(f'{args.out}split_{"-".join(np.array(args.iblocks, dtype=str))}_{intramat}',intra_md.T, delimiter=" ", fmt = ['%i', '%i', '%i', '%i', '%2.6f', '%.6e', '%2.6f' ])
    print(f"Finished job")


if args.type=="group":
    
    print("Group intramat_rc with intramat inter-domain_rc")
    print("")

    intra1 = args.input
    intra2 = args.input2

    #read intramats
    intra_rc = np.loadtxt(intra1, unpack=True)
    intra_domain_rc = np.loadtxt(intra2, unpack=True)

    if intra_rc.shape!=intra_domain_rc.shape:
        print("intramats of input 1 and 2 must have the same dimensions (they should be of the same system!)")
        exit()

    dim = int(np.sqrt(len(intra_rc[0])))
    domain_mask = np.full((dim, dim), False)
    full_blocks = [0, *args.iblocks, dim+1]

    for i, _ in enumerate([0, *args.iblocks]):

        print(f"Group {intra1} and {intra2} at atoms {full_blocks[i]} and {full_blocks[i+1]}")
        map = np.array([ True if x >= full_blocks[i] and x < full_blocks[i+1] else False for x in range(dim)])
        map = map * map[:,np.newaxis]
        domain_mask = np.logical_or(domain_mask, map)

    domain_mask_linear = domain_mask.reshape(dim**2)
    intra_rc[4] = np.where(domain_mask_linear, intra_rc[4], intra_domain_rc[4])
    intra_rc[5] = np.where(domain_mask_linear, intra_rc[5], intra_domain_rc[5])
    intra_rc[6] = np.where(domain_mask_linear, intra_rc[6], intra_domain_rc[6])

    if '/' in intra1:
        intra1 = intra1.split('/')[-1]

    np.savetxt(f'{args.out}group_{"-".join(np.array(args.iblocks, dtype=str))}_{intra1}',intra_rc.T, delimiter=" ", fmt = ['%i', '%i', '%i', '%i', '%2.6f', '%.6e', '%2.6f' ])
    print(f"Finished job")

    