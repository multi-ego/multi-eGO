import numpy as np 
import sys

i_block = int(sys.argv[1])

if len(sys.argv)<3:
    print("Error")
    exit()

elif len(sys.argv)==3:
    intra_md = np.loadtxt(sys.argv[2], unpack=True)
    print(f"Dividing {intra_md} at atom {i_block}")
    dim = int(np.sqrt(len(intra_md[0])))

    map_i = np.array([ True if i <= i_block - 1 else False for i in range(dim)])
    map_j = np.array([ True if i >  i_block - 1 else False for i in range(dim)])
    map_i = map_i * map_i[:,np.newaxis]
    map_j = map_j * map_j[:,np.newaxis]
    map = np.logical_or(map_j,map_i)
    map = map.reshape(dim**2)

    intra_md[4] = np.where(map, intra_md[4], 0.)
    intra_md[5] = np.where(map, intra_md[5], 0.)

    np.savetxt("split_intramat.ndx",intra_md.T, delimiter=" ", fmt = ['%i', '%i', '%i', '%i', '%2.6f', '%.6e', '%2.6f' ])
    print(f"Finished job ")

elif len(sys.argv)==4:
    intra_rc = np.loadtxt(sys.argv[2], unpack=True)
    intra_domain_rc = np.loadtxt(sys.argv[3], unpack=True)
    print("Group {intra_rc} with {intra_domain_rc} at atom {i_block}")
    if(intra_rc.shape != intra_domain_rc.shape): 
        print("input intramats have different dimensions")
        exit()
    dim = int(np.sqrt(len(intra_rc[0])))

    map_i = np.array([ True if i <= i_block - 1 else False for i in range(dim)])
    map_j = np.array([ True if i >  i_block - 1 else False for i in range(dim)])
    map_i = map_i * map_i[:,np.newaxis]
    map_j = map_j * map_j[:,np.newaxis]
    map = np.logical_or(map_j,map_i)
    map = map.reshape(dim**2)

    intra_rc[4] = np.where(map, intra_rc[4], intra_domain_rc[4])
    intra_rc[5] = np.where(map, intra_rc[5], intra_domain_rc[5])
    intra_rc[6] = np.where(map, intra_rc[6], intra_domain_rc[6])

    np.savetxt("grouped_intramat.ndx",intra_rc.T, delimiter=" ",fmt = ['%i', '%i', '%i', '%i', '%2.6f', '%.6e', '%2.6f' ])
    print("Finished job")

