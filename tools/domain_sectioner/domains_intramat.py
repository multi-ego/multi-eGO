import numpy as np 
import sys

i_block = int(sys.argv[1])

if len(sys.argv)<3:
    print("Error")
    exit()

elif len(sys.argv)==3:
    print("reading")
    intra_md = np.loadtxt(sys.argv[2], unpack=True)
    print("read")
    dim = int(np.sqrt(len(intra_md[0])))

    map_i = np.array([ True if i <= i_block - 1 else False for i in range(dim)])
    map_j = np.array([ True if i >  i_block - 1 else False for i in range(dim)])
    map_i = map_i * map_i[:,np.newaxis]
    map_j = map_j * map_j[:,np.newaxis]
    map = np.logical_or(map_j,map_i)
    map = map.reshape(dim**2)
    intra_md[4] = np.where(map, intra_md[4], 0)
    intra_md[5] = np.where(map, intra_md[5], 0)
    intra_md[6] = np.where(map, intra_md[6], 0)
    #intra_md[7] = np.where(map, intra_md[7], 0)
    #intra_md[8] = np.where(map, intra_md[8], 0)

    np.savetxt("split_intramat.ndx",intra_md.T, delimiter=" ", fmt = ['%4i', '%4i', '%4i', '%4i', '%2.6f', '%2.6f', '%2.6f' ])

elif len(sys.argv)==4:
    print("reading")
    intra_rc = np.loadtxt(sys.argv[2], unpack=True)
    intra_domain_rc = np.loadtxt(sys.argv[3], unpack=True)
    if(intra_rc.shape != intra_domain_rc.shape): 
        print("input intramats have different dimensions")
        exit()
    print("read")
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
    #intra_rc[7] = np.where(map, intra_rc[7], intra_domain_rc[7])
    #intra_rc[8] = np.where(map, intra_rc[8], intra_domain_rc[8])

    np.savetxt("grouped_intramat.ndx",intra_rc.T, delimiter=" ",fmt = ['%4i', '%4i', '%4i', '%4i', '%2.6f', '%2.6f', '%2.6' ])

