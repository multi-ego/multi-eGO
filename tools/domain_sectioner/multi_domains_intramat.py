import numpy as np 
import sys
import itertools

intramat = sys.argv[1]
blocks = [ int(x) for x in sys.argv[2:] ]

if len(sys.argv)<3:
    print("Error")
    exit()

# block_pairs = itertools.product(blocks, repeat=2)
index_list = []
intra_md = np.loadtxt(intramat, unpack=True)
dim = int(np.sqrt(len(intra_md[0])))
domain_mask = np.full((dim, dim), False)
full_blocks = [0, *blocks, dim]
for i, _ in enumerate([0, *blocks]):
    print(f"Dividing {intramat} at atoms {full_blocks[i]} and {full_blocks[i+1]}")
    map = np.array([ True if x > full_blocks[i] and x <= full_blocks[i+1] else False for x in range(dim)])
    # map_j = np.array([ True if i >  i_block - 1 else False for i in range(dim)])
    map = map * map[:,np.newaxis]
    # map = map.reshape(dim**2)

    domain_mask |= map

domain_mask_linear = domain_mask.reshape(dim**2)
intra_md[4] = np.where(domain_mask_linear, intra_md[4], 0.)
intra_md[5] = np.where(domain_mask_linear, intra_md[5], 0.)

prefix = './'
if '/' in intramat:
    intramat = intramat.split('/')[-1]
    # prefix = "/".join(intramat.split('/')[:-1])

print(f'split_{"-".join(np.array(blocks, dtype=str))}_{intramat}')
np.savetxt(f'split_{"-".join(np.array(blocks, dtype=str))}_{intramat}',intra_md.T, delimiter=" ", fmt = ['%i', '%i', '%i', '%i', '%2.6f', '%.6e', '%2.6f' ])
print(f"Finished job")
