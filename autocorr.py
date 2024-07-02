import numpy as np
from numba import njit
import h5py
from metdense2 import MetdenseDataset, MethCall

@njit
def process_chrom( chrom, num_dist_bins=17 ):
    corr_counts = np.zeros( ( num_dist_bins, chrom.ncells, 2, 2 ), dtype=np.int64 )

    for pos_idx_1 in range(chrom.posv.shape[0]):
        pos_idx_2 = pos_idx_1
        distbin_idx = 0
        cur_max_dist = 1
        while True:
            pos_idx_2 += 1
            dist = chrom.posv[pos_idx_2] - chrom.posv[pos_idx_1]

            # Find correct distbin
            while True:
                if dist <= cur_max_dist:
                    break
                cur_max_dist *= 2
                distbin_idx += 1
                if distbin_idx >= num_dist_bins:
                    break
            if distbin_idx >= num_dist_bins:
                break

            for cell_idx in range(chrom.ncells):
                call1 = chrom.get(pos_idx_1, cell_idx)
                if call1 == MethCall.nocall or call1 == MethCall.ambig:
                    continue
                call2 = chrom.get(pos_idx_2, cell_idx)
                if call2 == MethCall.nocall or call2 == MethCall.ambig:
                    continue
                corr_counts[ distbin_idx, cell_idx, call1-1, call2-1 ] += 1

        if pos_idx_1 % 1000 == 0:
            print( pos_idx_1, "/", chrom.posv.shape[0] )

    return corr_counts

md = MetdenseDataset( "scnmt_data__CpG__filt.metdense" )

corr_counts = process_chrom( md.chroms["10"] )

print( corr_counts.sum(1) )

with h5py.File('corr_counts_chr10.h5', 'w') as f:
    f.create_dataset( 'corr_counts', data=corr_counts )