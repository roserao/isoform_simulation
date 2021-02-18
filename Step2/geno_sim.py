import numpy as np
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='sample random genotype data from LD matrix')
    parser.add_argument('--gene', type=str,
                        help='generate extra files for GREML analysis')
    parser.add_argument('--size', type=int, default=0,
                        help="number of samples generated")
    parser.add_argument("--dir", type=str,
                        help="input and output directory")
    args = parser.parse_args()

    gene = args.gene
    nsample = args.size
    np.random.seed(105104394)

    ld = np.loadtxt(args.dir + "/" + gene + "_AFR.clean.ld", delimiter=' ', dtype = str)

    print(np.size(ld, 0))
    print(np.size(ld, 1))

    #ld = ld[:, :-1]
    ld = ld.astype(np.float)
    nsnp = np.size(ld, 0)
    mean = np.zeros(nsnp)

    geno_sample = np.random.multivariate_normal(mean, ld, size = nsample)
    np.savetxt(args.dir + "/" + gene + '_AFR_sample.txt', geno_sample, delimiter=',')

