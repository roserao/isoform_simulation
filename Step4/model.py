import numpy as np
import argparse
import sys
import os
import pandas as pd
from sklearn import covariance
from pysnptools.snpreader import Bed

class BaseModel:

    def __init__(self, gene, iso, sim):
        seed = 124
        np.random.seed(seed)
        self.gene = gene
        self.num_iso = iso
        self.num_sim = sim
        bfile = gene + "_AFR.clean"
        geno = Bed(bfile, count_A1=False).read().val
        self.n_ind, self.n_SNP = geno.shape
        print(geno.shape)
        f = np.sum(geno, axis=0) / (2 * self.n_ind)
        self.geno = (geno - 2 * f) / np.sqrt(2 * f * (1 - f))

    def glasso_iso_cov(self):
        # read in isoform expression data
        # run glasso to get adjusted isoform exptression covariance matrix
        cur_tpm = pd.read_csv(self.gene + "_tpm.csv")
        cur_tpm.drop(cur_tpm.columns[[0]], axis=1, inplace=True)
        cur_tpm = np.transpose(cur_tpm.to_numpy())
        cur_tpm_std = (cur_tpm - np.mean(cur_tpm, axis=0)) / np.std(cur_tpm, axis=0)
        md = covariance.GraphicalLassoCV().fit(cur_tpm_std)
        self.iso_cov = md.covariance_

        # self.iso_cov = np.loadtxt(file)
        # print("isoform expression covariance matrix:")
        print(self.iso_cov)

    def sample_SNP(self, snp_causal, snp_ld):

        if snp_causal >= 1:
            self.n_causal = int(snp_causal)
        else:
            self.n_causal = round(self.n_SNP * snp_causal)

        # sample causal SNPs index
        # assumption1: same causal SNP for each isoforms
        self.causal_idx = np.zeros([self.n_causal, self.num_sim], dtype=int)
        for sim_i in range(self.num_sim):
            if snp_ld == 1: # low ld: random choice
                self.causal_idx[:, sim_i] =  np.random.choice(self.n_SNP, size=self.n_causal, replace=False)
            else: # high ld: choose blocks of SNPs
                rand = np.random.choice((self.n_SNP + 1 - self.n_causal))
                self.causal_idx[:, sim_i] = list(range(rand, rand + self.n_causal))

    def sample_beta_eps(self, h2g):

        # assumption2: every SNP contributes equally to iso expression heritability
        # assumption3: every iso contributes equally to gene expression heritability
        # epsilon takes (1 - h2g) of covariance term
        eps_cov = (1 - h2g) * self.iso_cov
        beta_cov = h2g * self.iso_cov
        beta_cov = beta_cov / self.n_causal
        self.h2g_iso = h2g
        print("heritability for each isoform: " + str(self.h2g_iso))
        print("covariance matrix of beta:")
        print(beta_cov)
        print("covariance matrix of epsilon:")
        print(eps_cov)

        # # check positive semi-definite
        # if not np.all(np.linalg.eigvals(beta_cov) >= 0):
        #     os.rmdir(outdir)
        #     sys.exit("covariance matrix of beta is not positive-semidefinite")

        self.beta = np.zeros([self.n_causal, self.num_iso, self.num_sim], dtype=float)
        self.eps = np.zeros([self.n_ind, self.num_iso, self.num_sim], dtype=float)
        self.fail = np.zeros([self.num_sim])
        for sim_i in range(self.num_sim):

            # if not np.all(np.linalg.eigvals(eps_cov) >= 0):
            #     print("simulation" + str(sim_i) + "fails.")
            #     self.fail[sim_i] = 1
            #     continue

            for ind_j in range(self.n_ind):
                self.eps[ind_j, :, sim_i] = np.random.multivariate_normal(np.zeros(self.num_iso), eps_cov)
            for snp_j in range(self.n_causal):
                self.beta[snp_j, :, sim_i] = np.random.multivariate_normal(np.zeros(self.num_iso),
                                                                           beta_cov)  # simulate beta

    def sample_iso_expr(self, outdir):
        # geno = np.loadtxt(self.gene + "_AFR_sample.txt", delimiter=",", dtype=float)

        # simulate isoform expression data
        # exp = geno * beta + eps
        exp = np.zeros([self.n_ind, self.num_iso])
        g_full = np.zeros([self.n_ind, self.num_iso])
        heritability = np.zeros(self.num_sim)
        i = 0
        for sim_i in range(self.num_sim):
            # if self.fail[sim_i] == 1:
            #     continue
            cur_geno = self.geno[:, self.causal_idx[:, sim_i]]
            for iso_k in range(self.num_iso):

                cur_beta = self.beta[:, iso_k, sim_i]
                g = np.matmul(cur_geno, cur_beta)
                e = self.eps[:, iso_k, sim_i]
                g_full[:, iso_k] = g
                exp[:, iso_k] = g + e
            np.savetxt(outdir + "/" + self.gene + "_sim" + str(i) + ".csv", exp, delimiter=",")

            # sanity check
            exp_sum = np.sum(exp, axis=1)
            g_sum = np.sum(g_full, axis=1)
            # print("heritability:" + str(np.var(g_sum) / np.var(exp_sum)))
            # print(np.cov(np.transpose(g_full)))
            # print(np.cov(np.transpose(exp)))
            heritability[sim_i] = np.var(g_sum) / np.var(exp_sum)
            i = i + 1
        self.num_success = i
        print("In total, " + str(i) + " simulations succeed.")
        np.savetxt(outdir + "/h2g_real.csv", heritability, delimiter=",")
        print(np.mean(heritability))
        print(np.var(heritability))

    def sanity_check(self, outdir):
        # check whether isoform expression covariance matrix of the simulated data
        # is similar to the real matrix

        sim_cov = np.zeros([self.num_success, self.num_iso * self.num_iso])
        for sim_i in range(self.num_success):
            exp = np.loadtxt(outdir + "/" + self.gene + "_sim" + str(sim_i) + ".csv", delimiter=",")
            cov = np.cov(np.transpose(exp))
            sim_cov[sim_i, :] = np.reshape(cov, self.num_iso * self.num_iso)

        sim_cov_mean = np.mean(sim_cov, axis = 0)
        sim_cov_mean = np.reshape(sim_cov_mean, (self.num_iso, self.num_iso))
        print(sim_cov_mean)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='eQTL Isoform Heritability Simulation')

    parser.add_argument('--gene', type=str,
                        help="gene name (start with ENSG0000)")
    parser.add_argument('--isoform', type=int,
                        help="number of mRNA isoforms from the gene: 2-10")

    parser.add_argument('--sim', type=int, default=1000,
                        help="number of simulations")

    parser.add_argument("--h2g", type=float, default=0.1,
                        help="controlled gene heritability")
    parser.add_argument("--snp_causal", type=float, default=1,
                        help="number of causal SNPs: 1, 2, 0.1, 0.2")
    parser.add_argument("--snp_ld", type=int, default=1,
                        help="low 1: completely random; high 2: block of SNPs")

    args = parser.parse_args()

    os.chdir("./" + args.gene)
    outdir = "../iso_exp"
    #os.chdir("/u/scratch/r/roserao/data/geno/" + str(args.isoform) + "/" + args.gene)
    #outdir = "/u/scratch/r/roserao/data/iso_exp/1"

    model = BaseModel(args.gene, args.isoform, args.sim)
    model.glasso_iso_cov()
    model.sample_SNP(args.snp_causal, args.snp_ld)
    model.sample_beta_eps(args.h2g)
    model.sample_iso_expr(outdir)
    model.sanity_check(outdir)










