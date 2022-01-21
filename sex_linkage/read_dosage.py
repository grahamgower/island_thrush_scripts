#!/usr/bin/env python3

import sys
import os.path
import numpy as np
from sklearn import decomposition, cluster, linear_model
from scipy.stats import binom, chi2, pearsonr

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

def parse_seqtk_comp(fn, min_length):
    c_len = []
    c_GC = []
    with open(fn) as f:
        for line in f:
            fields = line.split()
            contig = fields[0]
            length = int(fields[1])
            if length < min_length:
                continue
            A, C, G, T = [int(x) for x in fields[2:6]]
            GC = (C+G) / (A+C+G+T)

            c_len.append(length)
            c_GC.append(GC)
    return np.array(c_len, dtype=int), np.array(c_GC)

def parse_idxstats1(fn, min_length):
    n = []
    L = []
    contigs = []
    with open(fn) as f:
        for line in f:
            fields = line.split()
            contig = fields[0]
            length = int(fields[1])
            if length < min_length:
                continue
            n_mapped = int(fields[2])
            n.append(n_mapped)
            L.append(length)
            contigs.append(contig)
    return np.array(n, dtype=int), L, contigs

def parse_idxstats(filelist, min_length, min_reads):
    samples = []
    N = []
    L = None
    for fn in filelist:
        sample = os.path.basename(fn)
        if sample.endswith(".txt"):
            sample = sample[:-len(".txt")]
        if sample.endswith(".idxstats"):
            sample = sample[:-len(".idxstats")]
        n, L2, contigs = parse_idxstats1(fn, min_length)
        if np.sum(n) < min_reads:
            continue
        if L is None:
            L = L2
        else:
            if np.any(L != L2):
                raise Exception("{}: reference mismatch".format(fn))
        samples.append(sample)
        N.append(n)
    return (np.array(samples, dtype=str), np.array(contigs, dtype=str),
            np.array(N, dtype=int), np.array(L, dtype=int))

def chr_predict(N, L, eps, GC=None):
    M = N / np.sum(N,axis=1).reshape(-1,1) / (L / np.sum(L))

    if GC is not None:
        # To remove variation between contigs due to GC content,
        # regress observations against GC, and do PCA on the residuals.
        GC2 = GC.reshape(-1,1)
        for row in range(N.shape[0]):
            lm = linear_model.LinearRegression()
            lm.fit(GC2, M[row])
            M[row] -= lm.predict(GC2)

    M -= np.mean(M, axis=0) # mean centering
    M_t = np.transpose(M)
    pca = decomposition.PCA(n_components=2)
    pc = pca.fit_transform(M_t)

    cl = cluster.DBSCAN(eps)
    clusters = cl.fit_predict(pc)

    clusters[np.where(np.bitwise_and(clusters!=0, clusters!=1))] = 2

    return pc, clusters

def sex_predict(N, L, clusters):
    i_aut = np.where(clusters == 0)[0]
    i_sex = np.where(clusters == 1)[0]

    Laut = np.sum(L[i_aut])
    Lsex = np.sum(L[i_sex])
    Naut = np.sum(N[:,i_aut], axis=1)
    Nsex = np.sum(N[:,i_sex], axis=1)

    E1 = Lsex / (Lsex+2*Laut)     # 1 copy of sex chrom
    E2 = 2*Lsex / (2*Lsex+2*Laut) # 2 copies of sex chrom

    ll1 = binom.logpmf(Nsex, Nsex+Naut, E1)
    ll2 = binom.logpmf(Nsex, Nsex+Naut, E2)

    Rx = Nsex / (Naut+Nsex)
    Rx[np.where(ll1>ll2)] /= 2*E1
    Rx[np.where(ll1<=ll2)] /= E2

    sex = np.zeros(len(Rx), dtype=int)
    alpha = 0.001
    sex[np.where(chi2.sf(2*(ll1-ll2), 1) < alpha)] = 1 # 1 copy
    sex[np.where(chi2.sf(2*(ll2-ll1), 1) < alpha)] = 2 # 2 copies

    # Suspicious samples, they probably violate the model in some way.
    sex[np.where(np.bitwise_and(Rx>0.6, Rx<0.8))] = 0

    return Rx, sex

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Read dosage analyser.")
    parser.add_argument("-l", "--min-length", type=int, default=int(1e5),
            help="Exclude contigs with length shorter than this [%(default)s].")
    parser.add_argument("-n", "--min-reads", type=int, default=int(1e6),
            help="Exclude samples with fewer reads than this [%(default)s].")
    parser.add_argument("-p", "--plot-file", type=str,
            help="Plot contig PCA and read dosage histogram to specified file.")
    parser.add_argument("-a", "--plot-all-hist", default=False, action="store_true",
            help="Plot separate read dosage histograms for each sample [%(default)s].")
    parser.add_argument("-S", "--sex-chromosomes", default="XY",
            help="The sex chromosome system. XY, ZW, etc.")
    parser.add_argument("-s", "--samples", metavar="samples.txt", default=None,
            help="Output file for sample table. Tab separated fields are:"
                    "sample, Rx, sex.")
    parser.add_argument("-c", "--contigs", metavar="contigs.txt", default=None,
            help="Output file for contig table. Tab separated fields are:"
                    "contig, PC1, PC2, chrom-assignment.")
    parser.add_argument("-e", "--DBSCAN-eps", default=0.2, type=float,
            help="DBSCAN clustering epsilon parameter.")
    parser.add_argument("-q", "--seqtk_fn", metavar="seqtk-comp.txt",
            help="nucleotide composition of reference, as output by `seqtk comp'.")
    parser.add_argument("-t", "--title", help="Plot title.")
    parser.add_argument("infiles", metavar="*.idxstats", nargs="+",
            help="Input files, as output by `samtools idxstats'.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    samples, contigs, N, L = parse_idxstats(args.infiles, args.min_length, args.min_reads)

    GC = None
    if args.seqtk_fn is not None:
        L2, GC = parse_seqtk_comp(args.seqtk_fn, args.min_length)
        if np.any(L != L2):
            raise Exception("{}: reference mismatch with idxstats files".format(args.seqtk_fn))

    pc, clusters = chr_predict(N, L, args.DBSCAN_eps, GC)
    Rx, sex = sex_predict(N, L, clusters)

    if args.contigs is not None:
        contigtab = {0:"A", 1:args.sex_chromosomes[0]}
        with open(args.contigs, "w") as f:
            print("contig", "pc1", "pc2", "chrom", sep="\t", file=f)
            for contig, (pc1, pc2), i in zip(contigs, pc, clusters):
                print(contig, pc1, pc2, contigtab.get(i, "U"), sep="\t", file=f)

    if args.samples is not None:
        sextab = {0:"U", 1:args.sex_chromosomes,
                2:args.sex_chromosomes[0]+args.sex_chromosomes[0]}
        with open(args.samples, "w") as f:
            print("sample", "Rx", "sex", sep="\t", file=f)
            for sample, r, x in zip(samples, Rx, sex):
                print(sample, r, sextab[x], sep="\t", file=f)

    if args.plot_file is None:
        exit(0)

    pdf = PdfPages(args.plot_file)
    fig_w, fig_h = plt.figaspect(9.0/16.0)

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    markers = "oX1"
    chr_labels = ["Autosomal", "Chr "+args.sex_chromosomes[0], "Unassigned"]

    for i in (0,1,2):
        j = np.where(clusters == i)[0]
        if len(j) == 0:
            continue
        x, y = zip(*pc[j])
        ax1.scatter(x, y, edgecolor="black", lw=.5, marker=markers[i], label=chr_labels[i])


    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")
    if args.title:
        ax1.set_title(args.title)
    else:
        ax1.set_title("PCA (contig read dosage)")
    ax1.legend()

    plt.tight_layout()
    pdf.savefig(figure=fig1)
    plt.close(fig1)

    if args.plot_all_hist:
        fig2 = plt.figure(figsize=(fig_w, fig_h))
        gs2 = gridspec.GridSpec(1, 1)
        ax2 = fig2.add_subplot(gs2[0])

        ax2.hist(Rx, bins=40)

        xlim = ax2.get_xlim()
        ylim = ax2.get_ylim()

        ax2.vlines(0.5, ylim[0], ylim[1]*2, linestyle=':')
        ax2.vlines(1.0, ylim[0], ylim[1]*2, linestyle=':')
        ax2.set_ylim(ylim)
        ax2.set_xlim(min(0,xlim[0]), max(1.5,xlim[1]))
        ax2.set_xlabel("Read dosage")
        ax2.set_ylabel("Sample counts")
        ax2.set_title("Sample histogram (Chr {} read dosage)".format(args.sex_chromosomes[0]))

        plt.tight_layout()
        pdf.savefig(figure=fig2)
        plt.close(fig2)

    M = N / np.sum(N,axis=1).reshape(-1,1) / (L / np.sum(L))
    M_mean = np.mean(np.transpose(M), axis=1)
    #r, p = pearsonr(GC, M_mean)
    #print(r, p)

    #if args.seqtk_fn is not None:
    if False:

        fig3 = plt.figure(figsize=(fig_w, fig_h))
        gs3 = gridspec.GridSpec(1, 1)
        ax3 = fig3.add_subplot(gs3[0])

        for i in (0,1,2):
            j = np.where(clusters == i)[0]
            if len(j) == 0:
                continue
            #lm = linear_model.LinearRegression()
            #lm.fit(GC[j].reshape(-1,1), M_mean[j])
            #r2 = lm.score(GC[j].reshape(-1,1), M_mean[j])
            r, p = pearsonr(GC[j], M_mean[j])
            ax3.scatter(GC[j], M_mean[j], edgecolor="black", lw=.5, marker=markers[i], label=chr_labels[i] + "; $R^2={:.2g}$; $p={:.2g}$".format(r*r, p))

        ax3.set_xlabel("GC")
        ax3.set_ylabel("mean read dosage")
        ax3.legend()

        plt.tight_layout()
        pdf.savefig(figure=fig3)
        plt.close(fig3)

    if args.plot_all_hist:
        bins = np.linspace(0.4, 1.5, 100)

        for i, sample in enumerate(samples):
            fig = plt.figure(figsize=(fig_w, fig_h))
            gs = gridspec.GridSpec(1, 1)
            ax = fig.add_subplot(gs[0])

            ax.hist(M[i], bins=bins, density=True)

            ax.set_xlim([0.4,1.5])
            ax.set_xlabel("Read dosage")
            ax.set_ylabel("log(Contig density)")
            ax.set_title("Histogram of $\\mathbf{{M}}_{{i}}$ for {})".format(sample))

            plt.yscale('log', nonposy='clip')

            plt.tight_layout()
            pdf.savefig(figure=fig)
            plt.close(fig)

    pdf.close()
