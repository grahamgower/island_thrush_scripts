#!/usr/bin/env python3

import sys
import numpy as np
from sklearn import decomposition, cluster

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

def parse_fai(fn):
    with open(fn) as f:
        for line in f:
            fields = line.split()
            chrom = fields[0]
            size = int(fields[1])
            if chrom.startswith("chr"):
                chrom = chrom[len("chr"):]
            yield chrom, size

LNK_AUT = 0
LNK_SEX1 = 1
LNK_SEX2 = 2
def sim_contigs(fn, sex_chroms=(82.5e6, 7.5e6), minlen=int(1e5)):
    """
    sex_chroms -- Tuple of the two sex chromosome lengths. The first should
                  be the chromosome found in the homogametic sex
                  (e.g. X in mammals, Z in birds).
    """
    empir_contigs = np.array([sz for _,sz in parse_fai(fn)], dtype=int)
    n_contigs = len(empir_contigs)
    genome_len = np.sum(empir_contigs)
    p_sex1 = sex_chroms[0] / genome_len
    p_sex2 = sex_chroms[1] / genome_len

    sim_contigs = empir_contigs[np.random.randint(n_contigs, size=n_contigs)]
    x = np.random.random(size=n_contigs)
    chr_type = np.full(n_contigs, LNK_AUT, dtype=int) # autosomal chrom
    chr_type[np.where(x < p_sex1)] = LNK_SEX1 # linked to sex chrom 1
    chr_type[np.where(x > 1-p_sex2)] = LNK_SEX2 # linked to sex chrom 2

    ifilt = np.where(sim_contigs > minlen)

    return sim_contigs[ifilt], chr_type[ifilt]

def sim_read_counts(contigs, chr_type,
                    n_reads=1e5, n_inds=74, sex_ratio=50/(50+22), hom_ref=False):
    """
    sex_ratio -- The proportion of individuals that have two copies of the
                 first sex chromosome.
    hom_ref   -- The reference is homogametic, so all reads from the
                 second sex chromosome map to the first.
    """

    n_contigs = len(contigs)

    j_aut = np.where(chr_type == 0)
    j_sc1 = np.where(chr_type == 1)
    j_sc2 = np.where(chr_type == 2)

    glen_sex1 = 2*np.sum(contigs[j_aut]) + 2*np.sum(contigs[j_sc1])
    glen_sex2 = (2*np.sum(contigs[j_aut]) + np.sum(contigs[j_sc1])
                    + np.sum(contigs[j_sc2]))

    # Build multinomial probability vectors.
    p_mn_sex1 = np.empty(n_contigs)
    p_mn_sex2 = np.empty(n_contigs)

    p_mn_sex1[j_aut] = 2 *contigs[j_aut] / glen_sex1
    p_mn_sex2[j_aut] = 2 *contigs[j_aut] / glen_sex2

    p_mn_sex1[j_sc1] = 2 *contigs[j_sc1] / glen_sex1
    p_mn_sex2[j_sc1] = contigs[j_sc1] / glen_sex2

    p_mn_sex1[j_sc2] = 0
    if hom_ref:
        p_mn_sex2[j_sc2] = 0
    else:
        p_mn_sex2[j_sc2] = contigs[j_sc2] / glen_sex2

    p_mn_sex1 /= np.sum(p_mn_sex1)
    p_mn_sex2 /= np.sum(p_mn_sex2)

    sex = np.random.binomial(1, sex_ratio, size=n_inds)
    n_sex1  = np.sum(sex)
    n_sex2  = len(sex) - n_sex1

    # Draw read counts.
    counts = np.empty((n_inds, n_contigs))
    counts[np.where(sex==1)] = np.random.multinomial(n_reads, p_mn_sex1,
                                                size=n_sex1)
    counts[np.where(sex==0)] = np.random.multinomial(n_reads, p_mn_sex2,
                                                size=n_sex2)

    return np.array(counts, dtype=int), sex

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} contig-ref.fai out.pdf".format(sys.argv[0]))
        exit(1)

    fai_fn = sys.argv[1]
    out_fn = sys.argv[2]

    contigs, chr_lnk = sim_contigs(fai_fn)
    counts, sex = sim_read_counts(contigs, chr_lnk)

    # Filter out contigs with no reads mapped to them.
    contig_filt = np.where(np.sum(counts,axis=0) > 0)[0]
    counts = counts[:,contig_filt]
    contigs = contigs[contig_filt]
    chr_lnk = chr_lnk[contig_filt]

    tot_reads_per_sample = np.sum(counts,axis=1).reshape(counts.shape[0],1)
    Rx = (counts / tot_reads_per_sample) / (contigs / np.sum(contigs))
    Rx -= np.mean(Rx, axis=0) # mean centering

    pca = decomposition.PCA(n_components=2)
    pc = pca.fit_transform(np.transpose(Rx))

    pdf = PdfPages(out_fn)
    fig_w, fig_h = plt.figaspect(9.0/16.0)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])

    markers = "oXd1"
    labels = ["Autosomal", "Sex chr1", "Sex chr2", "Outlier"]

    show_truth=False
    if show_truth:
        for lnk in (LNK_AUT, LNK_SEX1, LNK_SEX2):
            d = pc[np.where(chr_lnk==lnk)]
            if len(d) == 0:
                continue
            x, y = zip(*d)
            ax.scatter(x, y, edgecolor="black", lw=.5, marker=markers[lnk], label=labels[lnk])

    else:
        #cl = cluster.KMeans(n_clusters=2)
        cl = cluster.DBSCAN() #eps=0.1, min_samples=5)
        clusters = cl.fit_predict(pc)
        for i in np.unique(clusters):
            d = pc[np.where(clusters == i)]
            if len(d) == 0:
                continue
            x, y = zip(*d)
            if i > 3:
                i = 3
            ax.scatter(x, y, edgecolor="black", lw=.5, marker=markers[i], label=labels[i])


    plt.legend()

    plt.tight_layout()
    pdf.savefig(figure=fig)
    pdf.close()
