import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import numpy as np

def mad(x, y):
    # make this correspond to a "consistent" version of the estimator https://en.wikipedia.org/wiki/Median_absolute_deviation#Relation_to_standard_deviation
    return 1.4826 * np.abs(np.log2(x+1) - np.log2(y+1)).median()

def main():
    d = {'sal_bc' : {}, 'kal_bc' : {}}
    t = pd.read_table('seqc_samples.tsv')
    for sn in set(t['id']):
        # these are for sample A
        if sn.startswith('A'):
            continue
        d['sal_bc'][sn] = pd.read_table(sn + '/quant_bc/quant.sf').set_index('Name')['TPM'].values
        d['kal_bc'][sn] = pd.read_table(sn + '/kquant_bc/abundance.tsv').set_index('target_id')['tpm'].values

    import itertools
    dat_s = []
    dat_k = []
    min_dat = np.inf
    max_dat = -np.inf
    for k in ['sal_bc', 'kal_bc']:
        x = pd.DataFrame(d[k])
        for s1, s2 in itertools.combinations(x.columns, 2):
            a = x.loc[:,s1]
            b = x.loc[:,s2]
            v = mad(a,b)
            if v > max_dat:
                max_dat = v
            if v < min_dat:
                min_dat = v
            kind = "within" if s1.split('_')[1] == s2.split('_')[1] else "across"
            if k == 'sal_bc':
                dat_s.append((s1, s2, v))
            elif k == 'kal_bc':
                dat_k.append((s1, s2, v))

    df1 = pd.DataFrame(dat_s, columns=['center 1', 'center 2', 'mad']).pivot('center 1', 'center 2', 'mad')
    df2 = pd.DataFrame(dat_k, columns=['center 1', 'center 2', 'mad']).pivot('center 1', 'center 2', 'mad')

    ax = matplotlib.pyplot.axes()
    hms = sns.heatmap(df1, vmin=min_dat, vmax=max_dat, fmt='.2f', annot=True, linewidths=.1, ax=ax)
    ax.set_title("Salmon")
    matplotlib.pyplot.tight_layout()
    hms.get_figure().savefig("salmon_hmap_median.png")
   
    matplotlib.pyplot.clf()
    matplotlib.pyplot.cla()
    ax = matplotlib.pyplot.axes()
    hmk = sns.heatmap(df2, vmin=min_dat, vmax=max_dat, fmt='.2f', annot=True, linewidths=.1, ax=ax)
    ax.set_title("kallisto")
    matplotlib.pyplot.tight_layout()
    hmk.get_figure().savefig("kallisto_hmap_median.png")


if __name__ == "__main__":
    main()
