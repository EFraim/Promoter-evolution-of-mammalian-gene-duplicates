import pickle
import numpy as np
from row_matchers import one_to_one_matches
from pandas import concat
from scipy.stats import mannwhitneyu
export=pickle.load(open('/data/db/import/save/mouse-export.pkl', 'rb'))
export.loc[export["Branch point"] == "Mus musculus strain reference (CL57BL6)", "Branch point"] = "Mus musculus"

branch_order = ["Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Gnathostomata", "Euteleostomi", "Sarcopterygii", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Boreoeutheria", "Euarchontoglires", "Glires", "Rodentia", "Myomorpha", "Muroidea", "Muridae", "Murinae", "Mus", "Mus musculus"]
exprcols=[c for c in export.columns if c.endswith('-tpm')]
exprnorm = export.assign(avgexpr = np.log(export[exprcols].mean(axis=1)))[export['mut-pct'].notna()]


exprbranchnormed=concat([
    concat(concat([less, full]) for less, full in (one_to_one_matches(exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 0)], exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 2)], 'avgexpr', np.log(2.0)) for b in branch_order)),
    concat(concat([less, full]) for less, full in (one_to_one_matches(exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 0)], exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 1)], 'avgexpr', np.log(2.0)) for b in branch_order)),
    concat(concat([less, full]) for less, full in (one_to_one_matches(exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 1)], exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 2)], 'avgexpr', np.log(2.0)) for b in branch_order))])

exprdsnormed=concat([
    concat(concat([less, full]) for less, full in (one_to_one_matches(exprnorm[(exprnorm['dSbin'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 0)], exprnorm[(exprnorm['dSbin'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 2)], 'avgexpr', np.log(2.0)) for b in sorted(exprnorm['dSbin'].unique()))),
    concat(concat([less, full]) for less, full in (one_to_one_matches(exprnorm[(exprnorm['dSbin'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 0)], exprnorm[(exprnorm['dSbin'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 1)], 'avgexpr', np.log(2.0)) for b in sorted(exprnorm['dSbin'].unique()))),
    concat(concat([less, full]) for less, full in (one_to_one_matches(exprnorm[(exprnorm['dSbin'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 1)], exprnorm[(exprnorm['dSbin'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 2)], 'avgexpr', np.log(2.0)) for b in sorted(exprnorm['dSbin'].unique())))])

dsSkip = set()
branchSkip = set()

for x in sorted(exprdsnormed['dSbin'].unique()):
    if sum((exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 2)&(exprdsnormed['mut-pct'].notna()))>0 and sum((exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 0)&(exprdsnormed['mut-pct'].notna()))>0:
        print('Group sizes: %d, %d'%(sum((exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 2)&(exprdsnormed['mut-pct'].notna())), sum((exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 0)&(exprdsnormed['mut-pct'].notna()))))
        print("Mann-Whitney test of TFBM (abs) on dSbin %d CGI-full vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 2)]['mutual'], exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 0)]['mutual'], alternative='greater').pvalue*len(exprdsnormed['dSbin'].unique())))
        print("Mann-Whitney test of TFBM (abs) on dSbin %d CGI-full vs. mixed (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 2)]['mutual'], exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 1)]['mutual'], alternative='greater').pvalue*len(exprdsnormed['dSbin'].unique())))
        print("Mann-Whitney test of TFBM (abs) on dSbin %d mixed vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 1)]['mutual'], exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 0)]['mutual'], alternative='greater').pvalue*len(exprdsnormed['dSbin'].unique())))        
        print("Mann-Whitney test of TFBM (pct) on dSbin %d CGI-full vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 2)]['mut-pct'], exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 0)]['mut-pct'], alternative='greater').pvalue*len(exprdsnormed['dSbin'].unique())))
        print("Mann-Whitney test of TFBM (pct) on dSbin %d CGI-full vs. mixed (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 2)]['mut-pct'], exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 1)]['mut-pct'], alternative='greater').pvalue*len(exprdsnormed['dSbin'].unique())))
        print("Mann-Whitney test of TFBM (pct) on dSbin %d mixed vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 1)]['mut-pct'], exprdsnormed[(exprdsnormed['dSbin'] == x)&(exprdsnormed['CpG-ness'] == 0)]['mut-pct'], alternative='greater').pvalue*len(exprdsnormed['dSbin'].unique())))
    else:
        print('Skipping %d'%x)
        dsSkip.add(x)

for x in branch_order:
    if sum((exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 2)&(exprbranchnormed['mut-pct'].notna()))>0 and sum((exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 0)&(exprbranchnormed['mut-pct'].notna()))>0:
        print('Group sizes: %d, %d'%(sum((exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 2)&(exprbranchnormed['mut-pct'].notna())), sum((exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 0)&(exprbranchnormed['mut-pct'].notna()))))
        print("Mann-Whitney test of TFBM (abs) on Branch point %s CGI-full vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 2)]['mutual'].dropna(), exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 0)]['mutual'].dropna(), alternative='greater').pvalue*len(branch_order)))
        print("Mann-Whitney test of TFBM (abs) on Branch point %s CGI-full vs. mixed (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 2)]['mutual'].dropna(), exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 1)]['mutual'].dropna(), alternative='greater').pvalue*len(branch_order)))
        print("Mann-Whitney test of TFBM (abs) on Branch point %s mixed vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 1)]['mutual'].dropna(), exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 0)]['mutual'].dropna(), alternative='greater').pvalue*len(branch_order)))
        print("Mann-Whitney test of TFBM (pct) on Branch point %s CGI-full vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 2)]['mut-pct'].dropna(), exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 0)]['mut-pct'].dropna(), alternative='greater').pvalue*len(branch_order)))
        print("Mann-Whitney test of TFBM (pct) on Branch point %s CGI-full vs. mixed (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 2)]['mut-pct'].dropna(), exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 1)]['mut-pct'].dropna(), alternative='greater').pvalue*len(branch_order)))
        print("Mann-Whitney test of TFBM (pct) on Branch point %s mixed vs. CGI-less (mouse), Expression and FDR controlled: %s"%(x, mannwhitneyu(exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 1)]['mut-pct'].dropna(), exprbranchnormed[(exprbranchnormed['Branch point'] == x)&(exprbranchnormed['CpG-ness'] == 0)]['mut-pct'].dropna(), alternative='greater').pvalue*len(branch_order)))
    else:
        print('Skipping %s'%x)
        branchSkip.add(x)


import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

tfbm_scores = pickle.load(open('/data/db/import/save/mouse-tfbm-scores.pkl', 'rb'))
tfbm_scores_pct = pickle.load(open('/data/db/import/save/mouse-tfbm-scores-pct.pkl', 'rb'))

from seaborn import boxplot, despine
import matplotlib.pyplot as plt
from matplotlib import pyplot
import pandas as pd
fig, ax = pyplot.subplots(figsize=(15,8))
pyplot.xticks(rotation=45)
g=boxplot(ax=ax, x="Branch point", y="mutual", hue="CpG-ness", data=pd.concat([exprbranchnormed[(~exprbranchnormed['Branch point'].isin(branchSkip))&((exprbranchnormed['factors-g1'] >= 2) | (exprbranchnormed['factors-g2'] >= 2))], pd.DataFrame.from_dict({'mutual':tfbm_scores}).assign(**{'CpG-ness':3, 'Branch point':'baseline'})]), order=[*[b for b in branch_order if b not in branchSkip], 'baseline'], palette=['#0b9dcf', '#926AA5', '#ef893d', '#b080e0'], dodge=True, fliersize=0.0)
handles, _ = g.get_legend_handles_labels()
g.legend(handles, ['CGI-less', 'mixed', 'CGI-full', 'baseline'])
despine()
plt.savefig('mouse-tfbm-mutual-by-cgi-bp-with-baseline-expr-ctrl-mixed.pdf', dpi=600)
fig, ax = pyplot.subplots(figsize=(15,8))
pyplot.xticks(rotation=45)
g=boxplot(ax=ax, x="Branch point", y="mut-pct", hue="CpG-ness", data=pd.concat([exprbranchnormed[(~exprbranchnormed['Branch point'].isin(branchSkip))&((exprbranchnormed['factors-g1'] >= 2) | (exprbranchnormed['factors-g2'] >= 2))], pd.DataFrame.from_dict({'mut-pct':[x*100 for x in tfbm_scores_pct]}).assign(**{'CpG-ness':3, 'Branch point':'baseline'})]), order=[*[b for b in branch_order if b not in branchSkip], 'baseline'], palette=['#0b9dcf', '#926AA5', '#ef893d', '#b080e0'], dodge=True, fliersize=0.0)
handles, _ = g.get_legend_handles_labels()
g.legend(handles, ['CGI-less', 'mixed', 'CGI-full', 'baseline'])
despine()
plt.savefig('mouse-tfbm-mutpct-by-cgi-bp-with-baseline-expr-ctrl-mixed.pdf', dpi=600)
fig, ax = pyplot.subplots(figsize=(15,8))
g=boxplot(ax=ax, x="dSbin", y="mutual", hue="CpG-ness", data=pd.concat([exprdsnormed[(~exprdsnormed['dSbin'].isin(dsSkip))&((exprdsnormed['factors-g1'] >= 2) | (exprdsnormed['factors-g2'] >= 2))], pd.DataFrame.from_dict({'mutual':tfbm_scores}).assign(**{'CpG-ness':3, 'dSbin': 'baseline'})]), order=[*[b for b in sorted(exprdsnormed[((exprdsnormed['factors-g1'] >= 2) | (exprdsnormed['factors-g2'] >= 2))]['dSbin'].unique()) if b not in dsSkip], 'baseline'], palette=['#0b9dcf', '#926AA5', '#ef893d', '#b080e0'], dodge=True, fliersize=0.0)
handles, _ = g.get_legend_handles_labels()
g.legend(handles, ['CGI-less', 'mixed', 'CGI-full', 'baseline'])
despine()
plt.savefig('mouse-tfbm-mutual-by-cgi-ds-with-baseline-expr-ctrl-mixed.pdf', dpi=600)
fig, ax = pyplot.subplots(figsize=(15,8))
g=boxplot(ax=ax, x="dSbin", y="mut-pct", hue="CpG-ness", data=pd.concat([exprdsnormed[(~exprdsnormed['dSbin'].isin(dsSkip))&((exprdsnormed['factors-g1'] >= 2) | (exprdsnormed['factors-g2'] >= 2))], pd.DataFrame.from_dict({'mut-pct':[x*100 for x in tfbm_scores_pct]}).assign(**{'CpG-ness':3, 'dSbin': 'baseline'})]), order=[*[b for b in sorted(exprdsnormed[((exprdsnormed['factors-g1'] >= 2) | (exprdsnormed['factors-g2'] >= 2))]['dSbin'].unique()) if b not in dsSkip], 'baseline'], palette=['#0b9dcf', '#926AA5', '#ef893d', '#b080e0'], dodge=True, fliersize=0.0)
handles, _ = g.get_legend_handles_labels()
g.legend(handles, ['CGI-less', 'mixed', 'CGI-full', 'baseline'])
despine()
plt.savefig('mouse-tfbm-mutpct-by-cgi-ds-with-baseline-expr-ctrl-mixed.pdf', dpi=600)
