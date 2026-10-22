import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager as fm
try: fm.fontManager.addfont('/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf')
except Exception: pass
plt.rcParams['font.family']=['Liberation Sans','Arial','DejaVu Sans']
d=pd.read_csv('data/processed/ASD_module_trait.csv').set_index('Module')
order="ME31 ME4 ME17 ME19 ME8 ME10 ME27 ME13 ME23 ME7 ME22 ME2 ME20 ME0 ME26 ME25 ME15 ME16 ME32 ME6 ME28 ME3 ME11 ME12 ME5 ME14 ME30 ME29 ME1 ME9 ME21 ME18 ME24".split()
assert set(order)==set(d.index) and len(order)==33
r=d.loc[order,'Correlation'].values; p=d.loc[order,'Pvalue'].values
cmap=LinearSegmentedColormap.from_list('bwr',['#0B6BFF','#FFFFFF','#FF1A0B'])
def fmt_p(x):
    s='%.1g'%x
    return s.replace('e-0','e-')
fig=plt.figure(figsize=(4.08,6.12))
ax=fig.add_axes([0.17,0.08,0.52,0.90]); cax=fig.add_axes([0.80,0.30,0.04,0.45])
im=ax.imshow(r.reshape(-1,1),cmap=cmap,vmin=-1,vmax=1,aspect='auto')
for i,(ri,pi) in enumerate(zip(r,p)):
    ax.text(0,i,f'{ri:.2g} (p = {fmt_p(pi)})',ha='center',va='center',fontsize=6.5)
ax.set_yticks(range(33)); ax.set_yticklabels(order,fontsize=7)
ax.set_xticks([0]); ax.set_xticklabels(['ASD diagnosis'],fontsize=8)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_linewidth(0.8)
cb=fig.colorbar(im,cax=cax); cb.ax.tick_params(labelsize=7); cb.set_label('Pearson r',fontsize=7.5)
fig.savefig('figures/export_600dpi/Figure3_module_trait_heatmap.png',dpi=600,facecolor='white')
