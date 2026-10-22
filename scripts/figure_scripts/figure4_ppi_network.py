import pandas as pd, networkx as nx, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
e=pd.read_csv('results/ppi/ME11_ppi_network.csv')
G=nx.Graph()
for a,b,s in e.itertuples(index=False):
    if a!=b: G.add_edge(a,b,score=max(s,G[a][b]['score']) if G.has_edge(a,b) else s)
print(G.number_of_nodes(),G.number_of_edges())
hubs={'TMSB4X','RAP1B','C1D','ZCCHC10'}
import numpy as np
comps=sorted(nx.connected_components(G),key=len,reverse=True)
pos={}
gc=G.subgraph(comps[0]); p=nx.kamada_kawai_layout(gc)
xy=np.array(list(p.values())); mn,mx=xy.min(0),xy.max(0)
for n,(x,y) in p.items(): pos[n]=((x-mn[0])/(mx[0]-mn[0])*0.70,(y-mn[1])/(mx[1]-mn[1]))
small=comps[1:]; cols=2; cell=1.0/ ((len(small)+cols-1)//cols)
for i,c in enumerate(small):
    sg=G.subgraph(c); p=nx.kamada_kawai_layout(sg) if len(c)>2 else {n:(j*1.0-0.5,0) for j,n in enumerate(c)}
    cx=0.82+0.13*(i%cols); cy=1-cell*(i//cols)-cell/2
    for n,(x,y) in p.items(): pos[n]=(cx+x*0.045,cy+y*0.045*cell*6)
print(len(comps),[len(c) for c in comps][:5])
plt.rcParams['font.family']='DejaVu Sans'
fig,ax=plt.subplots(figsize=(9,7.5),dpi=200)
hi=[(u,v) for u,v,d in G.edges(data=True) if d['score']>=700]; lo=[(u,v) for u,v,d in G.edges(data=True) if d['score']<700]
nx.draw_networkx_edges(G,pos,edgelist=lo,edge_color='#C9CED3',width=0.6,ax=ax)
nx.draw_networkx_edges(G,pos,edgelist=hi,edge_color='#7F8C8D',width=1.1,ax=ax)
other=[n for n in G if n not in hubs]
nx.draw_networkx_nodes(G,pos,nodelist=other,node_color='#85C1E9',node_size=60,linewidths=0,ax=ax)
nx.draw_networkx_nodes(G,pos,nodelist=[n for n in hubs if n in G],node_color='#E74C3C',node_size=260,edgecolors='black',linewidths=0.6,ax=ax)
gc_nodes=list(comps[0]); cxy=np.mean([pos[n] for n in gc_nodes],axis=0)
manual={'RAP1B':(-0.16,-0.10),'C1D':(0.12,-0.13),'ZCCHC10':(0.13,0.10),'TMSB4X':(-0.05,0.10)}
for n in hubs:
    if n in G:
        x,y=pos[n]; dx,dy=manual[n]
        ax.annotate(n,(x,y),xytext=(x+dx,y+dy),ha='center',va='center',fontsize=10,fontweight='bold',zorder=10,
                    arrowprops=dict(arrowstyle='-',lw=0.8,color='black',shrinkA=2,shrinkB=7),
                    bbox=dict(fc='white',ec='#888888',lw=0.5,boxstyle='round,pad=0.25'))

from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],color='#C9CED3',lw=1,label='score 400–699'),Line2D([0],[0],color='#7F8C8D',lw=1.5,label='score ≥ 700')],loc='lower right',frameon=False,fontsize=9)
ax.axis('off'); fig.tight_layout(); fig.savefig('figures/export_600dpi/Figure4_ME11_PPI_network.png',dpi=360); 
