import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Polygon, Ellipse, PathPatch
from matplotlib.path import Path
from matplotlib import font_manager as fm
fm.fontManager.addfont('/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf')
fm.fontManager.addfont('/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf')
plt.rcParams['font.family']='Liberation Serif'
W,H=1600,980
fig=plt.figure(figsize=(6.3,6.3*H/W)); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,W); ax.set_ylim(H,0); ax.axis('off')
PT=6.3*72/W   # points per layout unit
BG='#FDFDFE'; BD='#DEDEE0'; INK='#1E1E1E'; ARW='#232323'
rng=np.random.default_rng(11)
def box(x0,y0,x1,y1,fc=BG,ec=BD,lw=2): ax.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fc=fc,ec=ec,lw=lw*PT*1.4,zorder=1))
def text(cx,y,lines,size=22,bold=False):
    for i,l in enumerate(lines): ax.text(cx,y+i*size*1.25,l,ha='center',va='top',fontsize=size*PT,fontweight='bold' if bold else 'normal',color=INK,zorder=5)
def vline(x,y0,y1,arrow=True):
    ax.plot([x,x],[y0,y1-(14 if arrow else 0)],color=ARW,lw=3*PT,solid_capstyle='butt',zorder=3)
    if arrow: ax.add_patch(Polygon([(x,y1),(x-8,y1-16),(x+8,y1-16)],fc=ARW,ec='none',zorder=3))
def hline(x0,x1,y):
    ax.plot([x0,x1-14],[y,y],color=ARW,lw=3*PT,solid_capstyle='butt',zorder=3)
    ax.add_patch(Polygon([(x1,y),(x1-16,y-8),(x1-16,y+8)],fc=ARW,ec='none',zorder=3))
# ---------- icons ----------
def tube(cx,cy,h,fill,dark):
    w=h*0.30; x0=cx-w/2; y0=cy-h/2
    ax.add_patch(FancyBboxPatch((x0,y0+h*0.08),w,h*0.92,boxstyle=f"round,pad=0,rounding_size={w/2}",fc='white',ec='#5E6670',lw=1.6*PT*1.4,zorder=3))
    ax.add_patch(FancyBboxPatch((x0+w*0.1,y0+h*0.30),w*0.8,h*0.66,boxstyle=f"round,pad=0,rounding_size={w*0.4}",fc=fill,ec='none',zorder=3))
    ax.add_patch(Rectangle((x0+w*0.13,y0+h*0.32),w*0.28,h*0.5,fc='white',ec='none',alpha=0.85,zorder=4))
    ax.add_patch(Rectangle((x0-w*0.08,y0),w*1.16,h*0.1,fc='#5E6670',ec='none',zorder=4))
    ax.add_patch(Rectangle((x0+w*0.1,y0+h*0.30),w*0.8,h*0.035,fc=dark,ec='none',zorder=4))
def microarray(cx,cy,s):
    P=np.array([[-0.9,0.25],[0.25,-0.35],[0.95,-0.15],[-0.2,0.45]])*s+[cx,cy]
    ax.add_patch(Polygon(P,fc='#E9EBEE',ec='#8A9099',lw=1.5*PT*1.4,zorder=3))
    ax.add_patch(Polygon(P+[0,0.07*s],fc='#B9BEC5',ec='none',zorder=2))
    cols=['#5BC0BE','#3A6EA5','#8E6CC2','#E07A5F','#E9A23B','#66B3A6','#4E79A7']
    for i,c in enumerate(cols):
        x=cx-0.45*s+i*0.15*s; t=np.linspace(0,1,40)
        ax.plot(x+0.06*s*np.sin(t*9+i),cy-0.55*s+t*0.75*s,color=c,lw=2.6*PT*1.4,solid_capstyle='round',zorder=4)
def network(cx,cy,wd,ht,n=30,hub=True):
    gx,gy=np.meshgrid(np.linspace(-1,1,10),np.linspace(-1,1,4)); g=np.c_[gx.ravel(),gy.ravel()]
    g=g[rng.choice(len(g),n,replace=False)]+rng.uniform(-0.08,0.08,(n,2))
    pts=g*[wd/2,ht/2]+[cx,cy]
    for i in range(n):
        d=np.hypot(*(pts-pts[i]).T); nb=np.argsort(d)[1:3]
        for j in nb: ax.plot(*pts[[i,j]].T,color='#4A4A4A',lw=1*PT*1.4,zorder=3)
    for k,p in enumerate(pts):
        ax.add_patch(Ellipse(p,13,10,fc='white',ec='#3A3A3A',lw=1*PT*1.4,zorder=4))
    if hub:
        k=np.argmin(np.hypot(*(pts-[cx,cy]).T)); c=pts[k]
        for j in np.argsort(np.hypot(*(pts-c).T))[1:5]: ax.plot(*np.array([c,pts[j]]).T,color='#4A4A4A',lw=1*PT*1.4,zorder=3)
        ax.add_patch(Ellipse(c,15,12,fc='#C94C4C',ec='#7A2626',lw=1*PT*1.4,zorder=5))
def ppi(cx,cy,s):
    hub=(cx,cy); cl=[((-0.55,-0.25),'#F4C6E3'),((0.0,-0.55),'#C6EBE3'),((0.5,-0.25),'#DCEFC8'),((0.45,0.35),'#F6F0C8'),((-0.1,0.45),'#F2DDB6')]
    for (dx,dy),col in cl:
        ccx,ccy=cx+dx*s,cy+dy*s; ax.add_patch(Circle((ccx,ccy),0.2*s,fc=col,ec='none',alpha=0.8,zorder=2))
        nodes=[(ccx+0.12*s*np.cos(a),ccy+0.12*s*np.sin(a)) for a in np.linspace(0,2*np.pi,6,endpoint=False)]
        for i in range(6): ax.plot([nodes[i][0],nodes[(i+2)%6][0]],[nodes[i][1],nodes[(i+2)%6][1]],ls=(0,(2,2)),color='#8A8A8A',lw=0.8*PT*1.4,zorder=3)
        for k in (0,2,4): ax.plot([hub[0],nodes[k][0]],[hub[1],nodes[k][1]],color='#7A7A7A',lw=0.9*PT*1.4,zorder=3)
        for nx_,ny_ in nodes: ax.add_patch(Circle((nx_,ny_),0.032*s,fc=col,ec='#9A9A9A',lw=0.8*PT*1.4,zorder=4))
    ax.add_patch(FancyBboxPatch((cx-0.08*s,cy-0.03*s),0.16*s,0.06*s,boxstyle=f"round,pad=0,rounding_size={0.03*s}",fc='#D9534F',ec='none',zorder=5))
def adhd_panel(cx,cy,s):
    tube(cx-0.38*s,cy+0.05*s,0.62*s,'#2C4A7A','#1D3557')
    for i,x in enumerate([cx+0.02*s,cx+0.32*s]):
        ax.plot([x+0.08*s]*2,[cy-0.42*s,cy+0.42*s],color='#555',lw=1.3*PT*1.4,zorder=3)
        for yy in (cy-0.42*s,cy+0.42*s): ax.plot([x+0.05*s,x+0.11*s],[yy,yy],color='#555',lw=1.3*PT*1.4,zorder=3)
        ax.add_patch(Rectangle((x,cy-0.2*s),0.16*s,0.38*s,fc='#C3D7F1',ec='#4D7FC4',lw=1.3*PT*1.4,zorder=4))
        ax.plot([x,x+0.16*s],[cy-0.02*s+i*0.01*s]*2,color='#2F6DB5',lw=2*PT*1.4,zorder=5)
def splice(cx,cy,s):
    x=np.linspace(-1,1,300)*s*0.85
    for side,col in [(-1,'#3BB8C4'),(1,'#7ACC3E')]:
        m=(x*side>0.18*s); xx=cx+x[m]; yy=cy+0.06*s*np.sin(x[m]/s*14)
        ax.plot(xx,yy,color=col,lw=2.6*PT*1.4,zorder=3)
        for k in range(0,len(xx),9): ax.plot([xx[k],xx[k]],[yy[k]-6,yy[k]+6],color=col,lw=1*PT*1.4,zorder=3)
    for dx,dy,r in [(-0.12,0.02,0.13),(0.05,-0.05,0.15),(0.14,0.08,0.12),(-0.02,0.12,0.11),(0.0,0.0,0.12)]:
        ax.add_patch(Circle((cx+dx*s,cy+dy*s),r*s,fc='#A9BCC4',ec='#7F949D',lw=0.8*PT*1.4,alpha=0.9,zorder=4))
    t=np.linspace(0,2*np.pi,60)
    for dx,dy in [(-0.12,-0.16),(0.06,-0.24),(0.22,-0.12)]:
        ax.plot(cx+dx*s+0.07*s*np.cos(t),cy+dy*s+0.07*s*np.sin(t),color='#D33F49',lw=2*PT*1.4,zorder=5)
# ---------- layout (identical to current graphical abstract) ----------
y0,y1=40,270; b1=(150,y0,370,y1); b2=(450,y0,690,y1); b3=(770,y0,1450,y1)
for b in (b1,b2,b3): box(*b)
tube((b1[0]+b1[2])/2,y0+72,112,'#B83A3A','#8E2626'); text((b1[0]+b1[2])/2,y0+150,["ASD","(n=238)"])
microarray((b2[0]+b2[2])/2,y0+95,90); text((b2[0]+b2[2])/2,y0+165,["Microarray"])
cx3=(b3[0]+b3[2])/2; network(cx3,y0+78,360,100,n=30)
ax.add_patch(Ellipse((cx3-110,y0+173),16,13,fc='#C94C4C',ec='#7A2626',lw=1*PT*1.4,zorder=5))
ax.text(cx3+8,y0+172,"ME11 (FDR=0.0008)",ha='center',va='center',fontsize=22*PT,color=INK)
text(cx3,y0+192,["Co-expression network (32 modules)"])
ym=(y0+y1)/2; hline(b1[2]+8,b2[0]-6,ym); hline(b2[2]+8,b3[0]-6,ym)
r0,r1=380,680; pw=440; gap=50; xs=[(W-3*pw-2*gap)/2+i*(pw+gap) for i in range(3)]; cent=[x+pw/2 for x in xs]
labels=[["Not centrally clustered (p=0.37)"],["ADHD (n=76, RNA-seq)","no dysregulation (p=0.87)"],["mRNA splicing (top-ranked pathways;","p = 0.0003, FDR = 0.085)"]]
for x,l in zip(xs,labels): box(x,r0,x+pw,r1); text(x+pw/2,r0+225,l)
ppi(cent[0],r0+118,140); adhd_panel(cent[1],r0+112,200); splice(cent[2],r0+112,200)
bus1=r0-45; vline(cx3,y1,bus1,False); ax.plot([cent[0],cent[2]],[bus1,bus1],color=ARW,lw=3*PT,zorder=3)
for c in cent: vline(c,bus1,r0-4)
c0,c1=820,935; cw=660; cx=W/2
box(cx-cw/2,c0,cx+cw/2,c1,fc='#EEF4FB',ec='#8CA5CD')
text(cx,c0+33,["ME11 = ASD-associated module, not shared with ADHD","mRNA splicing is a lead (not FDR-significant)"],size=24,bold=True)
bus2=r1+70
for c in cent: vline(c,r1,bus2,False)
ax.plot([cent[0],cent[2]],[bus2,bus2],color=ARW,lw=3*PT,zorder=3); vline(cx,bus2,c0-4)
fig.savefig('figures/graphical_abstract.png',dpi=600,facecolor='white')
fig.savefig('figures/graphical_abstract.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
fig.savefig('figures/graphical_abstract_preview.png',dpi=254,facecolor='white')
