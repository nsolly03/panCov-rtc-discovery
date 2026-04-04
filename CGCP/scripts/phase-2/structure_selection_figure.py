#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np, os

plt.rcParams.update({
    'font.family':'DejaVu Sans','font.size':11,
    'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':False,'axes.linewidth':1.2,
    'xtick.direction':'out','ytick.direction':'out',
    'xtick.major.width':1.2,'ytick.major.width':1.2,
    'xtick.major.size':5,'ytick.major.size':5,
    'figure.facecolor':'white','axes.facecolor':'white',
    'legend.framealpha':1.0,'legend.edgecolor':'#CCCCCC',
    'pdf.fonttype':42,'ps.fonttype':42,
})

COL={'sel_edge':'#1A5276','sel_fill':'#D6EAF8',
     'alt_edge':'#717D7E','alt_fill':'#EAECEE',
     'fail_edge':'#922B21','fail_fill':'#FADBD8',
     'af3_edge':'#B7770D','af3_fill':'#FEF9E7',
     'rfree':'#CB4335','annot':'#2E4057'}

MK={'X-ray':'D','Cryo-EM':'^','AlphaFold3':'*'}

DATA={
 'NSP12-NSP7':{'selected':'7C2K','note':'Selected: best ternary complex\n(NSP12+NSP7+NSP8)','structures':[
    {'id':'7BV2','res':2.50,'method':'Cryo-EM','rfree':None,'selected':False},
    {'id':'7C2K','res':2.93,'method':'Cryo-EM','rfree':None,'selected':True},
    {'id':'6NUR','res':3.10,'method':'Cryo-EM','rfree':None,'selected':False}]},
 'NSP12-NSP8':{'selected':'7C2K','note':'Selected: two NSP8 chains\nfully resolved','structures':[
    {'id':'7BV2','res':2.50,'method':'Cryo-EM','rfree':None,'selected':False},
    {'id':'7C2K','res':2.93,'method':'Cryo-EM','rfree':None,'selected':True},
    {'id':'6NUR','res':3.10,'method':'Cryo-EM','rfree':None,'selected':False}]},
 'NSP9-NSP12':{'selected':'8SQK','note':'Only structure capturing\nNiRAN-NSP9 interface','structures':[
    {'id':'8SQK','res':3.01,'method':'Cryo-EM','rfree':None,'selected':True}]},
 'NSP10-NSP16':{'selected':'6W4H','note':'Selected: highest resolution\n+ lowest R-free','structures':[
    {'id':'6W4H','res':1.80,'method':'X-ray','rfree':0.070,'selected':True},
    {'id':'6WKQ','res':1.98,'method':'X-ray','rfree':0.097,'selected':False},
    {'id':'6WVN','res':2.00,'method':'X-ray','rfree':0.097,'selected':False}]},
 'NSP7-NSP8':{'selected':'AF3 ModeB','note':'AF3 ModeB: PHE92=6.78 A\nCrystal ModeA: PHE92=35 A (fail)','structures':[
    {'id':'7BV2\n(ModeA)','res':2.50,'method':'Cryo-EM','rfree':None,'selected':False,'failed':True},
    {'id':'6NUR\n(ModeA)','res':3.10,'method':'Cryo-EM','rfree':None,'selected':False,'failed':True},
    {'id':'AF3\nModeB','res':None,'method':'AlphaFold3','rfree':None,'selected':True}]},
 'NSP10-NSP14':{'selected':'7DIY','note':'Only available structure\nfor exonuclease complex','structures':[
    {'id':'7DIY','res':2.69,'method':'X-ray','rfree':0.264,'selected':True}]},
 'NSP13-Helicase':{'selected':'7NIO','note':'Selected: highest resolution\nX-ray homodimer','structures':[
    {'id':'7NIO','res':2.20,'method':'X-ray','rfree':0.287,'selected':True},
    {'id':'7CXM','res':2.90,'method':'Cryo-EM','rfree':None,'selected':False}]},
 'NSP12-NSP13':{'selected':'7RDY','note':'Selected: captures\npolymerase-helicase junction','structures':[
    {'id':'7CXM','res':2.90,'method':'Cryo-EM','rfree':None,'selected':False},
    {'id':'7RDY','res':3.10,'method':'Cryo-EM','rfree':None,'selected':True},
    {'id':'6XEZ','res':3.50,'method':'Cryo-EM','rfree':None,'selected':False}]},
 'NSP15':{'selected':'9HH5','note':'Selected: lowest R-free\namong all available structures','structures':[
    {'id':'6WLC','res':1.82,'method':'X-ray','rfree':0.195,'selected':False},
    {'id':'6WXC','res':1.85,'method':'X-ray','rfree':0.194,'selected':False},
    {'id':'6W01','res':1.90,'method':'X-ray','rfree':0.185,'selected':False},
    {'id':'9HH5','res':2.08,'method':'X-ray','rfree':0.023,'selected':True},
    {'id':'6VWW','res':2.20,'method':'X-ray','rfree':0.178,'selected':False},
    {'id':'2H85','res':2.60,'method':'X-ray','rfree':0.219,'selected':False},
    {'id':'7TJ2','res':3.20,'method':'Cryo-EM','rfree':None,'selected':False}]},
}

ORDER=['NSP12-NSP7','NSP12-NSP8','NSP9-NSP12',
       'NSP10-NSP16','NSP7-NSP8','NSP10-NSP14',
       'NSP13-Helicase','NSP12-NSP13','NSP15']

def get_col(s):
    if s.get('failed'):    return COL['fail_edge'],COL['fail_fill']
    if s.get('selected'):
        if s['method']=='AlphaFold3': return COL['af3_edge'],COL['af3_fill']
        return COL['sel_edge'],COL['sel_fill']
    return COL['alt_edge'],COL['alt_fill']

def draw_panel(ax, iface, data):
    structs=data['structures']
    n=len(structs)
    xs=np.arange(n)
    bw=0.55
    for i,s in enumerate(structs):
        ec,fc=get_col(s)
        lw=1.8 if s.get('selected') else 1.0
        res=s['res']
        if res is not None:
            ax.bar(i,res,width=bw,facecolor=fc,edgecolor=ec,
                   linewidth=lw,zorder=3,clip_on=False)
            ax.text(i,res+0.06,f"{res:.2f}",ha='center',va='bottom',
                    fontsize=8,color=ec,
                    fontweight='bold' if s.get('selected') else 'normal',
                    clip_on=False)
        else:
            ax.bar(i,3.6,width=bw,facecolor=COL['af3_fill'],
                   edgecolor=COL['af3_edge'],linewidth=1.8,
                   hatch='///',alpha=0.75,zorder=3,clip_on=False)
            ax.text(i,1.8,'Predicted\nmodel',ha='center',va='center',
                    fontsize=7.5,color=COL['af3_edge'],
                    fontweight='bold',clip_on=False)
        # method marker
        mk_y=(res if res else 3.6)+0.38
        mk=MK.get(s['method'],'o')
        mc=(COL['sel_edge'] if s.get('selected') and not s.get('failed')
            else COL['fail_edge'] if s.get('failed')
            else COL['af3_edge'] if s['method']=='AlphaFold3'
            else COL['alt_edge'])
        ms=10 if mk=='*' else 6
        ax.plot(i,mk_y,marker=mk,color=mc,markersize=ms,
                zorder=5,clip_on=False,
                markeredgewidth=0.8,markeredgecolor='white')
        if s.get('selected') and s['method']!='AlphaFold3':
            ax.text(i,(res or 0)+0.58,'Selected',
                    ha='center',va='bottom',fontsize=7,
                    color=COL['sel_edge'],fontweight='bold',clip_on=False)
        if s.get('failed'):
            ax.text(i,(res or 0)*0.5,'X',ha='center',va='center',
                    fontsize=16,color=COL['fail_edge'],
                    fontweight='bold',alpha=0.55,clip_on=False)
    # R-free twin
    rf=[(i,s['rfree']) for i,s in enumerate(structs) if s.get('rfree') is not None]
    if rf:
        ax2=ax.twinx()
        ax2.spines['right'].set_visible(True)
        ax2.spines['right'].set_linewidth(0.9)
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.set_ylim(0,0.40)
        ax2.set_yticks([0.0,0.1,0.2,0.3])
        ax2.tick_params(axis='y',labelsize=7,colors=COL['rfree'],
                        direction='out',length=4,width=0.9)
        ax2.set_ylabel('R-free',fontsize=8,color=COL['rfree'],
                       fontweight='bold',labelpad=2)
        xi=[v[0] for v in rf]; yi=[v[1] for v in rf]
        ax2.scatter(xi,yi,marker='D',color=COL['rfree'],s=45,
                    zorder=6,edgecolors='white',linewidths=0.8)
        for x,y in rf:
            ax2.text(x+0.18,y+0.008,f"{y:.3f}",fontsize=6.5,
                     color=COL['rfree'],va='bottom')
    ax.set_xticks(xs)
    ax.set_xticklabels([s['id'] for s in structs],fontsize=8.5,ha='center')
    ax.tick_params(axis='x',length=0,pad=6)
    ax.tick_params(axis='y',direction='out',length=5,width=1.2,labelsize=8)
    ax.set_ylim(0,4.7)
    ax.set_yticks([1.0,1.5,2.0,2.5,3.0,3.5,4.0])
    ax.set_ylabel('Resolution (A)',fontsize=9,fontweight='bold')
    ax.spines['left'].set_position(('outward',4))
    ax.spines['bottom'].set_position(('outward',4))
    ax.set_xlim(-0.5,max(n-0.5,0.5))
    ax.annotate('',xy=(-.45,1.2),xytext=(-.45,2.6),
                arrowprops=dict(arrowstyle='->',color='#AAAAAA',lw=0.9),
                annotation_clip=False)
    ax.text(-.45,1.9,'Better',fontsize=6,color='#AAAAAA',
            ha='center',va='center',rotation=90,clip_on=False)
    ax.set_title(iface,fontsize=11,fontweight='bold',pad=10)
    ax.text(0.98,0.97,data['note'],transform=ax.transAxes,
            fontsize=7,ha='right',va='top',color=COL['annot'],
            bbox=dict(boxstyle='round,pad=0.35',facecolor='#F4F6F7',
                      edgecolor='#CACFD2',linewidth=0.7))

fig,axes=plt.subplots(3,3,figsize=(17,14))
fig.subplots_adjust(hspace=0.65,wspace=0.58,
                    left=0.06,right=0.96,top=0.90,bottom=0.07)

for ax,iface in zip(axes.flat,ORDER):
    draw_panel(ax,iface,DATA[iface])

handles=[
    mpatches.Patch(facecolor=COL['sel_fill'],edgecolor=COL['sel_edge'],
                   linewidth=1.5,label='Selected structure'),
    mpatches.Patch(facecolor=COL['alt_fill'],edgecolor=COL['alt_edge'],
                   linewidth=1.0,label='Alternative (not selected)'),
    mpatches.Patch(facecolor=COL['fail_fill'],edgecolor=COL['fail_edge'],
                   linewidth=1.0,label='Failed interface criterion'),
    mpatches.Patch(facecolor=COL['af3_fill'],edgecolor=COL['af3_edge'],
                   linewidth=1.5,hatch='///',label='AlphaFold3 (predicted)'),
    mlines.Line2D([],[],marker='D',color=COL['rfree'],linestyle='None',
                  markersize=7,label='R-free value (X-ray only)'),
    mlines.Line2D([],[],marker='D',color=COL['alt_edge'],linestyle='None',
                  markersize=6,label='X-ray crystallography'),
    mlines.Line2D([],[],marker='^',color=COL['alt_edge'],linestyle='None',
                  markersize=6,label='Cryo-EM'),
    mlines.Line2D([],[],marker='*',color=COL['af3_edge'],linestyle='None',
                  markersize=9,label='AlphaFold3'),
]

fig.legend(handles=handles,loc='upper center',ncol=4,fontsize=9,
           frameon=True,facecolor='white',edgecolor='#CCCCCC',
           bbox_to_anchor=(0.5,0.965),
           columnspacing=1.5,handlelength=1.6,handletextpad=0.5)

fig.suptitle('Structural Template Selection for CGCP Pharmacophore Pipeline',
             fontsize=14,fontweight='bold',y=1.00)

OUT=os.path.expanduser(
    '~/projects/rtc-pan-coronavirus/CGCP/02-deep-dive/structure_selection')
os.makedirs(OUT,exist_ok=True)

for ext in ['png','pdf']:
    p=os.path.join(OUT,f'Fig_StructureSelection_AllInterfaces.{ext}')
    fig.savefig(p,dpi=300,bbox_inches='tight',facecolor='white')
    print(f"Saved: {p}")
plt.close()
print("Done.")
