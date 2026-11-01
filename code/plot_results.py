#!/usr/bin/env python3
"""Publication figures from saved numerical outputs; does not refit models."""
from pathlib import Path
import io
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import locked_core as core

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
COL={'sat':'#245F91','recent':'#93ADBF','linear':'#E1AA60','transfer':'#C76B47','obs':'#26333D'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
    'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.linewidth':0.8,'xtick.direction':'out','ytick.direction':'out',
    'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})

def save(fig,name):
    for ext in ['png','pdf','svg']:
        buffer=io.BytesIO()
        fig.savefig(buffer,format=ext,dpi=300,bbox_inches='tight',pad_inches=.12)
        payload=buffer.getvalue()
        if not payload:
            raise RuntimeError(f'Empty figure export: {name}.{ext}')
        (OUT/f'{name}.{ext}').write_bytes(payload)
    plt.close(fig)

def style(ax):
    ax.grid(axis='y',color='#E5E9EC',lw=.7,zorder=0)
    ax.set_axisbelow(True)

def main():
    summary=pd.read_csv(ROOT/'results/multicell_summary.csv').set_index('cell')
    new=summary.loc[[c for c in summary.index if c not in ['CS2-21','CS2-8']]]
    fig,axes=plt.subplots(1,2,figsize=(14,4.7),gridspec_kw={'wspace':.28})
    x=np.arange(len(new));w=.22
    ax=axes[0]
    for j,(column,label,color) in enumerate([
        ('recent_state_capacity_MAE_mAh','Recent-state',COL['recent']),
        ('linear_capacity_MAE_mAh','Linear',COL['linear']),
        ('saturating_capacity_MAE_mAh','Saturating',COL['sat'])]):
        vals=new[column].values
        bars=ax.bar(x+(j-1)*w,vals,w,color=color,label=label,zorder=3)
        ax.bar_label(bars,labels=[f'{v:.1f}' for v in vals],padding=3,fontsize=8,zorder=6)
    ax.set_xticks(x,[f'{c}\n'+('0.5C' if c in ['CS2-33','CS2-34'] else '1C') for c in new.index])
    ax.set_ylim(0,new[['recent_state_capacity_MAE_mAh','linear_capacity_MAE_mAh','saturating_capacity_MAE_mAh']].max().max()*1.24);ax.set_ylabel('Capacity MAE (mAh)')
    ax.set_title('(a) Frozen capacity prediction | cycles 61–101',loc='left',pad=14)
    ax.legend(frameon=False,ncol=3,loc='upper left',fontsize=8.5)
    ax.axvline(1.5,color='#BCC6CC',lw=.8,ls=':');style(ax)
    ax=axes[1]
    cells=summary.index.tolist();xx=np.arange(len(cells))
    mean=summary['cycle_mean_voltage_MAE_mV'].values
    low=summary['cycle_mean_MAE_CI95_low_mV'].values;high=summary['cycle_mean_MAE_CI95_high_mV'].values
    ax.errorbar(xx,mean,yerr=[mean-low,high-mean],fmt='none',ecolor=COL['sat'],capsize=4,lw=1.4,zorder=3)
    ax.scatter(xx,mean,s=40,c=[COL['recent']]*2+[COL['sat']]*len(new),zorder=4)
    for j in range(len(cells)):
        ax.text(j,high[j]+1.1,f'{mean[j]:.2f}',ha='center',va='bottom',fontsize=8.5)
    ax.set_xticks(xx,[f'{c}\nn={int(summary.loc[c,"confirmation_cycles"])}' for c in cells],fontsize=8.7)
    ax.set_ylabel('Mean cycle voltage MAE (mV)');ax.set_ylim(0,max(42,float(np.max(high))*1.16))
    ax.axvline(1.5,color='#BCC6CC',lw=.8,ls=':');style(ax)
    ax.set_title('(b) Cell-specific voltage | cycles 81–101',loc='left',pad=14)
    fig.subplots_adjust(left=.065,right=.985,bottom=.21,top=.87)
    fig.text(.5,.027,'Voltage whiskers: cycle-mean 95% moving-block bootstrap CI (20,000 draws, block length 5).',ha='center',fontsize=9,color='#51616D')
    save(fig,'Fig01_multicell_validation')

    targets=summary.loc[[c for c in summary.index if c!='CS2-21']]
    fig,ax=plt.subplots(figsize=(11.4,4.6));x=np.arange(len(targets));w=.29
    for offset,key,label,color in [(-w/2,'pooled_voltage_MAE_mV','Cell-specific voltage map',COL['sat']),
                                (w/2,'CS2_21_transfer_MAE_mV','CS2-21 voltage map transfer',COL['transfer'])]:
        vals=targets[key].values;bars=ax.bar(x+offset,vals,w,color=color,zorder=3,label=label)
        ax.bar_label(bars,labels=[f'{v:.2f}' for v in vals],padding=4,fontsize=9)
    ax.set_xticks(x,[f'{c}\n'+('0.5C' if c in ['CS2-8','CS2-33','CS2-34'] else '1C') for c in targets.index])
    ax.set_ylabel('Pooled common-support voltage MAE (mV)')
    ax.set_ylim(0,max(320,targets.CS2_21_transfer_MAE_mV.max()*1.25));style(ax)
    ax.legend(frameon=False,loc='upper right',ncol=2,fontsize=9)
    ax.set_title('Cell-specific calibration and voltage-map transfer | cycles 81–101',loc='left',pad=15)
    fig.subplots_adjust(left=.10,right=.99,bottom=.24,top=.85)
    fig.text(.10,.05,'Both settings use each target cell’s own frozen capacity model; this is not zero-shot transfer of the whole model.',fontsize=8.7,color='#51616D')
    save(fig,'Fig02_voltage_map_transfer')

    fig,axes=plt.subplots(2,3,figsize=(12,7.0),sharex=True)
    for ax,cell in zip(axes.ravel(),new.index):
        folder=ROOT/'results'/cell
        observed=pd.read_csv(folder/'observed_capacity.csv')
        test=pd.read_csv(folder/'capacity_predictions.csv')
        ax.axvspan(0,60,color='#EDF1F4',zorder=0)
        ax.axvline(60,color='#97A7B1',lw=.8,ls='--')
        ax.scatter(observed.cycle,observed.capacity_mAh,s=11,color=COL['obs'],alpha=.72,zorder=4)
        for model,color in [('recent_state',COL['recent']),('linear',COL['linear']),('saturating',COL['sat'])]:
            ax.plot(test.cycle,test[f'{model}_mAh'],lw=1.7,color=color,zorder=3)
        ax.set_title(cell,loc='left');ax.set_xlim(0,103);style(ax)
        ax.set_ylabel('Capacity to 3.0 V (mAh)')
    for ax in axes[1]:ax.set_xlabel('Chronological cycle index')
    handles=[Line2D([],[],marker='o',ls='',color=COL['obs'],ms=4,label='Measured'),
             Line2D([],[],color=COL['recent'],lw=2,label='Recent-state'),
             Line2D([],[],color=COL['linear'],lw=2,label='Linear'),
             Line2D([],[],color=COL['sat'],lw=2,label='Saturating'),
             Patch(color='#EDF1F4',label='Development: cycles 1–60')]
    fig.legend(handles=handles,loc='upper center',ncol=5,frameon=False,bbox_to_anchor=(.5,1.005),fontsize=9)
    fig.subplots_adjust(top=.91,bottom=.12,left=.08,right=.985,hspace=.31,wspace=.30)
    fig.text(.5,.02,'Additional-cell capacity diagnostics use the same 3.0-V analysis endpoint and the fixed chronological horizon.',ha='center',fontsize=8.6,color='#51616D')
    save(fig,'Fig03_capacity_trajectories')

    fig,axes=plt.subplots(2,3,figsize=(12.4,7.0),sharey=True)
    for ax,cell in zip(axes.ravel(),new.index):
        d=pd.read_csv(ROOT/'results'/cell/'voltage_predictions.csv.gz')
        d=d[d.window=='confirmation_81_101']
        candidates=sorted(d.cycle.unique());representative=min(candidates,key=lambda c:abs(c-91))
        a=d[(d.cycle==representative)&(d.setting=='cell_specific')]
        b=d[(d.cycle==representative)&(d.setting=='CS2_21_voltage_map_transfer')]
        ax.plot(a.time_s/60,a.measured_voltage_V,color=COL['obs'],lw=1.7,label='Measured')
        ax.plot(a.time_s/60,a.predicted_voltage_V,color=COL['sat'],lw=1.5,ls='--',label='Cell-specific')
        ax.plot(b.time_s/60,b.predicted_voltage_V,color=COL['transfer'],lw=1.3,ls='-.',label='Transferred map')
        ax.set_title(f'{cell} | cycle {representative}',loc='left')
        ax.set_xlabel('Discharge time (min)');ax.set_ylabel('Loaded terminal voltage (V)');style(ax)
        ax.set_ylim(2.94,4.25)
    axes[0,0].legend(frameon=False,fontsize=8.5)
    fig.subplots_adjust(left=.08,right=.985,top=.94,bottom=.12,hspace=.32,wspace=.19)
    fig.text(.5,.025,'Shown cycle: nearest accepted cycle to 91, fixed by chronology. Curves show the full measured trajectory to 3.0 V.',ha='center',fontsize=8.8,color='#51616D')
    save(fig,'Fig04_voltage_trajectories')

if __name__=='__main__':main()
