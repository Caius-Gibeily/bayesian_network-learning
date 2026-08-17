# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 20:03:31 2026

@author: CGIBEIL
"""
import os
import pickle as pc
import numpy as np
import pandas as pd
import seaborn.objects as so
import matplotlib.pyplot as plt
plt.rcParams['lines.linewidth'] = 3    
plt.rcParams['font.size'] = 25
# Function scripts developed for this work
os.chdir("6_bayesian_networks/")
import BN_functions as bn
import BN_simulation as sim
os.chdir("..")
datadir = "data\\"
clips = [335,340,348, 351, 352]
dat = bn.collate_data(datadir = f"{datadir}clusterArray",
                  dir_adjBoot = f"{datadir}adjBoot_AIC_",
                  clips = clips, thresh = 0.54)


for c,clip in enumerate(clips):
    for j,dx in enumerate(["TD","ASD"]):
        dat_loop = dat[(dat["dx"]==dx) & (dat["clip"]==clip)]
        
        ratioarrays = bn.calc_state_changes(dat_loop["data"].iloc[0],
                                 dat_loop["model"].iloc[0],
                                 clip, resample="bootstrap",nboots = 500, vsamp = 1e3)
        bn.save_pickle(f"deltaP_{clip}_{dx}",ratioarrays)


ratioarrays = pd.DataFrame(ratioarrays)



def format_df(df,epochs):
    df_long = df.melt(
        id_vars=["p", "p_node", "c_node", "state"],
        var_name="child_state",
        value_name="delta_p"
    ).dropna()

    df_long["child_state"] = df_long["child_state"].astype(int)


    df_range = (
        df_long
        .groupby(["p_node", "c_node", "state", "child_state"])
        .agg(dmin=("delta_p", lambda x: x.quantile(0.01)),
             dmax=("delta_p", lambda x: x.quantile(0.99)),
             dmedian=("delta_p", lambda x: x.quantile(0.5)))
        .reset_index()
    )
    df_range["plot"] = (0 > df_range["dmin"]) & (0 < df_range["dmax"])
    
    df_range["y_label"] = df_range["p_node"] + "=" + df_range["state"].astype(str)
    y_levels = list(dict.fromkeys(df_range["y_label"]))
    y_map = {k: i for i, k in enumerate(y_levels)}
    df_range["y"] = df_range["y_label"].map(y_map)

    # columns = child node
    x_levels = epochs[1:]
    x_map = {k: i for i, k in enumerate(x_levels)}
    df_range["x"] = df_range["c_node"].map(x_map)
    return df_range,x_levels,y_levels

import matplotlib.pyplot as plt

#%%

def plot_deltaPs(df_range,summary,x_levels,y_levels,
                 sfactor=6,tile_halfwidth=0.8,scale = 1):
    

    
    fig, ax = plt.subplots(figsize=(10,20))
    
    # color per child state
    colors = {
        0: "red",
        1: "blue",
        2: "green",
        3: "orange"
    }
    
    # vertical offsets inside each tile
    offsets = {
        0: 0.25,
        1: 0.08,
        2: -0.08,
        3: -0.25
    }
    
    
    priorset = 0
    for i, row in df_range.iterrows():
        if row["plot"]: continue
        currset = row["p_node"] + row["c_node"] + str(row["state"])
        
        x = row["x"]
        y = row["y"]
    
        # draw tile
        if (priorset != currset):
            if summary[(summary["p_node"]==row["p_node"]) & (summary["c_node"]==row["c_node"]) & (summary["state"]==row["state"])]["comp"].item():
                rect = plt.Rectangle(
                    (x - 0.45, y - 0.45),
                    0.9, 0.9,
                    facecolor="white",
                    edgecolor="green",
                    linewidth=4  
                )
            else:
                rect = plt.Rectangle(
                    (x - 0.45, y - 0.45),
                    0.9, 0.9,
                    facecolor="white",
                    edgecolor="grey",
                    linewidth=4 
                )
            ax.add_patch(rect)
    
        ax.plot(
            [x, x],
            [y - 0.45, y + 0.45],
            color="black",
            linewidth=2.8,
            alpha=0.6
        )
        
        ax.plot(
            [x  + (0.5 / scale) * tile_halfwidth, 
             x  + (0.5 / scale) * tile_halfwidth],
            [y - 0.45, y + 0.45],
            color="black",
            linewidth=2.8,
            alpha=0.6,
            linestyle = "dotted"
        )
    
        ax.plot(
            [x  + (-0.5 / scale) * tile_halfwidth, 
             x  + (-0.5 / scale) * tile_halfwidth],
            [y - 0.45, y + 0.45],
            color="black",
            linewidth=2.8,
            alpha=0.6,
            linestyle = "dotted"
        )
        
    
    
        # normalize ΔP into tile coordinates
        xmin = x + (row["dmin"] / scale) * tile_halfwidth
        xmax = x + (row["dmax"] / scale) * tile_halfwidth
        xmed = x + (row["dmedian"] / scale) * tile_halfwidth
        y_offset = y + offsets[row["child_state"]]
    
        # draw horizontal range
        ax.plot(
            [xmin,xmax],
            [y_offset,y_offset],
            color=colors[row["child_state"]],
            linewidth=6
        )
        ax.scatter(xmed,
                   y_offset,
                   color = colors[row["child_state"]],s=100,zorder=5
                   )
        priorset = currset
    
    # axis formatting
    ax.set_xticks(range(0,len(x_levels)))
    ax.set_xticklabels(x_levels)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    ax.tick_params(axis='x', bottom=False, labelbottom=False,
                   top=True, labeltop=True)
    
    ax.set_yticks(range(0,len(y_levels)))
    ax.set_yticklabels(y_levels)
    
    ax.set_xlabel("Child node")
    ax.set_ylabel("Parent state")
    
    ax.invert_yaxis()
    
    for spine in ax.spines.values():
        spine.set_visible(False)
    #plt.tight_layout()

    




def prepare_deltaP_df(epochs,clip):
    
    ratioarrays_TD = bn.load_pickle(f"deltaP_{clip}_TD")
    ratioarrays_ASD = bn.load_pickle(f"deltaP_{clip}_ASD")
    df_TD = pd.DataFrame(ratioarrays_TD)
    df_ASD = pd.DataFrame(ratioarrays_ASD)
    
    df_range_TD,x_levels,y_levels = format_df(df_TD,epochs)
    df_range_ASD,_,_ = format_df(df_ASD,epochs)
    
    
    group_cols = ["p_node", "c_node", "state"]
    summary_TD = (
        df_range_TD
        .loc[df_range_TD.groupby(group_cols)["dmedian"].idxmax()]
        [group_cols + ["child_state", "plot"]]
        .rename(columns={
            "child_state": "max_state_TD",
            "plot": "TD_crosses_zero"
        })
    )
    
    summary_ASD = (
        df_range_ASD
        .loc[df_range_ASD.groupby(group_cols)["dmedian"].idxmax()]
        [group_cols + ["child_state", "plot"]]
        .rename(columns={
            "child_state": "max_state_ASD",
            "plot": "ASD_crosses_zero"
        })
    )
    summary = summary_TD.merge(summary_ASD, on=group_cols, how="inner")
    summary["shared"] = (
        summary["max_state_TD"] == summary["max_state_ASD"]
    )
    
    summary["both"] = (
        ~summary["TD_crosses_zero"] &
        ~summary["ASD_crosses_zero"]
    )
    
    summary["comp"] = summary["shared"] & summary["both"]

    return df_range_TD,df_range_ASD,summary_TD,summary,x_levels,y_levels

for c,clip in enumerate(clips):
    
    epochs = dat["data"].iloc[c*2].columns.tolist()
    
    df_range_TD,df_range_ASD,summary_TD,summary,x_levels,y_levels = prepare_deltaP_df(epochs,clip)
    for dx in ["TD","ASD"]:
        plot_deltaPs(eval(f"df_range_{dx}"), summary, x_levels, y_levels,tile_halfwidth=0.6)
        plt.savefig(f"deltaP_{clip}_{dx}.svg")