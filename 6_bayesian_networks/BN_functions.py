# -*- coding: utf-8 -*-
"""
@author: CGIBEIL

Bayesian Network Functions - to complement 
"""


import numpy as np
import pandas as pd
import pickle as pc
import string
import matplotlib.pyplot as plt
from pgmpy.models import DiscreteBayesianNetwork
import itertools
# Epoch labels
epochs = list(string.ascii_uppercase)
#%% Pickle functions 
def save_pickle(savedir,var):
    """
    Parameters
    ----------
    savedir : base directory to save pickle (.pkl)
    var : object to save

    Returns
    -------
    None.
    """
    
    with open(f"{savedir}.pkl", "wb") as f:
        pc.dump(var, f, protocol=pc.HIGHEST_PROTOCOL)
        
        
def load_pickle(datadir):
    """
    Parameters
    ----------
    datadir : directory + filename to load from

    Returns
    -------
    var : variable to assign unpickled file to 

    """
    with open(f"{datadir}.pkl","rb") as file:
        var = pc.load(file)
    return var
        
#%% Translation functions
def adj2DAG(adj):
    """

    Parameters
    ----------
    adj : input adjacency table (1 = edge, 
                                 0 = no edge)

    Returns
    -------
    dag : converted directed acyclic graph (dag) object

    """
    from pgmpy.base import DAG
    import string
    epochs = list(string.ascii_uppercase)[0:adj.shape[0]]

    edges = [(parent,child) for p,parent in enumerate(epochs) 
             for c,child in enumerate(epochs) if adj[p, c] == 1]
    dag = DAG()
    dag.add_nodes_from(epochs)  
    dag.add_edges_from(edges)
    return dag
    
def DAG2adj(model):    
    """
    Function: 
        Performs the inverse of the adj2DAG, taking a DAG object and converting it 
        into a DAG
    Parameters
    ----------
    model : network model (pgmpy.base.DAG or 
                           pgmpy.models.DiscreteBayesianNetwork)
    Returns
    -------
    adj : 

    """
    adj = np.zeros([len(model.nodes),len(model.nodes)])
    epochs = pd.Index(model.nodes)
    
    # Where 1 indicates an edge between two variables 
    for i,epoch in enumerate(model.edges()):
      
        adj[(epochs.get_loc(epoch[0]),
             epochs.get_loc(epoch[1]))] = 1                            
    return adj


#%% Bootstrap functions
def get_bootstrapped_samples(data,nboot=1e3,thresh=0.54):
    """
    Parameters
    ----------
    data : input state-matrix data (Pandas object)
    nboot : optional, number of bootstraps to run
    thresh : edge frequency threshold
    Returns
    -------
    adjAll : thresholded matrices
    adjs : bootstrapped adjacency matrices

    """
    nparts = data.shape[0]
    adjs = {}
    adjAll = np.zeros([data.shape[1],data.shape[1]])
    for i in range(nboot):
        bootID = np.random.choice(nparts,nparts,replace=True)
        sampleboot = data.iloc[bootID,:]
        
        _,adjs[i],_ = learn_structure(sampleboot,"AICScore")
    
    if thresh is None:
        return adjs
    else:
        
        adjAll = threshold_bootstrapped(adjs, thresh)
        return adjAll,adjs

def bootstrap_clips(rootdir,clips,nboot,score="AICScore"):
   """
    Parameters
    ----------
    datdir : directory and root file name

    without the  
    
    clips : list of movie clips
    
    nboot : number of bootstraps to run
    
    thresh : the threshold to apply on frequencies of edges across bootstrapped 
    networks
    
    score : Scoring method to use for structure learning. The default is "AICScore".

    Returns
    -------
    adjBootTD : All n bootstrapped adjacency matrices 
    for TD group (returned as dict)
    
    adjBootASD : All n bootstrapped adjacency matrices 
    for ASD group (returned as dict)

    """
   adjBootTD = {}
   adjBootASD = {}
   for i,m in enumerate(clips):
        adjBootTD[m] = {}
        adjBootASD[m] = {}
        
        arrayTD =  pd.read_csv(f"{rootdir}TD_clip{str(m)}.csv",
                               header=None)
        
        adjBootTD[m] = get_bootstrapped_samples(arrayTD,
                                                nboot,thresh=None)
        
        arrayASD =  pd.read_csv(f"{rootdir}ASD_clip{str(m)}.csv",
                                header=None)
        adjBootASD[m] = get_bootstrapped_samples(arrayASD,
                                                 nboot,thresh=None)

   return adjBootTD, adjBootASD 



def threshold_bootstrapped(adjBoots,clips,thresh):
    """
    Parameters
    ----------
    adjBoots : input dictionary of bootstrapped adjacency matrices
    thresh : edge frequency cutoff to use

    Returns
    -------
    adjAll : thresholded bootstrapped adjacency matrices
    """
    adj_all = {}
    for c,clip in enumerate(clips):
        
        adjs = list(adjBoots[clip].values())
        adj_all[clip] = np.mean(adjs,axis=0) >= thresh
    
    return adj_all

#%% Collate data
def collate_data(datadir,dir_adjBoot,clips,thresh):
    """
    Parameters
    ----------
    datadir : data directory of 
    state-epoch matrices
    clips : 
    cutoff : TYPE
        DESCRIPTION.

    Returns
    -------
    dat : TYPE
        DESCRIPTION.

    """
    dat = []
    
    
    adjBootTD = load_pickle(f"{dir_adjBoot}TD")
    adjBootASD = load_pickle(f"{dir_adjBoot}ASD")

    adjAllTD  =  threshold_bootstrapped(adjBootTD,clips,thresh)
    adjAllASD =  threshold_bootstrapped(adjBootASD,clips,thresh)
    
    adjAll = {"TD":adjAllTD,
              "ASD":adjAllASD}
    epochs = list(string.ascii_uppercase)
    for i, clip in enumerate(clips):
        
        for dx in ["TD","ASD"]:
            data =  pd.read_csv(f"{datadir}{dx}_clip{str(clip)}.csv",header=None).astype(str)

            data.columns = epochs[0:data.shape[1]]
         
            DAG_thresh = adj2DAG(adjAll[dx][clip])
            model = learn_parameters(data, DAG_thresh)
            
            dat.append({"clip": clip, "dx": dx, "model": model, 
                        "adj": adjAll[dx][clip], "data": data })
            
         
    dat = pd.DataFrame(dat)
    return dat



#%% NMI calculation and plotting
def calc_nmi(parent,child):
    """
    Parameters
    ----------
    parent : column of states in the upstream epoch
    child :  column of states in the downstream epoch

    Returns
    -------
    The pairwise normalised mutual information (NMI) between upstream and 
    downstream epochs
    """
    
    states_p, counts_p = np.unique(parent, return_counts=True)
    p_p = counts_p / counts_p.sum()
    
    states_c, counts_c = np.unique(child, return_counts=True)
    p_c = counts_c / counts_c.sum()
    
    # Mutual Information
    mi = 0.0
    for i, p_pi in enumerate(p_p):
        for j, p_cj in enumerate(p_c):
            joint = np.mean(
                (parent == states_p[i]) &
                (child == states_c[j])
            )
            if joint > 0:
                mi += joint * np.log2(joint / (p_pi * p_cj))
    
    # Entropies
    H_p = -np.sum(p_p * np.log2(np.clip(p_p, 1e-12, 1)))
    H_c = -np.sum(p_c * np.log2(np.clip(p_c, 1e-12, 1)))

    # Normalised MI
    denom = (H_p + H_c)

    if denom == 0:
        return 0.0
    return 2 * mi / denom


def calc_mi_model(data,model):
    """
    Parameters
    ----------
    data : input state-matrix dataframe (Pandas object)
    model : model graph structure (either a pgmpy.base.DAG or 
                                   pgmpy.models.DiscreteBayesianNetwork)

    Returns
    -------
    mis : all NMIs for each pair of epochs
    """
    # iterate through each node pair (parent_node,child_node)
    mis = {}
    for p, pnode in enumerate(model.nodes):
        for c,cnode in enumerate(model.nodes):
            if p == c:
                continue
            mis[(pnode,cnode)] = calc_nmi(data[pnode],data[cnode])
    return mis

def calc_mi_dist(mi_dict):
    """
    Function:
        Find the distances between pairs of epochs and 
        associate them with the NMI between them
    Parameters
    ----------
    mi_dict : dictionary of mutual information scores per pair of variables
    adjBoot : list of adjacency matrices (size [NxN] for N variables)

    Returns
    -------
    mi : TYPE
        DESCRIPTION.
    dists : TYPE
        DESCRIPTION.


    """
    
    pairs = list(mi_dict.keys())
    epochs = pd.Index(np.unique(list(mi_dict.keys())))
    mi = np.array(list(mi_dict.values()))
    dists = np.array([epochs.get_loc(pair[0])-
                      epochs.get_loc(pair[1]) for pair in pairs])
    
    mi = mi[dists>0]
    dists = dists[dists>0]

    return mi, dists

    
def find_ci_edges(model,data):
    
    """
    Function: 
        The recovered BN structure makes a statement regarding which variables 
        are conditionally separated (d-separation criteria). If a conditional dependency 
        between two unconnected variables is explained by the network structure, the
        two variables are said to be conditionally independent. Otherwise, the two variables are said to 
        have a conditional dependency not explained by the model. 
        
    Parameters
    ----------
    model : BN model of type pgmpy.models.DiscreteBayesianNetwork
    data : state-matrix 

    Returns
    -------
    ci_array : array of implied pairwise conditional independencies between epochs.
    Where an unexplained conditional dependency arises, this is marked as True
    """
    
    from pgmpy.metrics import implied_cis
    from pgmpy.estimators.CITests import chi_square
    epochs = data.columns

    cis = implied_cis(model,data,ci_test=chi_square)
    cis = cis[cis["p-value"]<0.05] 
    ci_array = np.zeros([len(model.nodes),len(model.nodes)])

    ci_array[epochs.get_indexer(cis["u"]),
             epochs.get_indexer(cis["v"])] = 1
    return ci_array
    

def plot_mi_adjacency(model, adj,ci_array, mi_dict, 
                      bootstrapped,cmap = "Blues", fig=None,ax=None):
    """
    Function: 
        We have derived (1) bootstrapped graphs (bootstrapped), (2) an edge frequency 
        threshold to incorporate edges that recur across more than the proportion 
        of bootstrapped samples set by the threshold, (3) the normalised mutual information 
        (NMI) and implied conditional indepedendencies between each pair of variables
        
        This function then plots these pieces of information together
    Parameters
    ----------
    model : Bayesian Network model (pgmpy.models.DiscreteBayesianNetwork)
    adj : Adjacency matrix
    ci_array : Implied conditional independencies (1 = implied independence given structure violated)
    mi_dict : Collection of normalised mutual information (NMI) values calculated per pair of epochs
    bootstrapped : Bootstrapped networks
    cmap : TYPE, optional
        "Blues" for TD group, "Greens" for ASD group
    fig : fig on which to plot
    ax : figure axes

    Returns
    -------
    None.

    """
    from matplotlib import patches
    
    nodes = list(model.nodes())
    n = len(nodes)
    adj_matrix = np.zeros((n, n))

    # Fill adjacency matrix with MI values
    for (parent, child), mi in mi_dict.items():
        if parent in nodes and child in nodes:
            i, j = nodes.index(parent), nodes.index(child)
            adj_matrix[i, j] = mi
    
    bootAdj = sum(bootstrapped.values())/len(bootstrapped)

    
    if fig == None:
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.set_aspect("equal")

    # Coordinates
    x, y = np.meshgrid(range(n), range(n))
    x, y, mi_values = x.flatten(), y.flatten(), adj_matrix.flatten()

    # Keep only upper triangle (lower triangle contains forbidden edges)
    mask = np.triu(np.ones((n, n)), k=1).flatten().astype(bool)
    x, y, mi_values = x[mask], y[mask], mi_values[mask]
    
    # bootstrapped edge frequencies
    bs_values = bootAdj.flatten()
    bs_values = bs_values[mask] 

    # Normalise circle sizes
    if np.nanmax(mi_values) > 0:
        sizes = 800 * (mi_values / np.nanmax(mi_values))
    else:
        sizes = np.zeros_like(mi_values)

    # Plot circles
    sc = ax.scatter(
        x,
        y,
        s=sizes,
        c=bs_values,
        cmap=cmap,
        edgecolors="none",
        alpha=0.9
    )

    # Plot shaded areas corresponding to incorporated edges at the selected cutoff
    parents, children = adj.nonzero()
 
    # Draw light gridlines for upper triangle
    for i in range(n):
        for j in range(i, n):
            ax.plot([j - 0.5, j + 0.5], [i - 0.5, i - 0.5], color="lightgray", lw=0.8)
            ax.plot([j - 0.5, j + 0.5], [i + 0.5, i + 0.5], color="lightgray", lw=0.8)
            ax.plot([j - 0.5, j - 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.8)
            ax.plot([j + 0.5, j + 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.8)

    # Hide lower triangle with white polygon
    triangle = plt.Polygon(
        [[-0.5, n - 0.5], [-0.5, -0.5], [n - 0.5, n - 0.5]],
        closed=True,
        facecolor="white",
        edgecolor="none",
        zorder=3,
    )
    ax.add_patch(triangle)

    # Add uniform bold border enclosing upper triangle
    border = patches.Polygon(
        [[-0.5, -0.5], [n - 0.5, -0.5], [n - 0.5, n - 0.5]],
        closed=False,
        fill=False,
        edgecolor="black",
        linewidth=2.5,
        capstyle='projecting',
        joinstyle='miter',
        zorder=4,
    )
    ax.add_patch(border)

    # Add bold diagonal
    ax.plot(
        [-0.5, n - 0.5],
        [-0.5, n - 0.5],
        color="black",
        lw=2.5,
        solid_capstyle="projecting",
        zorder=4,
    )

    
    # Ticks only on top
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    
    
    ax.set_xticklabels(epochs[:n], fontsize=35)
    ax.set_yticklabels(epochs[:n],fontsize=35)
    ax.xaxis.tick_top()

    # Clean up
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.invert_yaxis()

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, shrink=0.8)
    sc.set_clim(vmin=0, vmax=1)
    cbar.ax.tick_params(labelsize=20)
    
    rgba = np.zeros((*adj.shape, 4))
    rgba[adj == 1] = [1.0, 0.75, 0.8, 1.0]
    rgba[ci_array==1] = [0.42,0.65,0.73,1]
    ax.imshow(
        rgba,
        alpha=0.5
        )
    return bs_values, mi_values
    
    
def plot_mi_adj_allclips(dat,clips,adjBootTD,adjBootASD):
    mi_dat = []
    fig, axs = plt.subplots(2, 5, figsize=(50, 15)) 
    adjBoots = {"TD":adjBootTD,
                "ASD":adjBootASD}
    col = ["Blues","Greens"]
    for i,clip in enumerate(clips):
        for j,dx in enumerate(["TD","ASD"]):
             dat_loop = dat[(dat["dx"]==dx) & (dat["clip"]==clip)]
             mi_dict =  calc_mi_model(dat_loop["data"].iloc[0],
                                      dat_loop["model"].iloc[0])
             ci_array = find_ci_edges(dat_loop["model"].iloc[0], 
                                      dat_loop["data"].iloc[0])
             
             bs_values,mi_values = plot_mi_adjacency(dat_loop["model"].iloc[0], 
                               dat_loop["adj"].iloc[0],
                               ci_array,mi_dict, 
                               adjBoots[dx][clip],
                               cmap = col[j],
                               fig=fig,ax=axs[j,i])
             
             axs[j,i].text(s="Movie " + str(i+1), x = 0, y = len(dat_loop["adj"].iloc[0])*0.6)
             
             mi, dists =  calc_mi_dist(mi_dict)

             epochs = list(string.ascii_uppercase)
             included = []
             for e, edge in enumerate(mi_dict):
                if epochs.index(edge[0]) < epochs.index(edge[1]):
                    included.append(dat_loop["adj"].iloc[0][epochs.index(edge[0]),epochs.index(edge[1])])
                
            
             mi_dat.append({"clip":clip,"dx":dx,"mis":mi, 
                            "dists":dists, "included":included, "freqs":bs_values, "mi_vals":mi_values
                 })
    return mi_dat
  
#%% Network learning
def generate_forbidden_edges(columns):
    """
    Parameters
    ----------
    columns : column names of state-matrix. Epochs are labelled alphabetically 

    Returns
    -------
    forbidden : TYPE
        DESCRIPTION.

    """
    forbidden = []
    for i in range(len(columns)):
        for j in range(i+1, len(columns)):
            forbidden.append((columns[j], columns[i]))  # No future -> past edges
    return forbidden

def learn_structure(data,score="AICScore"): 
    """
    Parameters
    ----------
    data : input data (as str if directory or as numpy array/Pandas dataframe)
    score : TYPE, The default structure learning score to use. The default is 
    Akaike Information Criterion (AIC)

    Returns
    -------
    best_model : the learned structure
    adj : recovered adjacency matrix
    data : mirror output of data

    """
    
    from pgmpy.estimators import BIC,HillClimbSearch,AIC,BDeu
    from pgmpy.estimators import ExpertKnowledge



    epochs = list(string.ascii_uppercase)
    data.columns  = epochs[0:data.shape[1]]

    # Restrict structures to only 
    black_list= generate_forbidden_edges(data.columns)
    forbidden = ExpertKnowledge(forbidden_edges=black_list)
    if score == "BICScore":
        scoring_method = BIC(data)
    elif score == "AICScore":
        scoring_method = AIC(data)
    elif score == "BDeu":
        scoring_method = BDeu(data)

    state_names={data.columns[i]:np.unique(data.iloc[:,i]) for i in range(data.shape[1])}
    hc = HillClimbSearch(data,state_names=state_names)
    best_model = hc.estimate(scoring_method=scoring_method,
                             expert_knowledge=forbidden,show_progress=True) # ,tabu_length=10
    
    adj = DAG2adj(best_model)

    return best_model, adj,data


def learn_parameters(data,dag): 
    """
    Parameters
    ----------
    data : state-epoch matrix
    dag : directed acyclic graph object

    Returns
    -------
    model : model graph and parameters (pgmpy.models.DiscreteBayesianNetwork)
    """
    
    from pgmpy.estimators import MaximumLikelihoodEstimator, ExpectationMaximization,BayesianEstimator
    
    model = DiscreteBayesianNetwork(dag)
    model.fit(data, estimator=MaximumLikelihoodEstimator)

    return model

"""
def learn_parameters(data, dag):
    from pgmpy.estimators import MaximumLikelihoodEstimator, ExpectationMaximization,BayesianEstimator
    
    model = DiscreteBayesianNetwork(dag)
    model.add_nodes_from(data.columns)
    estimator = BayesianEstimator(model, data)
  
    m = model.fit(
    data, 
    estimator=BayesianEstimator, 
    prior_type='BDeu', 
    equivalent_sample_size=10
)

    return model_adj
"""

def compute_scores(dat,thresh):
    from pgmpy.metrics import correlation_score, log_likelihood_score,fisher_c
    from pgmpy.estimators.CITests import chi_square
    score = []
    for i,d in dat.iterrows():
        corr = correlation_score(d["model"], d["data"])
        lls = log_likelihood_score(d["model"], d["data"])
        fc = fisher_c(d["model"], d["data"], ci_test=chi_square, show_progress=False)
        #ss = structure_score(d["model"].iloc[0], d["data"].iloc[0],"aic-g")
        score.append({"clip":d["clip"],
                       "dx": d["dx"],
                       "thresh":thresh,
                       "corr":corr, 
                       "lls":lls, 
                       "fc": fc})
    return score

def collate_scores(datadir,clips,threshes):
    scores = []
    for i, thresh in enumerate(threshes):
        dat = collate_data(f"{datadir}clusterArray",
                          f"{datadir}adjBoot",clips,thresh)
        scores.extend(compute_scores(dat,thresh))
    return pd.DataFrame(scores)

def plot_edge_thresh(scores,threshes,thresh):
    import seaborn as sns
    metrics = ["corr", "lls", "fc"]
    metric_labels = ["Correlation score", "Log-likelihood score", "Fisher c"]
    groups = ["TD", "ASD"]

    fig, axs = plt.subplots(
        nrows=2,
        ncols=3,
        figsize=(30, 16),
        sharex=True
    )

    for i, dx in enumerate(groups):
        scores_dx = scores[scores["dx"] == dx]

        for j, metric in enumerate(metrics):
            ax = axs[i, j]

            sns.lineplot(
                data=scores_dx,
                x="thresh",
                y=metric,
                hue="clip",
                units="clip",
                estimator=None,
                ax=ax,
                linewidth=7,
                alpha=0.8,
                legend=False
            )
            
            if j == 0:
                ax.vlines(thresh,0,1,"red",linewidth=5,linestyle="--")
            if j == 1:
                ax.vlines(thresh,-1500,-300,"red",linewidth=5,linestyle="--")
                            
            if j == 2:
                ax.hlines(0.05,thresh,threshes[-1],"red",linewidth=5,linestyle="--")
                ax.vlines(thresh,0,1,"red",linewidth=5,linestyle="--")
            if i == 0:
                ax.set_title(metric_labels[j])
            if j == 0:
                ax.set_ylabel(dx)
            
            else:
                ax.set_ylabel("")
            if j == 1: 
                ax.set_xlabel("Edge frequency threshold")
            
            #ax.set_xticks(threshes[::4])

    plt.tight_layout()


#%% Model evaluation
def js_distance(p,q):
    """
    Compute Jensen-Shannon distance between two reference probability mass functions, p and q
    Parameters
    ----------
    p : list 
        first p.m.f.
    q : list
        second p.m.f.

    Returns
    ----------
    float 
        Jensen-Shannon distance

    """
    from scipy.spatial import distance
    if np.isnan(p).any() or np.isnan(q).any():
        return np.nan
    else: 
        return distance.jensenshannon(p, q)

def compare_recovered_ppds(model1,model2, fun=np.mean):
    """
    Compares posterior probability difference between each pair of equivalent nodes across two pgmpy.DiscreteBayesianNetwork models, model1 and model2

    Parameters
    ----------
    model1 : pgmpy.DiscreteBayesianNetwork
    model2 : pgmpy.DiscreteBayesianNetwork
    fun : function
        Summary function, default = np.mean

    Returns 
    ----------
    float
        summary difference of the two models
    """
    from pgmpy.inference import CausalInference
    variables = model1.nodes()
    variables2 = model2.nodes()
    if len(variables) != len(variables2):
        raise ValueError("Unequal number of variables")
    #model1Inf = VariableElimination(model1)
    #model2Inf = VariableElimination(model2)
    model1Inf = CausalInference(model1)
    model2Inf = CausalInference(model2)
    jsds = []
    for i, node in enumerate(variables):
        parents = model2.get_parents(node)
        par_states = [model2.states[k] for k in variables if k in parents]
        combs = itertools.product(*par_states)
        
        for c, comb in enumerate(combs):

            evidence = {parents[p]: state_val for p,state_val in enumerate(comb)}
            print(evidence)
            query_fit = model2Inf.query(variables=[node],do=evidence)
            query_ground = model1Inf.query(variables=[node],do=evidence)
            diff = len(query_fit.values) - len(query_ground.values)
            if diff == 0:   
                jsds.append(js_distance(query_fit.values, query_ground.values))
            elif diff < 0:
                query_fit.values = np.append(query_fit.values,np.repeat(0,abs(diff)))
                jsds.append(js_distance(query_fit.values, query_ground.values))
            elif diff > 0:
                query_ground.values = np.append(query_ground.values,np.repeat(0,abs(diff)))
                jsds.append(js_distance(query_fit.values, query_ground.values))                
    jsd_fn = fun(jsds)
    return jsd_fn


def hamming_distance(adj1,adj2,directed=True):
    """
    Compares structural Hamming distance between two adjacency matrices, adj1 and adj2

    Parameters
    ----------
    adj1 : pgmpy.DAG
    adj2 : pgmpy.DAG
    directed : Boolean
        whether the adjacency matrices represent digraphs

    Returns 
    ----------
    float
        Hamming distance [0,1]
    """
    if not directed:
        T = np.triu(adj1 + adj1.T, 1)
        P = np.triu(adj2 + adj2.T, 1)
        diff = (T != P).astype(int)
        possible = T.size // 2
        raw = diff.sum()
    else:
        diff = (adj1 != adj2).astype(int)
        raw = diff.sum()
        possible = adj1.size
    return raw / possible

import numpy as np
import pandas as pd


def _cpd_log_prob(cpd, values_by_var, alpha=1, fallback_prior=None):
    """
    Exact lookup of P(node = observed | parents = observed) from a fitted
    pgmpy TabularCPD, using the CPD's own stored array rather than sampling.

    Parameters
    ----------
    cpd : pgmpy.factors.discrete.TabularCPD
        The fitted CPD for a single node (as stored in model.get_cpds(node)).
    values_by_var : dict[str, str]
        Observed state (as a string) for every variable in cpd.variables
        (the node itself plus its parents, in whatever order pgmpy stored
        them — we read that order directly from cpd.variables so we never
        have to assume it matches epoch/column order elsewhere).
    alpha : float
        Pseudo-count used only for the fallback path below.
    fallback_prior : dict[str, float] or None
        Optional precomputed (state -> smoothed frequency) fallback, used
        only if the observed state combination was never seen in training
        (e.g. a rare state that didn't appear for any other participant in
        this LOO fold). If None, falls back to a small floor probability.

    Returns
    -------
    prob : float
        P(node = observed | parents = observed), a real probability read
        directly from the CPD table (not an importance-weighted estimate).
    """
    try:
        idx = []
        for var in cpd.variables:
            state_list = cpd.state_names[var]
            idx.append(state_list.index(values_by_var[var]))
        prob = float(cpd.values[tuple(idx)])
        # Guard against a structurally-present-but-zero-mass entry (can
        # happen with unsmoothed ML-fit CPDs on sparse data).
        if prob <= 0 or not np.isfinite(prob):
            raise ValueError
        return prob
    except (KeyError, ValueError):
        # Observed state (or a parent's observed state) was never seen
        # for this variable during training for this LOO fold. Fall back
        # to an alpha-smoothed frequency if provided, else a small floor.
        if fallback_prior is not None:
            node = cpd.variables[0]
            return float(fallback_prior.get(values_by_var[node], 1e-16))
        return 1e-16


def predictBayes(data, adj, model, cond="evidence", alpha=1, seed=None):
    """

    Parameters
    ----------
    data : state-matrix dataframe object (rows = participants, columns =
        epochs). Index need not be 0..n-1 on entry -- it is normalised
        internally.
    adj : adjacency matrix (used to build the M0/M2 structure via
        adj2DAG; ignored for "forder", which derives its own truncated
        structure from `model`).
    model : the fitted pgmpy DiscreteBayesianNetwork for this clip/group,
        used only to recover the DAG structure for the "forder" (M1)
        condition via DAG2adj(model).
    cond : "rand" (M0, state-frequency null model), "forder" (M1,
        first-order-edges-only model), or "evidence" (M2, full model).
    alpha : Dirichlet pseudo-count for the M0 / fallback frequency estimate.
    seed : unused (kept only for call-site compatibility with
        get_all_predictions, which still passes it).

    Returns
    -------
    pred : ndarray, shape (n_participants, n_epochs)
        Per-epoch log10 predictive probability for each held-out
        participant's actually-observed state.

    Notes
    -----
    Each participant's sequence is fully observed, so the joint probability
    of that sequence factorises exactly as the product (sum, in log space)
    of each node's own CPD evaluated at its actual parents' observed
    values. There is no latent variable to integrate over.

    """
    data = data.reset_index(drop=True)
    epochs = list(data.columns)
    n = data.shape[0]

    pred = np.zeros([n, len(epochs)])

    if cond == "forder":  # M1: first-order edges only
        base_adj = DAG2adj(model)
        base_adj[~np.eye(base_adj.shape[0], dtype=bool, k=1)] = 0
        dag = adj2DAG(base_adj)
    elif cond in ("evidence", "rand"):  # M2 (full) or M0 (null)
        dag = adj2DAG(adj)
    else:
        raise ValueError(f"Unrecognised cond: {cond!r}")

    for i in range(n):
        samp = data.iloc[i]
        data_train = data.drop(i)

        epoch_fallback = {}
        for ep in epochs:
            states, counts = np.unique(data_train[ep].astype(str), return_counts=True)
            smoothed = (counts + alpha) / (np.sum(counts) + states.size * alpha) # smoothed Dirichlet prior
            epoch_fallback[ep] = dict(zip(states, smoothed))

        if cond == "rand":
            model_rec = None 
        else:
            model_rec = learn_parameters(data_train, dag)

        row_pred = np.zeros(len(epochs))
        for node in epochs:
            obs_state = str(samp[node])

            if cond == "rand": 
                prob = epoch_fallback[node].get(obs_state, 1e-16)

            else:  
                cpd = model_rec.get_cpds(node)
                parents = cpd.variables[1:] 

                values_by_var = {node: obs_state}
                for p in parents:
                    values_by_var[p] = str(samp[p])

                prob = _cpd_log_prob(
                    cpd, values_by_var, alpha=alpha,
                    fallback_prior=epoch_fallback[node],
                )

            row_pred[epochs.index(node)] = np.log10(prob + 1e-16)

        pred[i, :] = row_pred

    return pred
    
def get_all_predictions(clips,data,nreps):
    
    """
    Function: 
        Collate predictive accuracies for each query participant across all epochs
        under three different models:
            M0 - null model, where predictions are drawn from a posterior categorical 
            distribution informed by state frequencies in each epoch
            M1 - first-order only model, predictions are drawn from a reduced model 
            containing only the first-order edges of each BN
            M2 - full model, containing all edges learned via the structure learning
            and bootstrapping steps
        The posterior probabilities of entering each state are then inferred 
        for each participant.
        The function using parallel programming to accelerate the runtime of the 
        function
    Downstream: 
        Bayes factors are computed between M1 and M0 and M2 and M0, providing the
        relative evidence favouring 
        the full (M2) or partial (M1) models in predicting subsequent states over
        the null model.
        
    Parameters
    ----------
    clips : List of movie clips to perform the prediction over
    data : dataframe object collating state-matrix, adjacency matrix, 
    model parameters, dx and clip id (Pandas object)
    
    nreps : Number of likelihood weighted samples to take

    Returns
    -------
    preds: A pandas dataframe of epoch-wise prediction probabilities for each clip
    and participant who viewed that clip 

    """
    from tqdm import tqdm
    preds = []
    for i,clip in enumerate(clips):
        
        it = tqdm(range(nreps), desc="Likelihood-weighted sampling")
        for rep in it:
            seed = np.random.randint(0,1e5)
            
            # TD #################################################################
            dat_loop = data[(data["dx"]=="TD") & (data["clip"]==clip)]
            preds.append({"clip": clip, "dx":"TD", "rep":rep, "cond":"evidence", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "evidence",
                                                   seed = seed),
                          "seed": seed})
            
            preds.append({"clip": clip, "dx":"TD", "rep":rep, "cond":"rand", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "rand",
                                                   seed = seed),
                          "seed": seed})
            
            preds.append({"clip": clip, "dx":"TD", "rep":rep, "cond":"forder", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "forder",
                                                   seed = seed),
                          "seed": seed})
          
            
            # ASD ############################################################
            dat_loop = data[(data["dx"]=="ASD") & (data["clip"]==clip)]
            preds.append({"clip": clip, "dx":"ASD", "rep":rep, "cond":"evidence", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "evidence",
                                                   seed = seed),
                          "seed": seed})
            
            
            preds.append({"clip": clip, "dx":"ASD", "rep":rep, "cond":"rand", 
                          "pred_array":  predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "rand",
                                                   seed = seed),
                          "seed": seed})
            
            preds.append({"clip": clip, "dx":"ASD", "rep":rep, "cond":"forder", 
                          "pred_array": predictBayes(dat_loop["data"].values[0],
                                                   dat_loop["adj"].values[0],
                                                   dat_loop["model"].values[0],
                                                   cond = "forder",
                                                   seed = seed),
                          "seed": seed})
            
        
            

            #################################################################
    preds = pd.DataFrame(preds)
    return preds

def plot_predictions(preds):
    """
    Plot cumulative distribution function of log-Bayes factors (BFs) across TD and ASD participants
    Parameters
    ----------
    preds : pandas.DataFrame

    Returns 
    ----------
    mediandat : pandas.DataFrame
        Median log-BFs by diagnostic group and clip
    """
    import seaborn as sns
    preds["pred_array_logsum"] = preds["pred_array"].apply(np.sum,axis=1)

    wide = (
        preds
        .pivot_table(
            index=["dx", "clip"],
            columns="cond",
            values="pred_array_logsum")
        )
        
    ratios = (
        wide
        .assign(
            
            ev_rand = wide["evidence"]-wide["rand"],
            forder_rand = wide["forder"]-wide["rand"],
            ev_forder = wide["evidence"]-wide["forder"])   
        .reset_index()
    )
        
    long = ratios.melt(
        id_vars=["dx", "clip"],
        value_vars=["ev_rand", "forder_rand", "ev_forder"],
        var_name="ratio_type",
        value_name="ratio"
    )

    clips = long["clip"].unique()
    types = long["ratio_type"].unique()


    fig, axs = plt.subplots(1,5, figsize=(125,25))

    # make axs iterable for 1-row / 1-col cases
    linestyles = {"ev_rand": "-", "forder_rand": "--","ev_forder":":"}
    palette = {"TD": "blue", "ASD": "green","forder":"grey","evidence":"blue"}
    mediandat = []
    for i, clip in enumerate(clips):
        ax = axs[i]
        for j, typ in enumerate(types):
            
            
            subset = long[(long["clip"] == clip) & (long["ratio_type"] == typ)]
            subset = subset.explode(["ratio"], ignore_index=False)
            #subset["ratio"] = 10**subset["ratio"]
            #subset["ratio"] = np.log10(np.add(subset["ratio"].astype(float),1e-8))
            sns.ecdfplot(
                data=subset,
                x="ratio",
                hue="dx",
                palette = palette,
                linestyle=linestyles[typ],
                linewidth=12,
                ax=ax
            )
            ax.hlines(0.5,-2,3,color="black",linewidth=4)
            ax.vlines(0,0,1,color="black",linewidth=4)
            ax.set_xlim([-2,2])
            ax.set_xlabel("")
            ax.minorticks_on()
            ax.grid(which='minor', linestyle=':', linewidth=2, color = "grey")
            ax.grid(which='major', linestyle='-', linewidth=2, color='black') 



            subset["ratio"] = subset["ratio"].astype(float)

            gmedian = subset.groupby("dx")["ratio"].median()
            glQ = subset.groupby("dx")["ratio"].apply(np.quantile,0.25)
            ghQ = subset.groupby("dx")["ratio"].apply(np.quantile,0.75)
            
            mediandat.append({"TDmed":gmedian["TD"],
                            "TDl":(glQ["TD"]),
                            "TDh":(ghQ["TD"]),
                            
                            "ASDmed":gmedian["ASD"],
                            "ASDl":(glQ["ASD"]),
                            "ASDh":(ghQ["ASD"]),
                            
                            "clip":clip,"cond":typ})
            
       

        if i > 0:
            ax.legend_.remove()
        if i == 4:
            ax.set_xlabel("Bayes factor")

    plt.tight_layout()
    return mediandat

#%% Causal inference

def calc_interventions(data,model,epochs,p,clip,vsamp=1e4):
    from pgmpy.inference import CausalInference
    from pgmpy.factors.discrete import TabularCPD
    model_infer = CausalInference(model)
        
    mvars = list(model.nodes)
    
    nrows = np.sum([len(model.states[node]) for node in mvars])
    ncols = len(mvars)
    
    ratioarray = np.full((nrows, ncols),np.nan)
    
    row_parent_idx = []
    row_labels = []
    row_state_boundaries = []
    c = 0
    
    ratioarray = []
    for n, node in enumerate(mvars):
        start_c = c
        for s,state in enumerate(model.states[node]):
            row_parent_idx.append(n)
            row_labels.append(f"{epochs[n]}={s+1}")
    
            for j, lnode in enumerate(mvars):
                if j <= n:
                    continue
    
                evidence = {node: state}
                ppds1 = model_infer.query([lnode], do=evidence)
                
                # empirical counts
                states, counts = np.unique(
                    data.iloc[:, n], return_counts=True
                )
                #counts = counts / counts.sum()
                counts[s] = 0
                counts = counts / counts.sum()

                counterfactual = TabularCPD(variable=node,
                          variable_card=len(model.states[node]),
                          values=[[counts[k]] for k in range(len(counts))],
                          state_names={node: model.states[node]})
                try:
                    ppds0_samp = model.simulate(n_samples=int(vsamp),virtual_intervention=[counterfactual])
                except:
                    print(counterfactual)
                states0,ppds0 = np.unique(ppds0_samp[lnode],return_counts=True)
                ppds0 = ppds0 / ppds0.sum()
                #ppd_states = list(ppds.state_names.values())[0]
                
                

                ppd1_states = list(ppds1.state_names.values())[0]
                
                ppd0_idx = {state: ppds0[i] for i, 
                            state in enumerate(states0)}
                
                ppds0 = [ppd0_idx[state] if state in states0 else 
                         0 for state in ppd1_states]

                  
                
                #max_idx = np.argmax(all_ppds)
                
                deltaP = ppds1.values - ppds0
          
                #statearray[c, j] = max_idx
                #ratioarray[c, j] = ppds1.values[max_idx]/ppds0[max_idx]#*np.sum(counts[~s]) #all_ppds[max_idx] / counts[max_idx]
                perturbations = {i: delta for i,delta in enumerate(deltaP)}
                ratioarray.append({"p":p, "p_node": node, "c_node":lnode, "state": s+1} | perturbations)
            c += 1
    
        end_c = c
        row_state_boundaries.append((start_c, end_c))
    

    row_parent_idx = np.array(row_parent_idx)
    
    
    return ratioarray,row_parent_idx,ncols,nrows,mvars,row_state_boundaries,epochs,row_labels

def calc_state_changes(data,model,clip,vsamp=1e4,resample="bootstrap",**kwargs):
    import string
    epochs = list(string.ascii_uppercase)
     
    ratioarrays = []
    if "nboots" in kwargs:
        nboots = kwargs["nboots"]
    else:
        nboots = 500
    # jackknife implementation
    if resample == "jackknife":
        for p in range(data.shape[0]):
            data_resampled = data.drop(p)
            ratioarray,row_parent_idx,ncols,nrows,mvars,row_state_boundaries,epochs,row_labels = calc_interventions(data_resampled,model,epochs,p,clip)
            ratioarrays.extend(ratioarray)

    elif resample == "bootstrap":
        for p in range(nboots):
            bootID = np.random.choice(data.shape[0],data.shape[0],replace=True)
            data_resampled = data.iloc[bootID,:]
            model_adj = learn_parameters(data_resampled, model)
            ratioarray,row_parent_idx,ncols,nrows,mvars,row_state_boundaries,epochs,row_labels = calc_interventions(data_resampled,model_adj,epochs,p,clip,vsamp=vsamp)
            ratioarrays.extend(ratioarray)

    return ratioarrays


"""    
def plotStateChanges(statearray,ratioarray,
                     row_parent_idx,ncols,nrows,mvars,
                     row_state_boundaries,
                     epochs,row_labels,mask_array=True,
                     show=True,fig = None,ax = None, **kwargs):
    from matplotlib import patches
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.lines import Line2D

    ###
    
    if show:
        x, y = np.meshgrid(range(ncols), range(nrows))
        x, y = x.flatten(), y.flatten()
        
        state_vals = statearray.flatten()
        ratio_vals = ratioarray.flatten()
        
        mask = (x > row_parent_idx[y]) & ~np.isnan(ratio_vals)
        
        x = x[mask]
        y = y[mask]
        state_vals = state_vals[mask]
        ratio_vals = ratio_vals[mask]
        
        ###
        ratio_vals = np.clip(ratio_vals,0,np.nanpercentile(ratio_vals, 99))
        sf = 30
        sizes = (sf * ratio_vals)**2 
        ##
        state_colors = ["red", "blue", "green", "orange"]
        cmap_states = ListedColormap(state_colors)
        
        n_states = len(state_colors)
        
        norm_states = BoundaryNorm(boundaries=np.arange(-0.5, n_states + 0.5),ncolors=n_states)
        
        if fig is None:
            fig, ax = plt.subplots(figsize=(10, 14))
            ax.set_aspect("equal")
            
        sc = ax.scatter(x, y,
            s=sizes,c=state_vals,cmap=cmap_states,norm=norm_states,edgecolors="none",alpha=0.9,zorder=2)
        
        for i in range(nrows):
            for j in range(row_parent_idx[i] + 1, ncols):
                ax.plot([j - 0.5, j + 0.5], [i - 0.5, i - 0.5], color="lightgray", lw=0.6)
                ax.plot([j - 0.5, j + 0.5], [i + 0.5, i + 0.5], color="lightgray", lw=0.6)
                ax.plot([j - 0.5, j - 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.6)
                ax.plot([j + 0.5, j + 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.6)
        
        v = len(mvars)-1
        for start, end in row_state_boundaries[:-1]:
            yline = end - 0.5
            ax.plot([ncols - v -0.5,len(mvars)],
                [yline, yline],color="black",
                lw=2.5,zorder=3)
            v-=1
        
     
        
        ax.set_xticks(range(ncols))
        ax.set_yticks(range(nrows))
        
        ax.set_xticklabels(epochs[0:len(mvars)])
        ax.set_yticklabels(row_labels)
        
        ax.xaxis.tick_top()
        ax.invert_yaxis()
        
        for spine in ax.spines.values():
            spine.set_visible(False)
        ###################
        color_handles = [
            Line2D(
                [0], [0], marker='o',
                linestyle='None',label=f"State {i+1}",
                markerfacecolor=state_colors[i],markeredgecolor='none',markersize=12)
            for i in range(n_states)
        ]
        
        color_legend = ax.legend(
            handles=color_handles,
            frameon=False,
            bbox_to_anchor=(1.02, 1.0),
            loc="upper left",
            handlelength=0,
            handletextpad=0.5,
            borderpad=0.2
        )
        
        ax.add_artist(color_legend) 
        
        ratio_levels = np.array([0.2,0.5,0.8,])
        
        ratio_levels = np.round(ratio_levels, 2)
        ##
        size_handles = [Line2D([0], [0],
                marker='o',linestyle='None',
                color='gray',alpha=0.8,
                markersize=np.sqrt((sf * r)**2),
                label=f"{r}")
            for r in ratio_levels
        ]
        
        ax.legend(handles=size_handles, title="Ratio",
            frameon=False, bbox_to_anchor=(1.02, 0.55),
            loc="lower left", handletextpad=1.0,
            labelspacing=1.2)
        ######################
        if "overlap" in kwargs:
            rgba = np.zeros((*statearray.shape, 4))
            rgba[kwargs["overlap"] == 1] = [1.0, 0.75, 0.8, 1.0]
    
            ax.imshow(rgba,
                alpha=0.5, zorder=1)
    
    
        plt.tight_layout()
"""


def format_df(df,epochs):
    """
    Parameters
    ----------
    df : pandas.DataFrame
        Wide-format input DataFrame containing:
    epochs : list or array-like
        List of epoch labels or identifiers. 
    Returns
    -------
    df_range : pandas.DataFrame
        Processed DataFrame containing aggregated metrics and plotting metadata
    x_levels : list
        List of unique child node levels used for the x-axis, derived from epochs[1:].
    y_levels : list
        List of unique unique y_label values maintaining original insertion order.
    """
    df_long = df.melt(
        id_vars=["p", "p_node", "c_node", "state"],
        var_name="child_state",
        value_name="delta_p").dropna()

    df_long["child_state"] = df_long["child_state"].astype(int)


    df_range = (df_long
        .groupby(["p_node", "c_node", "state", "child_state"])
        .agg(dmin=("delta_p", lambda x: x.quantile(0.01)),
             dmax=("delta_p", lambda x: x.quantile(0.99)),
             dmedian=("delta_p", lambda x: x.quantile(0.5)))
        .reset_index())
    df_range["plot"] = (0 > df_range["dmin"]) & (0 < df_range["dmax"])
    
    df_range["y_label"] = df_range["p_node"] + "=" + df_range["state"].astype(str)
    y_levels = list(dict.fromkeys(df_range["y_label"]))
    y_map = {k: i for i, k in enumerate(y_levels)}
    df_range["y"] = df_range["y_label"].map(y_map)


    x_levels = epochs[1:]
    x_map = {k: i for i, k in enumerate(x_levels)}
    df_range["x"] = df_range["c_node"].map(x_map)
    return df_range,x_levels,y_levels

import matplotlib.pyplot as plt

def prepare_deltaP_df(epochs,clip):
    
    ratioarrays_TD = load_pickle(f"data\\deltaP_{clip}_TD")
    ratioarrays_ASD = load_pickle(f"data\\deltaP_{clip}_ASD")
    df_TD = pd.DataFrame(ratioarrays_TD)
    df_ASD = pd.DataFrame(ratioarrays_ASD)
    
    df_range_TD,x_levels,y_levels = format_df(df_TD,epochs)
    df_range_ASD,_,_ = format_df(df_ASD,epochs)
    
    
    group_cols = ["p_node", "c_node", "state"]
    summary_TD = (df_range_TD
        .loc[df_range_TD.groupby(group_cols)["dmedian"].idxmax()][group_cols + ["child_state", "plot"]]
        .rename(columns={"child_state": "max_state_TD",
            "plot": "TD_crosses_zero"}))
    
    summary_ASD = (
        df_range_ASD
        .loc[df_range_ASD.groupby(group_cols)["dmedian"].idxmax()]
        [group_cols + ["child_state", "plot"]]
        .rename(columns={"child_state": "max_state_ASD",
                         "plot": "ASD_crosses_zero"}))
    summary = summary_TD.merge(summary_ASD, on=group_cols, how="inner")
    summary["shared"] = (summary["max_state_TD"] == summary["max_state_ASD"])
    
    summary["both"] = (~summary["TD_crosses_zero"] & ~summary["ASD_crosses_zero"])
    
    summary["comp"] = summary["shared"] & summary["both"]

    return df_range_TD,df_range_ASD,summary_TD,summary,x_levels,y_levels


#%%

def plot_deltaPs(df_range,summary,x_levels,y_levels,tile_halfwidth=0.8,scale = 1):
    """
    Plot changes in probability of each state given interventions on prior states, using the object generated from format_df()
    
    Parameters
    ----------
    df_range : pandas.DataFrame
        Processed DataFrame containing aggregated metrics and plotting metadata
    summary : 
        Summary dataframe containing information on parent and child state and whether the 99th % confidence interval crosses 0 
    x_levels : figure data; x positions for drawing rects 
    y_levels : figure data; y positions for drawing rects

    tile_halfwidth float
        determines the width of each tile rect
    scale float
        
    Returns
    -------
    None
    """

    
    fig, ax = plt.subplots(figsize=(10,20))
    
    # color per child state
    colors = {0: "red", 1: "blue", 2: "green", 3: "orange"}
    
    # vertical offsets inside each tile
    offsets = {0: 0.25, 1: 0.08, 2: -0.08, 3: -0.25}
    
    priorset = 0
    for i, row in df_range.iterrows():
        if row["plot"]: continue
        currset = row["p_node"] + row["c_node"] + str(row["state"])
        
        x = row["x"]
        y = row["y"]
    
        if (priorset != currset):
            if summary[(summary["p_node"]==row["p_node"]) & (summary["c_node"]==row["c_node"]) & (summary["state"]==row["state"])]["comp"].item():
                rect = plt.Rectangle((x - 0.45, y - 0.45), 0.9, 0.9, 
                                     facecolor="white", edgecolor="green", linewidth=4)
            else:
                rect = plt.Rectangle((x - 0.45, y - 0.45), 0.9, 0.9, 
                                     facecolor="white", edgecolor="grey", linewidth=4)
            ax.add_patch(rect)
    
        ax.plot([x, x], [y - 0.45, y + 0.45], color="black", linewidth=2.8, alpha=0.6)
        
        ax.plot([x  + (0.5 / scale) * tile_halfwidth, x + (0.5 / scale) * tile_halfwidth], [y - 0.45, y + 0.45],
            color="black", linewidth=2.8, alpha=0.6, linestyle = "dotted")
    
        ax.plot([x  + (-0.5 / scale) * tile_halfwidth, x + (-0.5 / scale) * tile_halfwidth], [y - 0.45, y + 0.45],
            color="black", linewidth=2.8, alpha=0.6, linestyle = "dotted")
        

        xmin = x + (row["dmin"] / scale) * tile_halfwidth
        xmax = x + (row["dmax"] / scale) * tile_halfwidth
        xmed = x + (row["dmedian"] / scale) * tile_halfwidth
        y_offset = y + offsets[row["child_state"]]
    
        # draw horizontal range
        ax.plot([xmin,xmax], [y_offset,y_offset], color=colors[row["child_state"]], linewidth=6)
        ax.scatter(xmed,y_offset, color = colors[row["child_state"]],s=100,zorder=5)
        priorset = currset

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

    



