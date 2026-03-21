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
def get_bootstrapped_samples(data,nboot=1000,thresh=0.54):
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
    mi_dict : TYPE
        DESCRIPTION.
    adjBoot : TYPE
        DESCRIPTION.

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
    mask = np.triu(np.ones((n, n)), k=0).flatten().astype(bool)
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
             
             plot_mi_adjacency(dat_loop["model"].iloc[0], 
                               dat_loop["adj"].iloc[0],
                               ci_array,mi_dict, 
                               adjBoots[dx][clip],
                               cmap = col[j],
                               fig=fig,ax=axs[j,i])
             
             axs[j,i].text(s="Movie " + str(i+1), x = 0, y = len(dat_loop["adj"].iloc[0])*0.6)
             
             mi, dists =  calc_mi_dist(mi_dict)
             mi_dat.append({"clip":clip,"dx":dx,"mis":mi, 
                            "dists":dists
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
    
    from pgmpy.estimators import MaximumLikelihoodEstimator, ExpectationMaximization
    
    model = DiscreteBayesianNetwork(dag)
    model.fit(data, estimator=MaximumLikelihoodEstimator)

    return model



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
    from scipy.spatial import distance
    if np.isnan(p).any() or np.isnan(q).any():
        return np.nan
    else: 
        return distance.jensenshannon(p, q)

def compare_recovered_ppds(model1,model2, fun=np.mean):
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

#%% Model prediction

def predictBayes(data, adj, model, size=10000, n_jobs=-1,
            cond = "evidence",alpha=1,seed = None):
    """
    Parameters
    ----------
    data : state-matrix dataframe object 
    adj : adjacency matrix 
    model : the BN model for which prediction is performed
    size : size of likelihood-weighted sample to pull for each prediction
    n_jobs : optional, specifies the number of computing cores to use
    cond : optional, which model setup to use, M0, M1 or M2
    alpha : the prior set on each state. By default = 1
    seed : optional, seed is set `get_all_predictions` for reproducibility

    Returns
    -------
    pred : set of all predictions for the query movie clip and group
    """
    from pgmpy.factors.discrete import State
    from joblib import Parallel, delayed
    from pgmpy.sampling import BayesianModelSampling
    

    pred = np.zeros([data.shape[0], data.shape[1]])
    model_preds= []
    if cond == "forder": # M1 condition with only first-order edges
        # Convert model adjacency matrix into a numpy array
        adj = DAG2adj(model)
        # Extract only the first-order edges from the model
        adj[~np.eye(adj.shape[0], dtype=bool,k=1)] = 0
        # Convert adjacency matrix back into DAG
        dag = adj2DAG(adj)

        # LOO parameter learning. Using the same bootstrapped structure, 
        # refit the model parameters leaving each participant out a time 
        for i in range(data.shape[0]): # for each participant
            
            data_train = data.drop(i) # drop that participant from the training 
            model_preds.append(learn_parameters(data_train, dag))

        
        
    elif cond == "evidence" or cond == "rand": # M2 or M0 conditions

        for i in range(data.shape[0]):
            data_train = data.drop(i)
            model_preds.append(learn_parameters(data_train, adj2DAG(adj)))

    def predict_row(i, samp, model_rec, adj_rec,data_train): # predict each 
        # participant's state sequence
        

        rng = np.random.default_rng(seed)
        
        inference = BayesianModelSampling(model_rec)
        
        row_pred = np.zeros(data.shape[1])
        
        epochs = data_train.columns

        for j,node in  enumerate(model_rec.nodes): # For each epoch, 

            if cond == "rand": # For M0, find the frequencies of all the states in the epoch
                states,prior_dist = np.unique(data_train.iloc[:,j],
                                          return_counts= True)
                #Specify alpha as the prior set on the null model for each state 
                #Dirichlet prior is the conjugate prior of the categorical distribution. The 
                #Updated probabilities are given for each state i as 
                #n_{state,i} + alpha / sum(n_{states}) + (alpha)
                
                null_posterior = (prior_dist+alpha)/(np.sum(prior_dist) + states.size * alpha)
                
                # which state did the participant enter
                
                # log-likelihood
                samp_state = str(samp.iloc[j])
                mask = states == samp_state
                idx = np.where(mask)[0]
                
                if len(idx) == 1:
                    row_pred[j] = np.log10(null_posterior[idx[0]]+1e-16)
                else:
                    row_pred[j] = np.log10(1e-16) 
                    
            else: # if not M0 (M1/M2)
                # for the first epoch with no history
                if j == 0: #len(parents) == 0:
                    
                    # generate a simulated dataframe given the model parameters.
                    df_sim = inference.likelihood_weighted_sample(
                         size=size, show_progress=False
                    )
                
                else:
                    # for subsequent epochs
                    history = epochs[0:j]
                    parent_states = samp.iloc[:]
                    evidence = [
                        State(var=str(pid), state=parent_states.iloc[n])
                        for n, pid in enumerate(history)
                    ]
                     
                    # generate likelihood-weighted samples after setting the participant's 
                    # state history 
    
                    df_sim = inference.likelihood_weighted_sample(
                         size=size, show_progress=False, evidence=evidence
                    )
                    
               
                samples = df_sim.iloc[:, j].astype(str)
                weights = df_sim["_weight"].values
                
                post_dist = {}
                for s, w in zip(samples, weights):
                    post_dist[s] = post_dist.get(s, 0.0) + w

                
                Z = sum(post_dist.values())
                if Z <= 0 or not np.isfinite(Z):
                    # Default to null if likelihoods low
                    states,prior_dist = np.unique(data_train.iloc[:,j],
                                              return_counts= True)
                    # as before for M0
                    null_posterior = (prior_dist+alpha)/(np.sum(prior_dist) + 
                                                         states.size * alpha)
                    
                    mask = states == samp_state
                    idx = np.where(mask)[0]
                    
                    if len(idx) == 1:
                        row_pred[j] = np.log10(null_posterior[idx[0]]+1e-16)
                    else:
                        row_pred[j] = np.log10(1e-16) 
                    
                    continue                
                
                else:
                    for s in post_dist:
                        post_dist[s] /= Z
                
                samp_state = str(samp.iloc[j])
                prob = post_dist.get(samp_state, 0.0)
  
                
                row_pred[j] = np.log10(prob+1e-16)

        return i, row_pred

    # Parallel processing 
    results = Parallel(n_jobs=n_jobs, verbose=False)(
        delayed(predict_row)(
            i, samp, model_preds[i], adj, data_train)
        
    for i, samp in data.iterrows()
    )

    for i, row_pred in results:

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
            seed = np.random.randint(0,10000)
            
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
          
            
            # AD ############################################################
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

def calcStateChanges(model,data,mask_array=True,
                     show=True,vsamp=1e4,thresh=0.2,**kwargs):
    

    from pgmpy.inference import CausalInference
    from pgmpy.factors.discrete import TabularCPD
    import string
    model_infer = CausalInference(model)
    
    epochs = list(string.ascii_uppercase)
    
    mvars = list(model.nodes)
    
    nrows = np.sum([len(model.states[node]) for node in mvars])
    ncols = len(mvars)
    
    statearray = np.full((nrows, ncols),np.nan)
    ratioarray = np.full((nrows, ncols),np.nan)
    
    row_parent_idx = []
    row_labels = []
    row_state_boundaries = []
    
    c = 0
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
                
                ppds0_samp = model.simulate(n_samples=int(vsamp),virtual_intervention=[counterfactual])
                states0,ppds0 = np.unique(ppds0_samp[lnode],return_counts=True)
                ppds0 = ppds0 / ppds0.sum()
                #ppd_states = list(ppds.state_names.values())[0]
    
                

                ppd1_states = list(ppds1.state_names.values())[0]
                
                ppd0_idx = {state: ppds0[i] for i, state in enumerate(states0)}

                ppds0 = [
                    ppd0_idx[state] if state in states0 else 0
                    for state in ppd1_states
                ]
                
                  
                
                #max_idx = np.argmax(all_ppds)
                
                deltaP = ppds1.values - ppds0
                #max_idx = np.argmax(ppds1.values)
                max_idx = np.argmax(deltaP)
                
                statearray[c, j] = max_idx
                #ratioarray[c, j] = ppds1.values[max_idx]/ppds0[max_idx]#*np.sum(counts[~s]) #all_ppds[max_idx] / counts[max_idx]
                ratioarray[c, j] = deltaP[max_idx]
            c += 1
    
        end_c = c
        row_state_boundaries.append((start_c, end_c))
    
    if mask_array:
        neutral = ratioarray < thresh #(ratioarray < 1.5) #& (ratioarray < 1.1)
        statearray[neutral] = np.nan
        ratioarray[neutral] = np.nan
    
    row_parent_idx = np.array(row_parent_idx)
    return statearray,ratioarray,row_parent_idx,ncols,nrows,mvars,row_state_boundaries,epochs,row_labels
    
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
        ratio_vals = np.clip(
            ratio_vals,
            0,
            np.nanpercentile(ratio_vals, 99)
        )
        sf = 30
        sizes = (sf * ratio_vals)**2 
        ##
        state_colors = ["red", "blue", "green", "orange"]
        cmap_states = ListedColormap(state_colors)
        
        n_states = len(state_colors)
        
        norm_states = BoundaryNorm(
            boundaries=np.arange(-0.5, n_states + 0.5),
            ncolors=n_states
        )
        
        if fig is None:
            fig, ax = plt.subplots(figsize=(10, 14))
            ax.set_aspect("equal")
            
        sc = ax.scatter(
            x,
            y,
            s=sizes,
            c=state_vals,
            cmap=cmap_states,
            norm=norm_states,
            edgecolors="none",
            alpha=0.9,
            zorder=2
        )
        
        for i in range(nrows):
            for j in range(row_parent_idx[i] + 1, ncols):
                ax.plot([j - 0.5, j + 0.5], [i - 0.5, i - 0.5], color="lightgray", lw=0.6)
                ax.plot([j - 0.5, j + 0.5], [i + 0.5, i + 0.5], color="lightgray", lw=0.6)
                ax.plot([j - 0.5, j - 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.6)
                ax.plot([j + 0.5, j + 0.5], [i - 0.5, i + 0.5], color="lightgray", lw=0.6)
        
        v = len(mvars)-1
        for start, end in row_state_boundaries[:-1]:
            yline = end - 0.5
            ax.plot(
                [ncols - v -0.5,len(mvars)],
                [yline, yline],
                color="black",
                lw=2.5,
                zorder=3
            )
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
                [0], [0],
                marker='o',
                linestyle='None',
                label=f"State {i+1}",
                markerfacecolor=state_colors[i],
                markeredgecolor='none',
                markersize=12
            )
            for i in range(n_states)
        ]
        """
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
        """
        ratio_levels = np.array([
            0.2,
            0.5,
            0.8,
        ])
        
        ratio_levels = np.round(ratio_levels, 2)
        ##
        size_handles = [
            Line2D(
                [0], [0],
                marker='o',
                linestyle='None',
                color='gray',
                alpha=0.8,
                markersize=np.sqrt(
                    (sf * r)**2 
                ),
                label=f"{r}"
            )
            for r in ratio_levels
        ]
        
        ax.legend(
            handles=size_handles,
            title="Ratio",
            frameon=False,
            bbox_to_anchor=(1.02, 0.55),
            loc="lower left",
            handletextpad=1.0,
            labelspacing=1.2
        )
        ######################
        if "overlap" in kwargs:
            rgba = np.zeros((*statearray.shape, 4))
            rgba[kwargs["overlap"] == 1] = [1.0, 0.75, 0.8, 1.0]
    
            ax.imshow(
                rgba,
                alpha=0.5, zorder=1
                )
    
    
        plt.tight_layout()


def quantify_cons_overlap(causal_infTD,causal_infASD,overlap):
    n_epochs = causal_infTD[3]
    overlap_col = np.zeros((n_epochs,n_epochs))
    td_col =np.zeros((n_epochs,n_epochs))
    asd_col = np.zeros((n_epochs,n_epochs))
    for j in range(n_epochs):
        overlap_col[j,:] += np.nansum(overlap[causal_infTD[2]==j,:],0)
        td_col[j,:] += np.nansum(causal_infTD[1][causal_infTD[2]==j,:]>0,0)
        asd_col[j,:] += np.nansum(causal_infASD[1][causal_infASD[2]==j,:]>0,0)
        
    td_col -= overlap_col
    asd_col -= overlap_col
    
    overlap_dist = np.zeros(n_epochs-1)
    td_dist = np.zeros(n_epochs-1)
    asd_dist = np.zeros(n_epochs-1)
    
    for j in range(n_epochs-1):
        overlap_dist[j] = np.sum(np.linalg.diagonal(overlap_col, offset=j+1))
        td_dist[j] = np.sum(np.linalg.diagonal(td_col, offset=j+1))
        asd_dist[j] = np.sum(np.linalg.diagonal(asd_col, offset=j+1))
    
    return n_epochs,overlap_dist,td_dist,asd_dist

def plotStateChanges_perclip(datadir,clips,venn=False):
    dist_deltas = []
    

    for c, clip in enumerate(clips):
        fig, axs = plt.subplots(1,2, figsize=(30,50)) 
        causal_infTD = load_pickle(f"{datadir}causal_infTD_{str(clip)}")
        causal_infASD = load_pickle(f"{datadir}causal_infASD_{str(clip)}")
        
        overlap = (causal_infTD[0]==causal_infASD[0]) & ~np.isnan(causal_infTD[0]) 
        
        plotStateChanges(*causal_infTD,overlap=overlap,
                            fig=fig,ax = axs[0])
        plotStateChanges(*causal_infTD,overlap=overlap,
                            fig=fig,ax = axs[1])
        if venn:
            from matplotlib_venn import venn2, venn2_circles
            n_overlap = np.sum(overlap)
            n_TD = np.sum(causal_infTD[1]>0)
            n_ASD = np.sum(causal_infASD[1]>0)
            subsets = (n_TD, n_ASD, n_overlap) 
            plt.figure(figsize=(15,15))
            venn2(subsets=subsets, set_labels=["TD","ASD"], set_colors=("blue", "green"))
            venn2_circles(subsets=subsets, linestyle="dashed", linewidth=5) 
            plt.show()
        n_epochs,overlap_dist,td_dist,asd_dist = quantify_cons_overlap(causal_infTD,
                                                              causal_infASD,overlap)
        dist_deltas.append({
            "clip":clip,"dist":np.arange(0,n_epochs-1,1), "overlap_dist":overlap_dist, "td_dist":td_dist,
            "asd_dist":asd_dist})
        plt.show()
    return dist_deltas
            
