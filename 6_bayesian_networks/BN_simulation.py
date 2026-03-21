# -*- coding: utf-8 -*-
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
import itertools
from joblib import Parallel, delayed
from pgmpy.models import DiscreteBayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.estimators import ExpertKnowledge
from pgmpy.inference import VariableElimination
from pgmpy.sampling import BayesianModelSampling
from pgmpy.factors.discrete import State

import matplotlib.pyplot as plt
os.chdir("C:\\Users\\cgibeil@emory.edu\\OneDrive - Emory\\Documents\\PhD\\Projects\\P1-CondProbs\\bn_analyses\\6_bayesian_networks\\")
import BN_functions as bn
import logging


class bnSimulation:
    
    def __init__(self, n_vars, spacing, num_parents, 
                 num_states, prior_imbalance = 0, alpha_imbalance = 0, 
                 base_conc = 10, num_samples = 100, seed = None, 
                 score = "AICScore"):
        
        # Variables of interest
        self.num_parents = num_parents
        self.num_states = num_states
        self.prior_imbalance = prior_imbalance
        self.alpha_imbalance = alpha_imbalance
        self.num_samples = num_samples
        
        # Fixed parameters
        self.n_vars = n_vars 
        self.spacing = spacing
        
        # Dirichlet concentration
        self.base_conc =  base_conc
        
        # For reproducibility
        self.seed = seed
        self.score = score

    def create_adjdes(self):
        
        self.adj = np.zeros((self.n_vars, self.n_vars), dtype=int)
        for child in range(self.n_vars):
            for p in range(1, self.num_parents + 1):
                parent = child - p * self.spacing
                if parent >= 0:
                    self.adj[parent, child] = 1
        return self.adj

    def build_alpha_vector(self):
    
        if self.alpha_imbalance == 0:
            weights = np.ones(self.num_states)
        else:
            weights = np.linspace(1,self.alpha_imbalance,self.num_states)
            weights = np.clip(weights, 1e-3, None)
        alpha = self.base_conc * weights
        return alpha


    def sample_cpts_from_adj(self):
        
        rng = np.random.default_rng(self.seed)
        cpts = dict()
        for child in range(self.n_vars):
            parents = np.where(self.adj[:, child] == 1)[0]
            if len(parents) == 0:
                alpha = self.build_alpha_vector()
                probs = rng.dirichlet(np.random.choice(alpha,size=len(alpha),replace=False))
                cpts[str(child)] = {(): probs}
            else:
    
                n_par_conf = self.num_states ** len(parents)
                cpt = {}
    
                for conf_idx in range(n_par_conf):
    
                    conf = []
                    temp = conf_idx
                    for _ in range(len(parents)):
                        conf.append(temp % self.num_states)
                        temp //= self.num_states
                    conf = tuple(conf[::-1])  # reverse to match order
                    alpha = self.build_alpha_vector()
                    probs = rng.dirichlet(np.random.choice(alpha,size=len(alpha),replace=False))
                    cpt[conf] = probs
                cpts[str(child)] = cpt
        self.cpts = cpts
        return self.cpts

    def create_sim_model(self, variables=None):
       
        if variables is None:
            variables = [str(i) for i in range(self.n_vars)]

        edges = [(variables[i], variables[j]) for i in range(self.n_vars) for j in range(self.n_vars) if self.adj[i, j] == 1]

        self.model_ground = DiscreteBayesianNetwork(edges)

        for i, var in enumerate(variables):
            parents = [variables[p] for p in np.where(self.adj[:, i] == 1)[0]]

            # If variable has no parents
            if len(parents) == 0:
                cpt = np.array(list(self.cpts[var].values())).T
                cpd = TabularCPD(
                    variable=var,
                    variable_card=self.num_states,
                    values=cpt
                )
            else:
                parent_cards = [self.num_states for _ in parents]
                cpt = np.array(list(self.cpts[var].values())).T
                cpd = TabularCPD(
                    variable=var,
                    variable_card=self.num_states,
                    values=cpt,
                    evidence=parents,
                    evidence_card=parent_cards
                )

            self.model_ground.add_cpds(cpd)

        self.model_ground.check_model()  
        return self.model_ground

    def simulate_dataset_from_cpts(self):
      
        rng = np.random.default_rng(self.seed)
        data = np.zeros((self.num_samples, self.n_vars), dtype=int)
        
        if self.prior_imbalance == 0:
            weights = np.ones(self.num_states)/self.num_states
        else:
            weights = np.linspace(1,self.prior_imbalance,self.num_states)
            weights /= sum(weights)
            
        for i in range(self.n_vars):
            parents = np.where(self.adj[:, i] == 1)[0]
            if len(parents) == 0:
                #probs = self.cpts[str(i)][()]
                samples = rng.choice(self.num_states, size=self.num_samples, p=weights)
                data[:, i] = samples
            else:
                # compute parent state tuples for all samples and sample accordingly
                parent_indices = parents.tolist()
                # extract parent states for each sample; shape (num_samples, n_parents)
                parent_states = data[:, parent_indices]
                # for each sample, map to tuple and use CPT
                for s in range(self.num_samples):
                    par_tuple = tuple(parent_states[s].tolist())
                    probs = self.cpts[str(i)][par_tuple]
                    data[s, i] = rng.choice(self.num_states, p=probs)
        self.df = pd.DataFrame(data, columns=[str(i) for i in range(self.n_vars)])
        return self.df
    
    def learn_structure_from_data(self, data = [], score="AICScore", show_progress=False, crossval = False):
    
        from pgmpy.estimators import HillClimbSearch, BIC, AIC
        # pgmpy requires string categories or ints; ensure df is suitable
        if len(data) == 0:
            data = self.df.copy().astype(int)
        
        #plt.imshow(data)
        #plt.show()

        state_names = {data.columns[i]: np.unique(data.iloc[:, i]) for i in range(self.n_vars)}
    
        black_list= bn.generate_forbidden_edges(data.columns)
        
        forbidden = ExpertKnowledge(forbidden_edges=black_list)
        
        hc = HillClimbSearch(data, state_names=state_names)
        if score == "BICScore":
            score_obj = BIC(data)
        elif score == "AICScore":
            score_obj = AIC(data)

        dag = hc.estimate(scoring_method=score_obj, expert_knowledge= forbidden, show_progress=False)
        # build adjacency
        
        if not crossval: 
            self.adj_rec = bn.create_adj(dag,data)
            self.model_rec = bn.learn_parameters(data,dag)
            
            return self.adj_rec,self.model_rec
        elif crossval:
            adj_rec = bn.create_adj(dag,data)
            model_rec = bn.learn_parameters(data,dag)
            
            return adj_rec,model_rec    
        
    def compare_recovered_ppds(self,fun=np.mean):
        variables = self.model_ground.nodes()
        model_groundInf = VariableElimination(self.model_ground)
        model_fitInf = VariableElimination(self.model_rec)
        jsds = []
        for i, node in enumerate(variables):
            parents = self.model_rec.get_parents(node)
            par_states = [self.model_rec.states[str(k)] for k in list(variables) if str(k) in parents]
            combs = itertools.product(*par_states)
            
            for c, comb in enumerate(combs):

                evidence = {parents[p]: state_val for p,state_val in enumerate(comb)}
                
                query_fit = model_fitInf.query(variables=[node],evidence=evidence)
                query_ground = model_groundInf.query(variables=[node],evidence=evidence)
                diff = len(query_fit.values) - len(query_ground.values)
                if diff == 0:   
                    jsds.append(bn.js_distance(query_fit.values, query_ground.values))
                elif diff < 0:
                    
                    g_states = list(query_ground.state_names.values())[0]
                    f_states = list(query_fit.state_names.values())[0]
                    
                    f_idx = {state: i for i, state in enumerate(f_states)}
                    
                    query_fit.values = [
                        query_fit.values[f_idx[state]] if state in f_idx else 0
                        for state in g_states
                    ]

                    jsds.append(bn.js_distance(query_fit.values, query_ground.values))
                elif diff > 0:
                    g_states = list(query_ground.state_names.values())[0]
                    f_states = list(query_fit.state_names.values())[0]
                    
                    g_idx = {state: i for i, state in enumerate(g_states)}
                    
                    query_ground.values = [
                        query_ground.values[g_idx[state]] if state in g_idx else 0
                        for state in f_states
                    ]
                    jsds.append(bn.js_distance(query_fit.values, query_ground.values))                
        jsd_fn = fun(jsds)
        return jsd_fn
    
    def hamming_distance(self,directed=True):

        if not directed:
            T = np.triu(self.adj + self.adj.T, 1)
            P = np.triu(self.adj_rec + self.adj_rec.T, 1)
            diff = (T != P).astype(int)
            possible = T.size // 2
            raw = diff.sum()
        else:
            diff = (self.adj != self.adj_rec).astype(int)
            raw = diff.sum()
            possible = self.adj.size
        return raw / possible
    
    def precision_recall(self):
        #plt.imshow(self.adj_rec)
        #plt.show()
        #plt.imshow(self.adj)
        #plt.show()
        
        triu_mask = np.triu(np.ones_like(self.adj, dtype=bool), k=1)
        
        ground_upper = (self.adj == 1) & triu_mask
        rec_upper = (self.adj_rec == 1) & triu_mask
        
        tp = np.sum(ground_upper & rec_upper)  # edges present in both
        fp = np.sum(~ground_upper & rec_upper)  # recovered but not ground
        fn = np.sum(ground_upper & ~rec_upper)  # ground but not recovered
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * ((precision*recall) / (precision + recall + 1e-8) )
        
        return precision, recall, f1
        
        
    def predict(self, size=50, n_jobs=-1, posterior_sampling=False):

        pred = np.zeros([self.num_samples, self.n_vars])
        if posterior_sampling:
            def predict_row(i, samp, model_rec, adj_rec, num_states, n_vars):
                
                
                from pgmpy.sampling import BayesianModelSampling  # re-import inside subprocess
                data_train = self.df.drop(i)
                _,model_rec = self.learn_structure_from_data(data = data_train,crossval=True)
                inference = BayesianModelSampling(model_rec)
                row_pred = np.zeros(n_vars)
        
                for j in range(n_vars):
                    parents = np.where(adj_rec[:, j] == 1)[0]
                    if len(parents) == 0:
                        continue
        
                    parent_indices = parents.tolist()
                    parent_states = samp.iloc[parent_indices]
        
                    evidence = [
                        State(var=str(pid), state=int(parent_states.iloc[n]))
                        for n, pid in enumerate(parent_indices)
                    ]
        
                    try:
                        df_sim = inference.likelihood_weighted_sample(
                            evidence=evidence, size=size, show_progress=False
                        )
                        pred_state = df_sim.iloc[:, j].mode().iloc[0]
                        row_pred[j] = float(pred_state == samp.iloc[j])
                    except Exception:
                        row_pred[j] = np.nan
        
                return i, row_pred
        
            results = Parallel(n_jobs=n_jobs, verbose=False)(
                delayed(predict_row)(
                    i, samp, self.model_rec, self.adj_rec, self.num_states, self.n_vars
                )
                for i, samp in self.df.iterrows()
            )
        
            for i, row_pred in results:
                pred[i, :] = row_pred
        
            return np.nansum(pred, axis=0)/self.num_samples
        else:
            inference = VariableElimination(self.model_rec)
            nodes = self.model_rec.nodes()
                        
            def predict_row(i, samp, model_rec, inference, nodes, n_vars):
                row_pred = np.zeros(n_vars)
                for j, node in enumerate(nodes):
                    parents = list(model_rec.get_parents(node))
                    if not parents:
                        cpd = model_rec.get_cpds(node)
                        most_prob_state = np.argmax(cpd.values)
                        row_pred[j] = int(samp[node] == most_prob_state)
                    else:
                        evidence = {p: samp[p] for p in parents}
                        q = inference.query([node], evidence=evidence, show_progress=False)
                      
                        pred_state = q.values.argmax()
                        row_pred[j] = int(samp[node] == pred_state)
                return i, row_pred
            
            results = Parallel(n_jobs=n_jobs, verbose=False)(
                delayed(predict_row)(i, samp, self.model_rec, inference, nodes, self.n_vars)
                for i, samp in self.df.iterrows()
            )
        
            for i, row_pred in results:
                pred[i, :] = row_pred
    
            return np.nansum(pred, axis=0)/self.num_samples
            
   
    def run(self):
        self.create_adjdes()
        self.build_alpha_vector()
        self.sample_cpts_from_adj()
        self.simulate_dataset_from_cpts()
        self.create_sim_model()
        
        self.learn_structure_from_data(score=self.score)

        precision,recall,f1 = self.precision_recall()

        
        hd = self.hamming_distance()
            
        ppd = self.compare_recovered_ppds()
        
        pred = self.predict()
        return self.df, self.adj, self.adj_rec, self.model_ground, self.model_rec, precision, recall, f1, hd, ppd, pred 
        
            

def run_single_experiment(n_vars = 15, spacing = 1, num_parents = 2, num_states = 3, 
                          prior_imbalance = 0, alpha_imbalance = 0, 
                          base_conc = 10, num_samples = 100, seed = None,return_data=False,score="AICScore"):

    sim = bnSimulation(n_vars=n_vars,
                   spacing=spacing,
                   num_parents=num_parents,
                   num_states=num_states,
                   prior_imbalance=prior_imbalance,
                   alpha_imbalance=alpha_imbalance,
                   base_conc=base_conc,
                   num_samples=num_samples,
                   seed=seed,score=score)
    if return_data:
        data, adj, adj_rec, model_ground, model_rec, precision, recall, f1, hd, ppd, pred = sim.run()
        return {
            "adj_rec": adj_rec,
            "model_ground": model_ground,
            "model_rec": model_rec,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "hd": (hd),
            "ppd": ppd,
            "pred": pred,
            "data":data
        }

    else:
        _,adj, adj_rec, model_ground, model_rec, precision, recall, f1, hd, ppd, pred = sim.run()
        return {
            "adj_rec": adj_rec,
            "model_ground": model_ground,
            "model_rec": model_rec,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "hd": (hd),
            "ppd": ppd,
            "pred": pred
        }
    
def bootstrap_parts(data,nboot,thresh):
    nparts = data.shape[0]
    #partstotal = np.arange(1,nparts,1)
    #nretain = round(nparts*pretain)
    adjs = {}

    for i in range(nboot):
        bootID = np.random.choice(nparts,nparts,replace=True)
        sampleboot = data.iloc[bootID,:]
        
        #sample = np.random.choice(partstotal,nretain,replace=False)
        #dataBoot = data.iloc[sample]
        _,adjs[i],_,_ = bn.learn_structure(sampleboot,"AICScore")
    
    adjAll = bn.threshold_bootstrapAdj(adjs, thresh)
    return adjs,adjAll


def plot_heatmaps(param_grid, results,xvar,yvar,zvar, xlabel, ylabel, zlabel, vmin,vmax, fixed=None, ngrid=80,fig=None,ax=None,interpolate=False):
    
    from scipy.interpolate import griddata
    import seaborn as sns

    plt.rcParams['figure.figsize'] = (8, 6) 
    plt.rcParams['lines.linewidth'] = 3    
    plt.rcParams['font.size'] = 40
    param_grid['num_participants'] = param_grid['num_samples']
    

    df = pd.pandas.DataFrame(results)
    
    mask = pd.Series(True, index=df.index)
    if fixed:
        for col, value in fixed.items():
            mask = mask & (df[col] == value)
        
        df = df[mask]
    agg = df.groupby([xvar, yvar], as_index=False)[zvar].mean()
    agg[zvar] = agg[zvar].apply(np.mean)


    if fig == None and ax == None:
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111)

    
    x = agg[xvar].values.astype(float)
    y = agg[yvar].values.astype(float)
    z = agg[zvar].values.astype(float)
    
    xshape = len(param_grid[xvar])
    yshape = len(param_grid[yvar])
    
    z_grid = np.reshape(z, [xshape, yshape]).T

    z = z_grid.T.ravel()
    
    # Create fine grid in XY-plane
    xi = np.linspace(x.min(), x.max(), ngrid)
    yi = np.linspace(y.min(), y.max(), ngrid)
    Xi, Yi = np.meshgrid(xi, yi)

    Zi = None

    Zi = griddata((x, y), z, (Xi, Yi), method="cubic")

    # mask NaNs for plotting (so no spurious artifacts)
    mask = np.isnan(Zi)
    Zi_masked = np.ma.array(Zi, mask=mask)

    cbar_kws = {"label": zlabel}
    
    if not interpolate:
        xshape = len(param_grid[xvar])
        yshape = len(param_grid[yvar])
        z = np.reshape(z,[xshape,yshape]).T
        sns.heatmap(z,cmap="inferno",cbar_kws=cbar_kws,ax=ax,vmin=vmin,vmax=vmax)
        ax.set_xticklabels(param_grid[xvar])
        ax.set_yticklabels(param_grid[yvar])
    else:
        sns.heatmap(Zi_masked,cmap="inferno",cbar_kws=cbar_kws,ax=ax,vmin=vmin,vmax=vmax)
        posx = np.round((np.array(param_grid[xvar])-min(param_grid[xvar]))/(max(param_grid[xvar])-min(param_grid[xvar])))
        ax.set_xticks(posx)
        ax.set_xticklabels(param_grid[xvar])
    
        posy = np.round((np.array(param_grid[yvar])-min(param_grid[yvar]))/(max(param_grid[yvar])-min(param_grid[yvar])))
        ax.set_yticks(posy)
        ax.set_yticklabels(param_grid[yvar])

    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def plot_mi_adjacency_sim(it,results_single_pd,scores):
    
    fig, axs = plt.subplots(2,4, figsize=(45, 15)) 
    for i, (npar, nstat, nepoch,nsamp) in enumerate(it):
        results_loop = results_single_pd[(results_single_pd["num_parents"]==npar) &
                                         (results_single_pd["num_states"]==nstat)]
        for j,sc in enumerate(scores):

            mi_dict = bn.calc_mi_model(results_loop["data"].iloc[j],results_loop["model_rec"].iloc[j])
            
            adjBoot = sum(results_loop["adjBoot"].iloc[j].values())/len(results_loop["adjBoot"].iloc[j]) > 0.54
            modelBoot = bn.learn_parameters(results_loop["data"].iloc[j], bn.adj2DAG(adjBoot))
            
            ci_array = bn.find_ci_edges(modelBoot, results_loop["data"].iloc[j])
            
            bn.plot_mi_adjacency(modelBoot,
                                 adjBoot,ci_array,mi_dict,results_loop["adjBoot"].iloc[j],
                                 fig=fig,ax=axs[j,i])
            axs[j,i].text(s="Score: " + sc[0:3], x = 0, y = len(results_loop["adj_rec"].iloc[j])*0.5)
            axs[j,i].text(s="N parents: " + str(npar), x = 0, 
                          y = len(results_loop["adj_rec"].iloc[j])*0.6)
            axs[j,i].text(s="N states: " + str(nstat), x = 0, 
                          y = len(results_loop["adj_rec"].iloc[j])*0.7)










