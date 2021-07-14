import pandas as pd
import cobra
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import linprog
import scipy
import mip
from copy import deepcopy

def create_stoichiometry_matrix(model):
    """ This creates a stoichiometry matrix"""
    metabolites = model.metabolites 
    reactions = model.reactions 
    S = np.zeros((len(metabolites), len(reactions)))
    
    met_id = dict()
    rec_id = dict()
    for i,reaction in enumerate(model.reactions):
        rec_id[reaction.id] = i
        for metabolite, stoich in reaction.metabolites.items():
            met_id[metabolite.id] = int(metabolites.index(metabolite))
            S[metabolites.index(metabolite), i] = stoich
    return S, met_id, rec_id 


class Model():
    def __init__(self, model, biomass_function):
        """ This is a new class of metabolic model, capable of flux balance analysis
        Attributes:
        models (list): CobraPy models of single organisms which will be used in construction
        biomass_reactions (list): List of strings containing the ids for the growth reactions
        """
        self.biomass_function = biomass_function
        self.model = model
        self.id = model.id
        # Compute stoichimetry_matrix
        S, met_id, rec_id = create_stoichiometry_matrix(model)
        self.num_reactions = S.shape[1]
        self.num_metabolites = S.shape[0]
        self.stoichiometry_matrix = scipy.sparse.csr_matrix(S)
        self.met_id = met_id
        self.rec_id = rec_id 
        # Set objective
        idx = self.rec_id[biomass_function]
        c = np.zeros(self.num_reactions)
        c[idx] = 1
        self.objective_c = c
        # Set bounds
        self._reset_bounds()
    
    @property
    def reactions(self):
        return self.model.reactions
    @property
    def exchanges(self):
        return self.model.exchanges
    @property
    def metabolites(self):
        return self.model.metabolites
    @property
    def medium(self):
        return self.model.medium

    def set_medium(self, medium):
        ex_ids = [ex.id for ex in self.exchanges]
        new_med = {}
        for key,val in medium.items():
            if key in ex_ids:
                new_med[key] = val
        self.model.medium = new_med
        self._reset_bounds()
        
    def optimize(self, disp=False):
        sol = linprog(-self.objective_c, A_eq=self.stoichiometry_matrix, b_eq=np.zeros(self.num_metabolites), bounds=self.bounds, method="highs", options={"disp":disp})
        sol["fun"] = -sol["fun"] # As we have to minimize
        return sol 
    
    def slim_optimize(self, disp=False):
        sol = self.optimize(disp=disp)
        return sol["fun"]

    def summary(self):
        sol = self.optimize()
        flux = sol["x"]
        ex_ids = [ex.id for ex in self.exchanges]
        fluxes = []
        for ex in ex_ids:
            idx = self.rec_id[ex]
            fluxes.append(flux[idx])
        summary_df = pd.DataFrame({"Exchange reaction": ex_ids, "Flux": fluxes})
        summary_df.sort_values(["Flux"], inplace=True)
        return summary_df

    def _reset_bounds(self):
        self.bounds = []
        for rec in self.model.reactions:
            self.bounds.append((rec.lower_bound, rec.upper_bound))

    def __add__(self, model2):
        """ Adding another model creates a community model """
        return CommunityModel([self,model2], [1.,1.])

class CommunityModel(Model):
    def __init__(self, models, weight, shared_exchanges=None):
        self.models = models
        self.id = "|".join([model.id for model in models])
        self.weight = weight
        if shared_exchanges == None:
            self.shared_exchanges = []
            for model in models:
                for rec in model.exchanges:
                    if rec.id not in self.shared_exchanges:
                        self.shared_exchanges.append(rec.id)
        else:
            self.shared_exchanges = shared_exchanges
        # Create community stoichimetry matrix with shuttel reactions!
        self._shifts = [len(self.shared_exchanges)]
        for i,model in enumerate(models):
            self._shifts.append(self._shifts[i] + model.num_reactions)
        S_EX = -np.eye(self._shifts[0])
        matrices = [S_EX] + [m.stoichiometry_matrix.todense() for m in models]
        S = scipy.linalg.block_diag(*matrices)
        self.num_reactions = S.shape[1]
        self.num_metabolites = S.shape[0]
        for i, id in enumerate(self.shared_exchanges):
            for j,model in enumerate(models):
                if id in model.rec_id:
                    S[i,self._shifts[j] + model.rec_id[id]] = 1
        self.stoichiometry_matrix = scipy.sparse.csr_matrix(S)
        # Cretae objective:
        self._weighted_objective(weight)
        # Create bounds
        self._reset_bounds()
    
    @property
    def reactions(self):
        return [model.reactions for model in self.models]
    @property
    def exchanges(self):
        return self.shared_exchanges
    @property
    def metabolites(self):
        return [model.metabolites for model in self.models]
    @property
    def medium(self):
        medium = {}
        for model in self.models:
            for key,val in model.medium:
                medium[key] = val
        return medium

    def _weighted_objective(self, weight):
        self.weight = weight
        self.objective_c = np.zeros(self._shifts[0])
        for i,model in enumerate(self.models):
            self.objective_c = np.append(self.objective_c, weight[i]*model.objective_c)
    def _reset_bounds(self):
        self.bounds = []
        for id in self.shared_exchanges:
            min_lower_bound = 0
            for model in self.models:
                if id in model.rec_id:
                    rec = model.reactions.get_by_id(id)
                    if rec.lower_bound < min_lower_bound:
                        min_lower_bound = rec.lower_bound 
            self.bounds.append((min_lower_bound, 1000))
        for model in self.models:
            self.bounds += model.bounds

    def set_medium(self, medium):
        for model in self.models:
            model.set_medium(medium)
        self._reset_bounds()

    def compute_alpha_weightings(self, alpha, maxbiomass=0.1):
        assert alpha <= 1 and alpha >= 0
        assert len(self.models) == 2
        alphas = np.array([alpha,1-alpha])
        # Alpha objective...
        c = np.zeros(self._shifts[0])
        for i, model in enumerate(self.models):
            c = np.append(c, model.objective_c)
        # Biomoss constraints
        c_mask = (c > 0)
        A_ub = np.zeros((2,len(c)))
        A_ub[0,c_mask] = np.array([1,0])
        A_ub[1,c_mask] = np.array([0,1])
        b_ub = alphas*maxbiomass
        sol = linprog(-c, A_eq=self.stoichiometry_matrix, b_eq=np.zeros(self.num_metabolites), A_ub = A_ub, b_ub=b_ub, bounds=self.bounds, method="highs")
        fluxes = sol["x"]
        growths = fluxes[c > 0]
        summary = self.summary(sol)
        return growths, fluxes,summary

    def summary(self, sol=None):
        if sol == None:
            sol = self.optimize()
        flux = sol["x"]
        ex_ids = self.shared_exchanges
        ex_flux = flux[:len(ex_ids)]
        df_ex = pd.DataFrame({"Exchanges": ex_ids, "Flux": ex_flux})
        df_ex.sort_values(["Flux"], inplace=True)
        for i,model in enumerate(self.models):
            shuttel_ids = self.shared_exchanges
            id = str(model.id) + " Shuttel Flux"
            df_ex[id] = 0.
            for sh in shuttel_ids:
                idx = model.rec_id[sh]
                df_ex[id][ex_ids.index(sh)] = flux[self._shifts[i] +idx]
        
        return df_ex
    
    def set_weights(self, weight):
        self._weighted_objective(weight)
    
    def get_model_growths(self):
        mask = self.objective_c != 0
        sol = self.optimize()
        flux = sol["x"]
        return flux[mask]
    

from mip import xsum, maximize, BINARY
import mip

class MIP_community_model():
    def __init__(self, model1, model2, weights = [1.,1.]):
        # Necessary information of model1
        self.S1 = model1.stoichiometry_matrix.todense()
        self.S1_dict = model1.rec_id 
        self.bounds1 = model1.bounds 
        self.obj1 = np.where(model1.objective_c > 0)[0][0]

        self.S2 = model2.stoichiometry_matrix.todense()
        self.S2_dict = model2.rec_id 
        self.bounds2 = model2.bounds
        self.obj2 = np.where(model2.objective_c > 0)[0][0]

        self.medium = dict(list(model1.medium.items()) + list(model2.medium.items()))

        self.comm_model = mip.Model()
        self.build_mip_model()
        self.comm_model.write("model.lp")

        self.weights = weights

    def reset_model(self):
        model = mip.Model()
        model.read("model.lp")
        for key in self.medium:
            x = model.var_by_name(key)
            self.shuttel_reactions[key] = x
        x1 = []
        for i in range(len(self.bounds1)):
            x = model.var_by_name("x1" + str(i))
            x1 += [x]
        x2 = []
        for i in range(len(self.bounds2)):
            x = model.var_by_name("x2" + str(i))
            x2 += [x]
        self.x1 = x1
        self.x2 = x2
        self.comm_model = model

    def build_mip_model(self):
        x_sh = []
        id1 = []
        id2 = []
        x_sh_dict = {}
        for key, val in self.medium.items():
            x = self.comm_model.add_var(lb=-val, ub=1000, name=key)
            x_sh +=[x]
            x_sh_dict[key] = x
            if key in self.S1_dict:
                id1 += [self.S1_dict[key]]
            else:
                id1 += [None]
            if key in self.S2_dict:
                id2 += [self.S2_dict[key]]
            else:
                id2 += [None]

        # Flux first model
        x1 = []
        for i, (lb, ub) in enumerate(self.bounds1): 
            x1 += [self.comm_model.add_var(lb = lb, ub=ub, name="x1" + str(i))]

        # Flux second model
        x2 = []
        for i, (lb, ub) in enumerate(self.bounds2): 
            x2 += [self.comm_model.add_var(lb = lb, ub=ub, name="x2" + str(i))]
        # Stoichiometry
        for i in range(self.S1.shape[0]):
            self.comm_model.add_constr(xsum(self.S1[i,j]*x1[j] for j in range(self.S1.shape[1])) == 0)

        for i in range(self.S2.shape[0]):
            self.comm_model.add_constr(xsum(self.S2[i,j]*x2[j] for j in range(self.S2.shape[1])) == 0)

        # Shuttel constraints
        for i,key in enumerate(self.medium):
            if id1[i] is not None and id2[i] is not None:
                idx1 = id1[i]
                idx2 = id2[i]
                self.comm_model.add_constr(-x_sh[i] + x1[idx1] + x2[idx2] == 0)
            elif id1[i] is not None:
                idx = id1[i]
                self.comm_model.add_constr(-x_sh[i] + x1[idx] == 0)
            else:
                idx = id2[i]
                self.comm_model.add_constr(-x_sh[i] + x2[idx] == 0)

        self.shuttel_reactions = x_sh_dict
        self.x1 = x1
        self.x2 = x2

    def compute_coopm(self, enforce_survival=True):
        minMBR = 0.1*self.optimize()
        # thetas
        thetas = []
        thetas_constraint = []
        for key,x in self.shuttel_reactions.items():
            V_min = -10
            if key == "EX_o2_e":
                V_min = -20
            if "_fe" in key:
                V_min = -0.1
            theta = self.comm_model.add_var(var_type=BINARY)
            thetas_constraint += [self.comm_model.add_constr(x + V_min*theta >= V_min)]
            thetas.append(theta)
        # Both must grow
        if enforce_survival:
            self.comm_model.add_constr(self.weights[0]*self.x1[self.obj1] >= minMBR)
            self.comm_model.add_constr(self.weights[1]*self.x2[self.obj2] >= minMBR)
        else:
            self.comm_model.add_constr(self.weights[0]*self.x1[self.obj1] + self.weights[1]*self.x2[self.obj2] >= minMBR)

        self.comm_model.objective = maximize(xsum(thetas))
        self.comm_model.optimize()

        coopm = dict()
        for key, x in self.shuttel_reactions.items():
            if x.x < 0:
                coopm[key] = abs(x.x)
        self.reset_model()
        return coopm
        

    def optimize(self):
        self.comm_model.objective = maximize(self.weights[0]*self.x1[self.obj1] + self.weights[1]*self.x2[self.obj2])
        self.comm_model.optimize()
        return self.weights[0]*self.x1[self.obj1].x + self.weights[1]*self.x2[self.obj2].x

    def optimize_single(self, model_idx):
        if model_idx == 0:
            self.comm_model.objective = maximize(self.x1[self.obj1])
            self.comm_model.optimize()
            return self.x1[self.obj1].x
        else:
            self.comm_model.objective = maximize(self.x2[self.obj2])
            self.comm_model.optimize()
            return self.x2[self.obj2].x

    def get_exchange_flux(self):
        dic1 ={}
        for key, val in self.S1_dict.items():
            if "EX_" in key:
                dic1[key] = self.x1[val].x
        dic2 ={}
        for key, val in self.S2_dict.items():
            if "EX_" in key:
                dic2[key] = self.x2[val].x
        return dic1, dic2

    def summary(self):
        self.optimize()
        ex1, ex2 = self.get_exchange_flux()
        shared_ex = list(set(ex1.keys()).intersection(set(ex2.keys())))
        interchange = {"DP flux":[], "SA flux":[]}
        index = shared_ex
        for key in index:
            interchange["DP flux"].append(ex1[key])
            interchange["SA flux"].append(ex2[key])
        df = pd.DataFrame(interchange)
        df.index = index 
        return df
             

    def set_medium(self, medium):
        for key in self.shuttel_reactions:
            if key in medium:
                self.shuttel_reactions[key].lb = -medium[key]
                self.medium[key] = -medium[key]
            else:
                self.shuttel_reactions[key].lb = 0.
                self.medium[key] = 0.
        
    def compute_alpha_weightings(self, alpha, maxbiomass=0.1):
        assert alpha <= 1 and alpha >= 0
        # Alpha objective...
        self.comm_model.add_constr(self.weights[0]*self.x1[self.obj1] <= alpha*maxbiomass)
        self.comm_model.add_constr(self.weights[1]*self.x2[self.obj2] <= (1-alpha)*maxbiomass)
        growth = self.optimize()
        summary = self.summary()
        self.reset_model()
        return growth, summary

    def compute_alpha_coopm(self, alpha, maxbiomass=0.1, enforce_survival=True):
        minMBR = maxbiomass*0.1
        thetas = []
        thetas_constraint = []
        for key,x in self.shuttel_reactions.items():
            V_min = -10
            if key == "EX_o2_e":
                V_min = -20
            if "_fe" in key:
                V_min = -0.1
            theta = self.comm_model.add_var(var_type=BINARY)
            thetas_constraint += [self.comm_model.add_constr(x + V_min*theta >= V_min)]
            thetas.append(theta)
        # Growth constraints
        self.comm_model.add_constr(self.weights[0]*self.x1[self.obj1] <= alpha*maxbiomass)
        self.comm_model.add_constr(self.weights[1]*self.x2[self.obj2] <= (1-alpha)*maxbiomass)
        # Both must grow
        if enforce_survival:
            self.comm_model.add_constr(self.weights[0]*self.x1[self.obj1] >= alpha*minMBR)
            self.comm_model.add_constr(self.weights[1]*self.x2[self.obj2] >= (1-alpha)*minMBR)
        else:
            self.comm_model.add_constr(self.weights[0]*self.x1[self.obj1] + self.weights[1]*self.x2[self.obj2] >= minMBR)

        self.comm_model.objective = maximize(xsum(thetas))
        self.comm_model.optimize()

        coopm = dict()
        for key, x in self.shuttel_reactions.items():
            if x.x < 0:
                coopm[key] = abs(x.x)
        self.reset_model()
        return coopm

def create_bag_of_react_model(models, biomass_functions, weights="auto"):
    # Bag of words model
    model = models[0] + models[1]

    if weights == "auto":
        biomass = []
        for model_i in models:
            biomass.append(model_i.slim_optimize())
        max_biomass = np.max(biomass)
        weights = [max_biomass/g for g in biomass]
    
    # Set objective as sum of growth
    objective = []
    for i, growth in enumerate(biomass_functions):
        v_i = weights[i] * model.reactions.get_by_id(growth).flux_expression
        objective.append(v_i)
    model.objective = sum(objective)
    return model

# Custum optimization

def optimize_coopm_community(model, MBR, biomass_reactions, weights, enforce_survival=True):
    model = model.copy()
    minMBR = 0.1*MBR
    medium = list(model.medium.keys())
    biomass = []
    for i,id in enumerate(biomass_reactions):
        biomass.append(weights[i]*model.reactions.get_by_id(id).flux_expression)
    # Binary variables: Theta_i 
    thetas = []
    for i in range(len(medium)):
        thetas.append(model.problem.Variable('theta_'+str(i), type="binary"))

    # Constraints for exchanges, which are turned of for theta_i = 1
    theta_constraints = []
    for i,id in enumerate(medium):
        reaction = model.reactions.get_by_id(id)
        min_bound = model.reactions.get_by_id(id).lower_bound
        cons = model.problem.Constraint(
            (reaction.flux_expression + min_bound*thetas[i]),
            lb=min_bound,
            ub=1000)
        theta_constraints.append(cons)

    # Constraints for growth rates, which must be at least 10% MBR
    if enforce_survival:
        constraint_growth = [model.problem.Constraint(
        f,
        lb=minMBR,
        ub=1000) for i,f in enumerate(biomass)]
    else:
        constraint_growth = model.problem.Constraint(
        sum(biomass),
        lb=minMBR,
        ub=1000)

    # Adding new variables and constraints.
    model.add_cons_vars(thetas)
    model.add_cons_vars(theta_constraints)
    model.add_cons_vars(constraint_growth)

    # Objevtive is maximising turned of exchanges, that is sum of theta_is
    objective = model.problem.Objective(sum(thetas), direction="max")
    model.objective = objective
    model.solver.update()

    sol = model.optimize()
    COOPM = dict()
    for id in medium:
        if sol.fluxes[id] < 0:
            COOPM[id] = abs(sol.fluxes[id])

    return COOPM
