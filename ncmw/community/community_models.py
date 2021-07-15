import cobra
import numpy as np
import scipy

from mip import xsum, maximize, BINARY
import mip

from abc import ABC, abstractmethod, abstractproperty
import typing

import sys, os
import tempfile

import pickle

import pandas as pd


def create_stoichiometry_matrix(model):
    """This creates a stoichiometry matrix"""
    metabolites = model.metabolites
    reactions = model.reactions
    S = np.zeros((len(metabolites), len(reactions)))

    met_id = dict()
    rec_id = dict()
    for i, reaction in enumerate(model.reactions):
        rec_id[reaction.id] = i
        for metabolite, stoich in reaction.metabolites.items():
            met_id[metabolite.id] = int(metabolites.index(metabolite))
            S[metabolites.index(metabolite), i] = stoich
    # Compress as sparse matrix
    S = scipy.sparse.csr_matrix(S)
    return S, met_id, rec_id


def get_biomass_reaction(model):
    """Return the biomass reaction of a model"""
    objective_str = str(list(model.objective.variables)[0])
    for rec in model.reactions:
        if rec.id in objective_str:
            return rec


class CommunityModel(ABC):
    """This gives the backbone of any community model."""

    @property
    def medium(self):
        return self._medium

    @medium.setter
    def medium(self, medium):
        return self._set_medium(medium)

    @abstractmethod
    def _set_medium(self, medium):
        """This method should set the medium properly"""
        pass

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, weights):
        return self._set_weights(weights)

    @abstractmethod
    def _set_weights(self, weights):
        """This method should set the weights properly"""
        pass

    @abstractmethod
    def slim_optimize(self):
        """This method returns the current objective value"""
        pass

    @abstractmethod
    def single_optimize(self, idx):
        """This method returns the current objective value"""
        pass

    @abstractmethod
    def optimize(self):
        """This method return the current objective value and additional information"""
        pass

    def computeCOOPM(self, MBR, fraction=0.1, enforce_survival=True):
        """This method computes the COOPM medium"""
        raise NotImplementedError(
            "This method is not implemented fro your current model"
        )

    def summary(self):
        """This method should report a summary of the model"""
        raise NotImplementedError(
            "This method is not implemented for your current model"
        )

    def save(self, path):
        """This saves the model using pickle. For other format overwrite this function"""
        with open(path + ".pkl", "wb+") as f:
            pickle.dump(self, f)

    def load(self, path):
        with open(path + ".pkl", "wb") as f:
            self = pickle.load(f)


class BagOfReactionsModel(CommunityModel):
    """This is a community model, which treats the community as a bag of all reactions that at least one the species has."""

    def __init__(self, models: list):
        self.community_model = cobra.Model(
            "Community__" + "".join([model.id + "__" for model in models])
        )
        for model in models:
            self.community_model += model
        self.models = models
        self._weights = np.ones(len(models))
        self._medium = self.community_model.medium
        self.biomass_reactions = [get_biomass_reaction(model) for model in models]
        self.objective = sum([f.flux_expression for f in self.biomass_reactions])
        self.community_model.objective = self.objective

    def _set_medium(self, medium):
        self.community_model.medium = medium
        self._medium = medium

    def _set_weights(self, weights):
        assert len(weights) == len(
            self.models
        ), "You need to specify for each species in the community a weights..."
        self._weights = weights
        self.objective = sum(
            [
                weights[i] * self.biomass_reactions[i].flux_expression
                for i in range(len(self.weights))
            ]
        )
        self.community_model.objective = self.objective

    def slim_optimize(self):
        return self.community_model.slim_optimize()

    def optimize(self):
        sol = self.community_model.optimize()
        total_growth = self.community_model.slim_optimize()
        single_growths = [sol[r.id] for r in self.biomass_reactions]
        return total_growth, single_growths

    def single_optimize(self, idx):
        weights = np.zeros(len(self.models))
        weights[idx] = 1.0
        old_weights = self._weights
        self._set_weights(weights)
        growth = self.slim_optimize()
        self._set_weights(old_weights)

        return growth

    def summary(self):
        return self.community_model.summary()

    def computeCOOPM(self, MBR, fraction=0.1, enforce_survival=True):
        model = self.community_model.copy()
        minMBR = fraction * MBR
        medium = list(model.medium.keys())
        biomass = [
            model.reactions.get_by_id(f.id).flux_expression
            for f in self.biomass_reactions
        ]

        # Binary variables: Theta_i
        thetas = []
        for i in range(len(medium)):
            thetas.append(model.problem.Variable("theta_" + str(i), type="binary"))

        # Constraints for exchanges, which are turned of for theta_i = 1
        theta_constraints = []
        for i, id in enumerate(medium):
            reaction = model.reactions.get_by_id(id)
            min_bound = model.reactions.get_by_id(id).lower_bound
            cons = model.problem.Constraint(
                (reaction.flux_expression + min_bound * thetas[i]),
                lb=min_bound,
                ub=1000,
            )
            theta_constraints.append(cons)

        # Constraints for growth rates, which must be at least 10% MBR
        if enforce_survival:
            constraint_growth = [
                model.problem.Constraint(f, lb=self.weights[i] * minMBR, ub=1000)
                for i, f in enumerate(biomass)
            ]
        else:
            constraint_growth = model.problem.Constraint(
                sum(biomass), lb=minMBR, ub=1000
            )

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


class ShuttleCommunityModel(CommunityModel):
    def __init__(self, models, shared_exchanges=None):
        # Set up model data
        self.models = models
        self.stoichiometry_matrixes = []
        self.rec_id_dicts = []
        self.met_id_dicts = []
        self._medium = []
        for model in models:
            S, met_id, rec_id = create_stoichiometry_matrix(model)
            self.stoichiometry_matrixes.append(S)
            self.rec_id_dicts.append(rec_id)
            self.met_id_dicts.append(met_id)
            self._medium += list(model.medium.items())
        self._medium = dict(self._medium)
        self._bounds = self._get_bounds()

        # Set shuttle reactions to all exchanges if not specified
        if shared_exchanges is None:
            self.shared_exchanges = []
            for model in models:
                for ex in model.exchanges:
                    if ex.id not in self.shared_exchanges:
                        self.shared_exchanges.append(ex.id)

        self._weights = np.ones(len(models))
        self.biomass_reactions = [get_biomass_reaction(model) for model in models]
        self.biomass_ids = [
            self.rec_id_dicts[i][self.biomass_reactions[i].id]
            for i in range(len(models))
        ]

        self.comm_model = mip.Model("Community Model")
        self.build_mip_model()

    def _set_weights(self, weights):
        self._weights = weights
        self.objective = xsum(
            [
                self.weights[i] * self.xs[i][self.biomass_ids[i]]
                for i in range(len(self.xs))
            ]
        )
        self.comm_model.objective = maximize(self.objective)

    def _set_medium(self, medium):
        self._medium = medium
        for key in self.shuttel_reactions:
            if key in medium:
                self.shuttel_reactions[key].lb = -medium[key]
                self.medium[key] = -medium[key]
            else:
                self.shuttel_reactions[key].lb = 0.0
                self.medium[key] = 0.0

    def optimize(self):
        self.comm_model.optimize()
        total_growth = self.objective.x
        single_growths = []
        for x, id in zip(self.xs, self.biomass_ids):
            single_growths.append(x[id].x)
        return total_growth, single_growths

    def single_optimize(self, idx):
        weights = np.zeros(len(self.models))
        weights[idx] = 1.0
        old_weights = self._weights
        self._set_weights(weights)
        growth = self.slim_optimize()
        self._set_weights(old_weights)
        return growth

    def slim_optimize(self):
        self.comm_model.optimize()
        return self.objective.x

    def _get_exchange_flux(self):
        exchanges = []
        for j, rec_id in enumerate(self.rec_id_dicts):
            ex_dict = {}
            for key, val in rec_id.items():
                if "EX_" in key:
                    ex_dict[key] = self.xs[j][val].x
            exchanges.append(ex_dict)

        return exchanges

    def summary(self):
        self.optimize()
        exchanges = self._get_exchange_flux()
        shared_ex = set(exchanges[0].keys()).intersection(
            *[set(ex.keys()) for ex in exchanges[1:]]
        )
        titles = [model.id + " flux" for model in self.models]
        columns = [[] for _ in range(len(self.models))]
        interchange = dict(zip(titles, columns))
        index = shared_ex
        for key in index:
            for i in range(len(titles)):
                interchange[titles[i]].append(exchanges[i][key])
        df = pd.DataFrame(interchange)
        df.index = index
        print("Objective: ", self.objective.x)
        for i in range(len(self.models)):
            print(
                self.models[i].id + " : ",
                self.xs[i][self.biomass_ids[i]].x,
                " with weights ",
                self.weights[i],
            )
        return df

    def _reset_model(self):
        model = mip.Model()
        model.read(self.path)
        for key in self.shared_exchanges:
            x = model.var_by_name(key)
            self.shuttle_reactions[key] = x
        for j, x in enumerate(self.xs):
            for i in range(len(self._bounds[i])):
                x_i = model.var_by_name(f"x{j}" + str(i))
                x[i] += [x_i]

        self.comm_model = model

    def build_mip_model(self):
        ids = [[] for _ in range(len(self.rec_id_dicts))]
        x_sh = []
        x_sh_dict = {}
        for key in self.shared_exchanges:
            if key in self._medium:
                val = self._medium[key]
            else:
                val = 0
            x = self.comm_model.add_var(lb=-val, ub=1000, name=key)
            x_sh_dict[key] = x
            x_sh += [x]
            for j, rec_id in enumerate(self.rec_id_dicts):
                if key in rec_id:
                    ids[j] += [rec_id[key]]
                else:
                    ids[j] += [None]

        # Flux first model
        xs = []
        for j, bounds in enumerate(self._bounds):
            xs_j = []
            for i, (lb, ub) in enumerate(bounds):
                xs_j += [self.comm_model.add_var(lb=lb, ub=ub, name=f"x{j}" + str(i))]
            xs.append(xs_j)

        # Stoichiometry
        for k, S in enumerate(self.stoichiometry_matrixes):
            for i in range(S.shape[0]):
                self.comm_model.add_constr(
                    xsum(S[i, j] * xs[k][j] for j in range(S.shape[1])) == 0
                )

        # Shuttel constraints
        for i, key in enumerate(self.shared_exchanges):
            exchanges = [
                xs[k][ids[k][i]] for k in range(len(ids)) if ids[k][i] is not None
            ]
            self.comm_model.add_constr(-x_sh[i] + xsum(exchanges) == 0)

        fd, path = tempfile.mkstemp()
        path = path + ".lp"
        self.comm_model.write(path)

        self.path = path
        self.shuttle_reactions = x_sh_dict
        self.xs = xs
        self.objective = xsum(
            [
                self.weights[i] * self.xs[i][self.biomass_ids[i]]
                for i in range(len(self.xs))
            ]
        )
        self.comm_model.objective = maximize(self.objective)

    def computeCOOPM(self, MBR, fraction=0.1, enforce_survival=True):
        minMBR = fraction * MBR
        # thetas
        thetas = []
        thetas_constraint = []
        for key, x in self.shuttel_reactions.items():
            # TODO SET THIS TO THE MEDIUM??? Thats actually wrong...
            V_min = -10
            if key == "EX_o2_e":
                V_min = -20
            if "_fe" in key:
                V_min = -0.1
            theta = self.comm_model.add_var(var_type=BINARY)
            thetas_constraint += [
                self.comm_model.add_constr(x + V_min * theta >= V_min)
            ]
            thetas.append(theta)
        # Both must grow
        if enforce_survival:
            for i in range(len(self.models)):
                self.comm_model.add_constr(
                    self.xs[i][self.biomass_ids[i]] >= self.weights[i] * minMBR
                )
        else:
            self.comm_model.add_constr(self.objective >= minMBR)

        self.comm_model.objective = maximize(xsum(thetas))
        self.comm_model.optimize()

        coopm = dict()
        for key, x in self.shuttel_reactions.items():
            if x.x < 0:
                coopm[key] = abs(x.x)
        self.reset_model()
        return coopm

    def compute_convex_combination(self, alphas, maxbiomass=0.1):
        assert sum(alphas) == 1, "The weights must sum to one!"
        assert len(alphas) == len(self.models), "Scpecify a weight for each model..."
        assert (
            max([self.single_optimize(i) for i in range(len(alphas))]) > maxbiomass
        ), "Each of the models must reach the maxbiomass..."
        # Alpha objective...
        for i in range(len(alphas)):
            self.comm_model.add_constr(
                self.weights[i] * self.x1[self.obj1] <= alphas[i] * maxbiomass
            )
        growth = self.optimize()
        summary = self.summary()
        self.reset_model()
        return growth, summary

    def compute_convex_COOPM(
        self, alphas, maxbiomass=0.1, fraction=0.1, enforce_survival=True
    ):
        assert sum(alphas) == 1, "The weights must sum to one!"
        assert len(alphas) == len(self.models), "Scpecify a weight for each model..."
        minMBR = maxbiomass * fraction
        thetas = []
        thetas_constraint = []
        for key, x in self.shuttel_reactions.items():
            V_min = -10
            if key == "EX_o2_e":
                V_min = -20
            if "_fe" in key:
                V_min = -0.1
            theta = self.comm_model.add_var(var_type=BINARY)
            thetas_constraint += [
                self.comm_model.add_constr(x + V_min * theta >= V_min)
            ]
            thetas.append(theta)
        # Growth constraints
        for i in range(len(alphas)):
            self.comm_model.add_constr(
                self.weights[i] * self.x1[self.obj1] <= alphas[i] * maxbiomass
            )
        # Both must grow
        if enforce_survival:
            for i in range(len(self.models)):
                self.comm_model.add_constr(
                    self.weights[i] * self.xs[i][self.biomass_ids[i]]
                    >= alphas[i] * minMBR
                )
        else:
            self.comm_model.add_constr(self.objective >= minMBR)

        self.comm_model.objective = maximize(xsum(thetas))
        self.comm_model.optimize()

        coopm = dict()
        for key, x in self.shuttel_reactions.items():
            if x.x < 0:
                coopm[key] = abs(x.x)
        self.reset_model()
        return coopm

    def _get_bounds(self):
        all_bounds = []
        for model in self.models:
            bounds = []
            for rec in model.reactions:
                bounds.append((rec.lower_bound, rec.upper_bound))
            all_bounds.append(bounds)
        return all_bounds

    def save(self, path):
        """This saves the model using pickle. For other format overwrite this function"""
        self.comm_model.write(path + ".lp")

    def load(self, path):
        self.path = path + ".lp"
        self._reset_model()
