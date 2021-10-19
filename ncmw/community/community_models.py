import cobra
import numpy as np
import scipy

from mip import xsum, maximize, BINARY
import mip

from abc import ABC, abstractmethod, abstractproperty
import typing
from warnings import warn

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
        """This method returns the current objective value, if the objective only
        accounts a single model"""
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

    @staticmethod
    def load(path):
        with open(path, "rb") as f:
            return pickle.load(f)


class BagOfReactionsModel(CommunityModel):
    """This is a community model, which treats the community as a bag of all reactions that at least one the species has."""

    def __init__(self, models: list):
        self.community_model = cobra.Model(
            "Community__" + "".join([model.id + "__" for model in models])
        )
        self._type = "bag"
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
                for i in range(len(self._weights))
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
                model.problem.Constraint(f, lb=self._weights[i] * minMBR, ub=1000)
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

    def compute_convex_combination(self, alphas, maxbiomass=0.1):
        assert sum(alphas) == 1, "The weights must sum to one!"
        assert len(alphas) == len(self.models), "Scpecify a weight for each model..."
        model = self.community_model.copy()
        biomass = [
            model.reactions.get_by_id(f.id).flux_expression
            for f in self.biomass_reactions
        ]
        constraint_growth = [
            model.problem.Constraint(
                self.weights[i] * f,
                lb=alphas[i] * maxbiomass,
                ub=alphas[i] * maxbiomass,
            )
            for i, f in enumerate(biomass)
        ]
        model.add_cons_vars(constraint_growth)
        model.solver.update()

        sol = model.optimize()
        return model.summary()

    def save(self, path):
        """This saves the model using pickle. For other format overwrite this
        function"""
        self.biomass_reactions = None
        self.objective = None

        with open(path + ".pkl", "wb+") as f:
            pickle.dump(self, f)

        self.biomass_reactions = [get_biomass_reaction(model) for model in self.models]
        self.objective = sum([f.flux_expression for f in self.biomass_reactions])
        self.community_model.objective = self.objective

    @staticmethod
    def load(path):
        with open(path, "rb") as f:
            model = pickle.load(f)
        model.biomass_reactions = [
            get_biomass_reaction(model) for model in model.models
        ]
        model.objective = sum([f.flux_expression for f in model.biomass_reactions])
        model.community_model.objective = model.objective
        return model


class ShuttleCommunityModelMIP(CommunityModel):
    def __init__(self, models, shared_exchanges=None):
        # Set up model data
        self.models = models
        self._type = "compartmentalized"
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
        self._open_exchanges_for_shuttles(-100)

    def _set_weights(self, weights):
        self._weights = weights
        self.objective = xsum(
            [
                self._weights[i] * self.xs[i][self.biomass_ids[i]]
                for i in range(len(self.xs))
            ]
        )
        self.comm_model.objective = maximize(self.objective)

    def _set_medium(self, medium):
        for key in self.shuttle_reactions:
            if key in medium:
                self.shuttle_reactions[key].lb = -medium[key]
            else:
                self.shuttle_reactions[key].lb = 0.0
        self._medium = medium

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
                    ex_dict[key] = np.round(self.xs[j][val].x, 5)
            exchanges.append(ex_dict)

        return exchanges

    def _open_exchanges_for_shuttles(self, bound=-10):
        for j, rec_id in enumerate(self.rec_id_dicts):
            for key, val in rec_id.items():
                if key in self.shared_exchanges:
                    self.xs[j][val].lb = bound

    def summary(self):
        self.optimize()
        exchanges = self._get_exchange_flux()
        keys = []
        for dic in exchanges:
            keys.extend(list(dic.keys()))
        keys = list(set(keys))

        vals = np.zeros((len(keys), len(exchanges) + 1))
        for i in range(len(keys)):
            for j in range(len(exchanges) + 1):
                if j < len(exchanges) and keys[i] in exchanges[j]:
                    vals[i, j] = exchanges[j][keys[i]]
                if j >= len(exchanges):
                    if keys[i] in self.shuttle_reactions:
                        vals[i, j] = np.round(self.shuttle_reactions[keys[i]].x, 5)

        df = pd.DataFrame(
            vals, index=keys, columns=[m.id for m in self.models] + ["Total exchange"]
        )
        # Only return nonzero
        df = df[df["Total exchange"] != 0].sort_values(["Total exchange"])

        print("Objective: ", self.objective.x)
        for i in range(len(self.models)):
            print(
                self.models[i].id + " : ",
                self.xs[i][self.biomass_ids[i]].x,
                " with weights ",
                self._weights[i],
            )
        return df

    def _reset_model(self):
        model = mip.Model()
        model.read(self.path)
        for key in self.shared_exchanges:
            x = model.var_by_name(key)
            self.shuttle_reactions[key] = x
        for j, x in enumerate(self.xs):
            for i in range(len(x)):
                x_i = model.var_by_name(f"x{j}" + str(i))
                x[i] = x_i

        self.comm_model = model
        self._set_medium(self._medium)
        self._set_weights(self._weights)
        self.objective = xsum(
            [
                self._weights[i] * self.xs[i][self.biomass_ids[i]]
                for i in range(len(self.xs))
            ]
        )
        self.comm_model.objective = maximize(self.objective)
        self._open_exchanges_for_shuttles(-100)

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
                self._weights[i] * self.xs[i][self.biomass_ids[i]]
                for i in range(len(self.xs))
            ]
        )
        self.comm_model.objective = maximize(self.objective)

    def computeCOOPM(self, MBR, fraction=0.1, enforce_survival=True, n_tries=2):
        minMBR = fraction * MBR
        # thetas
        thetas = []
        thetas_constraint = []
        for key in self.medium:
            x = self.shuttle_reactions[key]
            V_min = -10.0
            if key == "EX_o2_e":
                V_min = -20.0
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
                    self.xs[i][self.biomass_ids[i]] >= self._weights[i] * minMBR
                )
        else:
            self.comm_model.add_constr(self.objective >= minMBR)

        self.comm_model.objective = maximize(xsum(thetas))

        status = self.comm_model.optimize()
        print("Optimization Status: ", status)

        if self.comm_model.objective.x is not None:
            coopm = dict()
            for key, x in self.shuttle_reactions.items():
                if x.x < 0:
                    coopm[key] = abs(x.x)
            self._reset_model()
            return coopm
        elif n_tries >= 0:
            new_fraction = 0.1 * fraction
            warn(
                f"The optimization failed, maybe some of the model is unable to achive fraction*MBR. We reduce the fraction to {new_fraction}"
            )
            self._reset_model()
            return self.computeCOOPM(
                MBR,
                new_fraction,
                enforce_survival=enforce_survival,
                n_tries=n_tries - 1,
            )
        else:
            coopm = dict()
            self._reset_model()
            raise ValueError(
                "The optimization failed, maybe some of the model is unable to achive fraction*MBR."
            )

    def compute_convex_combination(self, alphas, maxbiomass=0.1):
        assert sum(alphas) == 1, "The weights must sum to one!"
        assert len(alphas) == len(self.models), "Scpecify a weight for each model..."
        assert (
            max([self.single_optimize(i) for i in range(len(alphas))]) > maxbiomass
        ), "Each of the models must reach the maxbiomass..."
        # Alpha objective...
        for i in range(len(alphas)):
            self.comm_model.add_constr(
                self._weights[i] * self.xs[i][self.biomass_ids[i]]
                <= alphas[i] * maxbiomass
            )
        growth = self.optimize()
        summary = self.summary()
        self._reset_model()
        return summary

    def compute_convex_COOPM(
        self, alphas, maxbiomass=0.1, fraction=0.1, enforce_survival=True
    ):
        assert sum(alphas) == 1, "The weights must sum to one!"
        assert len(alphas) == len(self.models), "Scpecify a weight for each model..."
        minMBR = maxbiomass * fraction
        thetas = []
        thetas_constraint = []
        for key, x in self.shuttle_reactions.items():
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
                self.weights[i] * self.xs[i][self.biomass_ids[i]]
                == alphas[i] * maxbiomass
            )
        # # Both must grow
        # if enforce_survival:
        #     for i in range(len(self.models)):
        #         self.comm_model.add_constr(
        #             self.weights[i] * self.xs[i][self.biomass_ids[i]]
        #             >= alphas[i] * minMBR
        #         )
        # else:
        #     self.comm_model.add_constr(self.objective >= minMBR)

        self.comm_model.objective = maximize(xsum(thetas))
        self.comm_model.optimize()

        coopm = dict()
        for key, x in self.shuttle_reactions.items():
            if x.x < 0:
                coopm[key] = abs(x.x)
        self._reset_model()
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
        self.comm_model.write(path + "_help.lp")
        self.path = path + "_help.lp"
        self.comm_model = None
        for key in self.shuttle_reactions:
            self.shuttle_reactions[key] = None
        for i in range(len(self.xs)):
            for j in range(len(self.xs[i])):
                self.xs[i][j] = None
        self.objective = None
        super().save(path)
        self._reset_model()

    @staticmethod
    def load(path):
        with open(path, "rb") as f:
            model = pickle.load(f)
            model._reset_model()
            return model


class ShuttleCommunityModel(BagOfReactionsModel):
    def __init__(self, models, shared_exchanges=None):
        self.models = models
        # For this reactions we impose shuttle reactions if possible.
        if shared_exchanges is None:
            self.shared_exchanges = []
            for model in models:
                for ex in model.exchanges:
                    if ex.id not in self.shared_exchanges:
                        self.shared_exchanges.append(ex.id)
        self._type = "compartmentalized"

        self.build_community()

        self._weights = np.ones(len(self.models))
        self.biomass_reactions = []
        self.biomass_id = []
        for m in self.models:
            rec = get_biomass_reaction(m)
            id = rec.id + "__" + m.id
            self.biomass_id.append(id)
            self.biomass_reactions.append(self.community_model.reactions.get_by_id(id))
        self._set_weights(self._weights)

    def _set_weights(self, weights):
        self._weights = weights
        objective = []
        for i, r in enumerate(self.biomass_reactions):
            objective.append(self._weights[i] * r.flux_expression)

        self.community_model.objective = sum(objective)
        self.community_model.solver.update()

    def optimize(self):
        df = self.community_model.optimize()
        total_growth = df.objective_value
        single_growths = []
        for r in self.biomass_reactions:
            single_growths.append(r.flux)
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
        total_growth, single_growths = self.optimize()
        exchanges = self.community_model.exchanges
        non_zero_ids = []
        non_zero_flux = []
        for ex in exchanges:
            if ex.flux != 0:
                non_zero_ids.append(ex.id)
                non_zero_flux.append(ex.flux)

        vals = np.zeros((len(non_zero_ids), len(self.models) + 1))
        shuttle_reactions_ids = [s.id for s in self.shuttle_reactions]
        for i in range(len(non_zero_ids)):
            met = non_zero_ids[i][3:]
            for j in range(len(self.models)):
                shuttle_id = "SH_" + met + "__" + self.models[j].id
                if shuttle_id in shuttle_reactions_ids:
                    shuttle_rec = self.community_model.reactions.get_by_id(shuttle_id)
                    vals[i, j] = shuttle_rec.flux
        vals[:, -1] = np.array(non_zero_flux)

        df = pd.DataFrame(
            vals,
            index=non_zero_ids,
            columns=[m.id for m in self.models] + ["Total exchange"],
        ).sort_values(["Total exchange"])

        print("Objective: ", total_growth)
        for i in range(len(self.models)):
            print(
                self.models[i].id + " : ",
                single_growths[i],
                " with weights ",
                self._weights[i],
            )
        return df

    def build_community(self):
        # Community model we will build
        community_model = cobra.Model("Community")

        # Collect new metabolite contained in each model
        new_metabolites = []
        for model in self.models:
            for old_metabolite in model.metabolites:
                new_metabolite = cobra.Metabolite()
                for key, val in old_metabolite.__dict__.items():
                    if (
                        key != "_id"
                        and key != "_model"
                        and key != "_reaction"
                        and key != "compartment"
                    ):
                        new_metabolite.__dict__[key] = val
                new_metabolite.id = old_metabolite.id + "__" + model.id
                new_metabolite.compartment = (
                    old_metabolite.compartment + "__" + model.id
                )
                new_metabolites.append(new_metabolite)
        # Add new metabolites to community model
        community_model.add_metabolites(new_metabolites)

        # Collect new reactions contained in each model
        new_reactions = []
        for model in self.models:
            for old_reaction in model.reactions:
                # Skip exchange reaction, they are no longer part of the internal models
                if old_reaction.id[:3] == "EX_":
                    continue
                new_reaction = cobra.Reaction()
                for key, val in old_reaction.__dict__.items():
                    if (
                        key != "_id"
                        and key != "_model"
                        and key != "_metabolites"
                        and key != "compartment"
                    ):
                        new_reaction.__dict__[key] = val
                new_reaction.id = old_reaction.id + "__" + model.id
                for met, stoch in old_reaction.metabolites.items():
                    new_met_id = met.id + "__" + model.id
                    new_met = community_model.metabolites.get_by_id(new_met_id)
                    new_reaction.add_metabolites({new_met: stoch})
                new_reactions.append(new_reaction)
        # Add new reactions to community model
        community_model.add_reactions(new_reactions)

        # Add new exchange reactions
        external_metabolites = []
        for model in self.models:
            # The union of external metabolites
            for met in model.metabolites:
                if "_e" == met.id[-2:] and met.id not in external_metabolites:
                    external_metabolites.append(met.id)
        new_external_metabolites = []
        for id in external_metabolites:
            for model in self.models:
                if id in [met.id for met in model.metabolites]:
                    old_metabolite = model.metabolites.get_by_id(id)
                    break
            new_metabolite = cobra.Metabolite()
            for key, val in old_metabolite.__dict__.items():
                if (
                    key != "_id"
                    and key != "_model"
                    and key != "_reaction"
                    and key != "compartment"
                ):
                    new_metabolite.__dict__[key] = val
            new_metabolite.id = old_metabolite.id
            new_metabolite.compartment = "external"
            new_external_metabolites.append(new_metabolite)
        community_model.add_metabolites(new_external_metabolites)
        for met in new_external_metabolites:
            community_model.add_boundary(met)

        # Add shuttle reactions
        for model in self.models:
            for met in new_external_metabolites:
                if met.id in [met.id for met in model.metabolites]:
                    met2 = community_model.metabolites.get_by_id(
                        met.id + "__" + model.id
                    )
                else:
                    continue
                shuttle_reaction = cobra.Reaction()
                id = met.id
                shuttle_reaction.id = f"SH_{id}__{model.id}"
                shuttle_reaction.name = f"Shuttle reaction for {id}"
                shuttle_reaction.lower_bound = -100
                shuttle_reaction.upper_bound = 1000
                shuttle_reaction.add_metabolites({met: 1, met2: -1})
                community_model.add_reaction(shuttle_reaction)
        # Add shuttle constraints!
        all_shuttles = [ex.id for ex in community_model.reactions if "SH_" == ex.id[:3]]
        shuttle_constraints = []
        for met in new_external_metabolites:
            present_shuttles = []
            for model in self.models:
                id = f"SH_{met.id}__{model.id}"
                if id in all_shuttles:
                    present_shuttles.append(
                        community_model.reactions.get_by_id(id).flux_expression
                    )
            exchange = community_model.reactions.get_by_id(
                f"EX_{met.id}"
            ).flux_expression
            cons = community_model.problem.Constraint(
                exchange - sum(present_shuttles),
                lb=0,
                ub=1000,
            )
            shuttle_constraints.append(cons)
        community_model.add_cons_vars(shuttle_constraints)

        # Set correct medium
        medium = {}
        for model in self.models:
            for key, val in model.medium.items():
                if key not in medium:
                    medium[key] = val

        self.community_model = community_model
        self.shuttle_reactions = [
            ex for ex in community_model.reactions if "SH_" == ex.id[:3]
        ]
        self._set_medium(medium)

    def save(self, path):
        """This saves the model using pickle. For other format overwrite this
        function"""
        self.biomass_reactions = None
        self.objective = None

        with open(path + ".pkl", "wb+") as f:
            pickle.dump(self, f)

        self.biomass_ids = [
            get_biomass_reaction(model).id + "__" + model.id for model in self.models
        ]
        self.biomass_reactions = [
            self.community_model.reactions.get_by_id(i) for i in self.biomass_ids
        ]
        self.objective = sum([f.flux_expression for f in self.biomass_reactions])
        self.community_model.objective = self.objective

    def load(path):
        with open(path, "rb") as f:
            model = pickle.load(f)
        model.biomass_ids = [
            get_biomass_reaction(model).id + "__" + model.id for model in model.models
        ]
        model.biomass_reactions = [
            model.community_model.reactions.get_by_id(i) for i in model.biomass_ids
        ]
        model.objective = sum([f.flux_expression for f in model.biomass_reactions])
        model.community_model.objective = model.objective
        return model
