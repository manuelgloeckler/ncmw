from ncmw.community.community_models import CommunityModel, create_stoichiometry_matrix, get_biomass_reaction
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
        self._open_exchanges_for_shuttles(-50)

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
        self._open_exchanges_for_shuttles(-50)

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

    def computeCOOPM(self, MBR, fraction=0.1, enforce_survival=False, n_tries=2):
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
