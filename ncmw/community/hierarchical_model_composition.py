from libsbml import *

import libsbml
import os
import cobra

class HierarchicalModel():
    """Hierarchical community model

    This class uses LibSBML and the comp package for SBML to generate a community model
    from single SBML models!
    
    NOTE: To use it for optimization it must be "flattened", i.e. the hierarchy is
    collapsed and a merged model is generated. The "flattended" model can then be used
    e.g. within COBRA!
    
    """
    def __init__(self, name="Community"):
        sbmlns = SBMLNamespaces(3,1)
        sbmlns.addPackageNamespace("comp",1)
        sbmlns.addPackageNamespace("fbc",2)
        document = SBMLDocument(sbmlns)
        document.setPackageRequired("comp", True)
        document.setPackageRequired("fbc", False)

        self.sbml_doc = document 
        self.model = self.sbml_doc.createModel(name)
        self.model.setName(name)
        self.model.getPlugin("fbc").setStrict(True)
        self.models = []
        self.models_ref = []

    def add_model(self,source:str, id:str=None):
        """Does add a new model to the community
        
        Args:
            id: Identifier used within the community
            source: Path to the corresponding file    
        """
        reader = SBMLReader()
        doc = reader.readSBML(source)
        self.models.append(doc.getModel())

        if id is None:
            id = self.models[-1].getId()
        else:
            self.models[-1].setId(id)

        self.models_ref.append(id)

        m1 = ExternalModelDefinition(3,1)
        m1.setId(id)
        m1.setSource(source)
        dplugin = self.sbml_doc.getPlugin("comp")
        dplugin.addExternalModelDefinition(m1)

        smodel = Submodel()
        smodel.setModelRef(id)
        smodel.setId(id)
        mplugin = self.model.getPlugin("comp")
        mplugin.addSubmodel(smodel)

    def getMedium(self):
        medium = dict()
        for m in self.models:
            parameters = m.getListOfParameters()
            for p in parameters:
                id = p.getId().replace("R_", "").replace("_lower_bound", "")
                if "EX_" in id and not "upper_bound" in id and id not in medium:
                    val = abs(p.getValue())
                    if val != 0:
                        medium[id] = abs(p.getValue())
        return medium


    def getIds(self):
        """ Returns the model ids"""
        ids = []
        for m in self.models:
            ids.append(m.getId())
        return ids

    def getSubmodel(self, idx):
        """Retruns the submodel given an index. The order is given by the order the
        models were added to the community!
          
        Args:
            idx: Integer identifier
        
        Returns:
            object: Submodel
        
        """
        return self.model.getPlugin("comp").getSubmodel(idx)

    def getModel(self,idx):
        """ Returns the original model"""
        return self.models[idx]

    def getReactionsFromModel(self, idx):
        """ Returns the reactions form the original model"""
        return [r.getId() for r in self.models[idx].getListOfReactions()]

    def getMetabolitesFromModel(self, idx):
        """ Returns the metabolites from the original model"""
        return [e.getId() for e in self.models[idx].getListOfSpecies()]

    def getCompartmentsFromModel(self,idx):
        """Retruns the compartments from the original model""" 
        return [c.getId() for c in self.models[idx].getListOfCompartments()]

    def getExchangesFromModel(self, idx):
        """ Returns all exchange reaction from the original model"""
        rec = self.getReactionsFromModel(idx)
        i = 0
        while i < len(rec):
            id = rec[i]
            if not ("EX_" in id and "_e" in id):
                del rec[i]
            else:
                i += 1
        return rec

    def getExternalMetabolites(self,idx):
        """Retruns all external metabolites from the original model""" 
        met = self.getMetabolitesFromModel(idx)
        i = 0
        while i < len(met):
            id = met[i]
            if not (id[-2:] == "_e"):
                del met[i]
            else:
                i += 1
        return met

    def getAllExchanges(self):
        """ Returns the union of all exchange reactions"""
        ex = []
        for idx in range(len(self.models)):
            ex.extend(self.getExchangesFromModel(idx))
        return list(set(ex))
    
    def getAllExternalMetabolites(self):
        """ Retruns the union of all external metabolites"""
        ex = []
        for idx in range(len(self.models)):
            ex.extend(self.getExternalMetabolites(idx))
        return list(set(ex))

    def getAllBiomassFunctons(self):
        """ Returns all biomass functions.""" 
        biomass_functions = []
        for ids, m in zip(self.models_ref, self.models):
            mfbc = m.getPlugin("fbc")
            o = mfbc.getObjective(0)
            reaction_id = o.getFluxObjective(0).getReaction()
            biomass_functions.append(ids + "__" + reaction_id)
        return biomass_functions

    def getModelRefContainMetabolites(self, met:str):
        refs = []
        for i in range(len(self.models)):
            if met in self.getMetabolitesFromModel(i):
                refs.append(self.models_ref[i])
        return refs

    def collapse_compartments(self, name:str, model_ref:str, compartments:list, **kwargs):
        """ Collapses the compartments to a single one!"""
        c = self.add_compartment(name, **kwargs)
        for comp in compartments:
            self.replace(c, model_ref, comp)
        return c

    def removeAllExchangesFromModel(self, idx):
        """Reomoves old exchange reactions """ 
        s1 = self.getSubmodel(idx)
        exchanges = self.getExchangesFromModel(idx)
        for e in exchanges:
            self.delete(s1, e)

    def addExternalEnvironment(self, shuttle_lower_bound=-50, shuttle_upper_bound=50):
        """ Adds an new external environment with new exchange reactions """
        _ = self.add_compartment("e", "external")
        mets = self.getAllExternalMetabolites()
        for m in mets:
            self.add_metabolite(m, "e")
            within_models = self.getModelRefContainMetabolites(m)
            _ = self.add_reaction("EX_" + m[2:], {m:1}, dict())
            for model in within_models:
                _ = self.add_reaction(f"SH__{model}_"+m, {m:1}, {f"{model}__{m}":1}, lower_bound=shuttle_lower_bound, upper_bound=shuttle_upper_bound)


    def delete(self,submodel, id_ref:str):
        """Deletes the corresponding id_ref from a given submodel

        NOTE: The element is deleted within the community not in the original model itself
        
        Args:
            submodel: Submodel object
            id_ref: Identifier to delete e.g. a reaction id 
        
        
        """
        deleted = Deletion()
        deleted.setIdRef(id_ref)
        submodel.addDeletion(deleted)

    def replace(self, element, model_ref:str, id_ref:str):
        """Replaces a element of the original one with a substitute within the community
    
        Args:
            element: Element to replace
            model_ref: Model reference for which we should replace 
            id_ref: Identifier of the object to replace
        """
        replaced = ReplacedElement()
        replaced.setIdRef(id_ref)
        replaced.setSubmodelRef(model_ref)
        element.getPlugin("comp").addReplacedElement(replaced)

    def setCommunityObjective(self, coefficients, biomass_reactions):
        """Community objectie of sum of individual growths is set as new objective
        
        Args:
            coefficients: Weights for each growth.
            biomass_reactions: Individual biomass reactions.
        
        
        """
        mfbc = self.model.getPlugin("fbc")
        objective = mfbc.createObjective()
        objective.setId("Growth")
        objective.setType("maximize")
        mfbc = self.model.getPlugin("fbc")
        mfbc.setActiveObjectiveId("Growth")

        for c, b in zip(coefficients, biomass_reactions):
            fluxObjective = objective.createFluxObjective()
            fluxObjective.setReaction(b)
            fluxObjective.setCoefficient(c)

    def add_compartment(self,id,name="", meta_id="", constant=True, size=1):
        """Adds a new compartment to the community
        
        Args:
            id: Identifier of the compartment
            name: Name of the compartment
            meta_id: Meta identfier
            constant: If it is ocnstant
            size: Size of the compartment
        
        Returns:
            object: Compartment
        
        """
        comp = self.model.createCompartment()
        comp.setId(id)
        comp.setName(name)
        comp.setMetaId(meta_id)
        comp.setConstant(constant)
        comp.setSize(size)
        return comp

    def add_parameter(self, id:str, value, constant=True):
        """A name parameter e.g. can be used to save the default medium within the model.
        
        Args:
            id: Identifer.
            value: Numerical value fo the parameter.
            constant: If it is constant.
        
        
        """
        p = Parameter(3,1)
        p.setId(id)
        p.setValue(value)
        p.setConstant(constant)
        self.model.addParameter(p)

    def add_metabolite(self, id:str, compartment:str, boundary:bool=False, substance_units:bool=True, constant:bool=True):
        """Adds a metabolites to the community
        
        Args:
            id: Identifier
            compartment: Compartment to which the metabolite is added
            boundary: If it is a bondary reaction
            substance_units: Units
            constant: If it is constant
        
        Returns:
            object: metabolite
        
        """
        species = self.model.createSpecies()
        species.setId(id)
        species.setCompartment(compartment)
        species.setHasOnlySubstanceUnits(substance_units)
        species.setBoundaryCondition(boundary)
        species.setConstant(constant)
        return species

    def add_reaction(self, id:str, reactants:dict, products:dict, lower_bound=-10, upper_bound=1000,reversible=True, constant=False):
        """Adds a reaction to the community
        
        Args:
            id: Identifer
            reactants: Dictionary of reactants
            products: Dictionary of products
            lower_bound: Lower flux bound
            upper_bound: Upper flux bound
            reversible: Reversibility
            constant: Constant?
        
        Returns:
            object: Reaction
        
        """
        if not isinstance(reactants, dict):
            raise ValueError()
        if not isinstance(products, dict):
            raise ValueError()
        reaction = self.model.createReaction()
        reaction.setId(id)
        reaction.setReversible(reversible)
        for key,val in reactants.items():
            reactant = reaction.createReactant()
            reactant.setSpecies(key)
            reactant.setStoichiometry(val)
            reactant.setConstant(constant)
        for key,val in products.items():
            product = reaction.createProduct()
            product.setSpecies(key)
            product.setStoichiometry(val)
            product.setConstant(constant)

        self.add_parameter(id+ "_lb", lower_bound)
        self.add_parameter(id+ "_ub", upper_bound)

        rplug = reaction.getPlugin("fbc")
        rplug.setLowerFluxBound(id+"_lb")
        rplug.setUpperFluxBound(id + "_ub")
        
        return reaction

    def convert_to_cobra(self, path="test.xml", delete_file_after=True):
        """ Converts the model to a cobra model. """
        self.save_flatten(path)
        model = cobra.io.read_sbml_model(path)
        if delete_file_after:
            os.remove(path)
        return model

    def reset(self):
        """ Resets the model, clears everythink"""
        self.__init__(self.model.getName())

    def save_flatten(self, path):
        """Runs the flatting algorithm to produce a 'cobra' readable flattened community
        model. And saves it as SBML.
        
        Args:
            path: Path where the file is saved
        
        
        """
        sbmldoc = self.sbml_doc

        props = ConversionProperties()
        props.addOption('flatten comp', True) 
        converter = libsbml.CompFlatteningConverter()
        converter.setDocument(sbmldoc)

        status = converter.convert()
        if status != 0:
            raise ValueError("Not a valid hierarchical model")

        writer  = SBMLWriter()
        writer.writeSBML(sbmldoc, path)        

    def save(self, path):
        writeSBMLToFile(self.sbml_doc, path)


def generate_hierarchical_community_model(model_paths, weights=None):
    """Generates a hierarchical_community_model writer
    
    Args:
        model_paths: Paths to the xml files
        model_ids: Model identifiers
    """

    model = HierarchicalModel()

    if weights is None:
        weights = [1 for _ in range(len(model_paths))]

    i = 0
    for path in model_paths:
        model.add_model(path)
        id = model.models[-1].getId()
        model.collapse_compartments(id + "_c", id, model.getCompartmentsFromModel(i))
        model.removeAllExchangesFromModel(i)
        i += 1
    model.addExternalEnvironment()
    model.setCommunityObjective(weights, model.getAllBiomassFunctons())

    return model
        