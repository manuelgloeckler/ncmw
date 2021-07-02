

def check_transport_summary(model):
    compartments = list(model.compartments.keys())

    def check_transport(reaction):
        react = reaction.reaction
        result = False
        for comp1 in compartments:
            for comp2 in compartments:
                result = result or (("_" + comp1) in react and  ("_" + comp2) in react)
        return result

    transport_reactions = list(filter(check_transport, model.reactions))
    num_transport_reactions = len(transport_reactions)