
def get_simulation(experiment, simulation_id):
    ''' If the simulation id is known, get the simulation objects. '''
    
    # Initialize the experiment objects
    simulations = None
    
    # Get the simulation object that matches the id 
    for simulation in experiment.simulations:
        if simulation.id == simulation_id:
            simulations = [simulation] 
            
    # If no id matched, get the first simulation
    if simulations==None and simulation_id!="All simulations":
        simulations = [experiment.simulations[0]]
    
    # Unless the id is "All experiments", then add them all
    if simulations==None and simulation_id=="All simulations":
        simulations = experiment.simulations
        
    # Return the experiment objects
    return simulations
