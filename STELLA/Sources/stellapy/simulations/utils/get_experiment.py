
def get_experiment(research, experiment_id):
    ''' If the experiment id is known, get the experiment objects. '''
    
    # Initialize the experiment objects
    experiments = None
    
    # Get the experiment object that matches the id 
    for experiment in research.experiments:
        if experiment.id == experiment_id:
            experiments = [experiment] 
            
    # If no id matched, get the first experiment
    if experiments==None and experiment_id!="All experiments":
        experiments = [research.experiments[0]]
    
    # Unless the id is "All experiments", then add them all
    if experiments==None and experiment_id=="All experiments":
        experiments = research.experiments
        
    # Return the experiment objects
    return experiments
