import sys
import imas,os

from vmec.actor import vmec
from vmec.common.runtime_settings import RunMode, DebugMode



class ExampleWorkflowManager:

    def __init__(self):

        self.vmec = vmec()
        self.input_entry = None
        self.output_entry = None

    def init_workflow(self):
        # INPUT/OUTPUT CONFIGURATION
        shot                = 131024
        run_in              = 1
        input_user_or_path  = 'public'
        input_database      = 'iter'
        run_out             = 10
        output_user_or_path = os.getenv('USER')
        output_database     = input_database

        # OPEN INPUT DATAFILE TO GET DATA FROM IMAS SCENARIO DATABASE
        print('=> Open input datafile')
        self.input_entry = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,input_database,shot,run_in,input_user_or_path)
        self.input_entry.open()

        # CREATE OUTPUT DATAFILE
        print('=> Create output datafile')
        self.output_entry = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,output_database,shot,run_out,output_user_or_path)
        self.output_entry.create()

        # # # # # # # # Initialization of ALL actors  # # # # # # # #
        code_parameters = self.vmec.get_code_parameters()
        code_parameters.parameters_path='./indata.xml'
        
        runtime_settings = self.vmec.get_runtime_settings()
        runtime_settings.run_mode = RunMode.STANDALONE
        runtime_settings.mpi.mpi_nodes = 1

        self.vmec.initialize(runtime_settings=runtime_settings, code_parameters=code_parameters)

    def execute_workflow(self):
        # READ INPUT IDSS FROM LOCAL DATABASE
        time_slice          = 200.
        print('=> Read input IDSs')
        input_equilibrium = self.input_entry.get_slice('equilibrium', time_slice, 1)
        input_equilibrium = imas.equilibrium()

        # EXECUTE PHYSICS CODE
        print('=> Execute physics code')
        output_equilibrium = self.vmec(input_equilibrium)

        # SAVE IDSS INTO OUTPUT FILE
        print('=> Export output IDSs to local database')
        self.output_entry.put(output_equilibrium)
        print('Done exporting.')

    def end_workflow(self):

        # Finalize ALL actors
        self.vmec.finalize()

        #other finalization actions
        self.input_entry.close()
        self.output_entry.close()

manager = ExampleWorkflowManager()

manager.init_workflow()
manager.execute_workflow()
manager.end_workflow()
