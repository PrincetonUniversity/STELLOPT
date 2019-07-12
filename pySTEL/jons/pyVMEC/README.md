pyVMEC is a simple interface to the VMEC output files and simple analysis tools/methods for these files.

calc_jll_w7x.py -- load vmec results from Reference Equilibria webservice and plot the parallel current density
check_if_ID_exists.py -- boolean check wether a specified ID for a VMEC webservice run already exists
download_vmec_wout.py -- download the wout files for all successful VMEC runs on the webservice (netCDF)
gen_mgrid.py -- create an mgrid file on the VMEC webservice server for main field currents of a specified experiment
get_shafranov_shift_betavar.py -- plot the variation depending on beta of several equilibrium quantities for a W7-X reference equilibrium
read_vmec.py -- native Python interface to VMEC wout files (netCDF); also calculates currents (currumnc...)
README.md -- this file
test_pyVMEC.py -- check wether read_vmec reads the netCDF file contents correctly (with pySTEL from LIBSTELL as reference)
test_read_vmec.py -- test script for debugging the current calculation routines
test_read_vmec_against_libstell.py -- check wether read_vmec calculates the currents correctly (with pySTEL from LIBSTELL as reference)
test_vmec_webservice_currents.py -- check for all successful VMEC runs on the webservice wether the getFourierCoefficients method returns the currents correctly
w7x_currents.py -- fetch the main field coil currents of W7-X for either a specified experiment or a specified time interval; can do averaging/calculation of standard deviation
write_currents.py -- export the currents of VMEC reference runs to netCDF files



for python3: compile "netCDF4" package for python3
download: https://pypi.python.org/packages/36/42/69979a44e2b9a2eea3a083c1c955451d084a563ac6627a90c01c996c3756/netCDF4-1.2.8.tar.gz#md5=20e3d07c833fc5b7c29f401a6ee3159b
unpack, then do:
python3 setup.py build
sudo python3 setup.py install

History

 2017-06-18/jons created repository
 2017-08-16/jons moved code from MPCDF gitlab to this repository
