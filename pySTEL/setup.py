#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='pySTEL',
	version = '1.0.0',
	description = 'Python library for interfacing with STELLOPT',
	long_description =	'This software package contains python '+ \
						'software for interfacing with the STELLOPT'+\
						'package of codes.',
	author = 'Samuel A. Lazerson',
	author_email = 'lazersos@gmail.com',
	url = 'https://github.com/PrincetonUniversity/STELLOPT',
	packages=['libstell'],
	scripts = ['VMECplot.py','FIELDLINESplot.py', 'bootsj_util.py', \
		'boozer_util.py','coils_util.py','fieldlines_util.py',\
		'focus_util.py','gist_util.py','make_mesh.py','nescoil_util.py',\
		'STELLOPT.py','stellopt_renorm.py','vmec2beams3d.py','vmec2focus.py',\
		'vmec_util.py','wall_util.py'],
	install_requires=['numpy','matplotlib','PyQt5','scipy', \
		'contourpy','PyVTK','numpy-stl','gmsh','vtk']
	)
