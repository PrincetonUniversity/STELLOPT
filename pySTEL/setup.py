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
	scripts = ['VMECplot.py','FIELDLINESplot.py','vmec_util.py',\
		'boozer_util.py','coils_util.py','fieldlines_util.py'],
	install_requires=['numpy','matplotlib','PyQt5','scipy', \
		'contourpy','vtk','numpy-stl','gmsh']
	)
