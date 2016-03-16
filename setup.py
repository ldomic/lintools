from setuptools import setup

setup(name='lintools',
      version='0.1-beta',
      description='Illustrates protein-ligand interactions',
      url='http://github.com/ldomic/lintools',
      download_url='http://github.com/ldomic/lintools/tarball/0.1b0',
      author='Laura Domicevica',
      author_email='ldomicevica@gmail.com',
      license='GPL',
      packages=['lintools'],
      install_requires=[
	     'MDAnalysis','openbabel','shapely'
      ],
      zip_safe=False)
