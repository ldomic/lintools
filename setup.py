
from setuptools import setup

setup(name='lintools',
      version='06.2017',
      description='Illustrates protein-ligand interactions',
      url='https://github.com/ldomic/lintools.git',
      author='Laura Domicevica',
      author_email='ldomicevica@gmail.com',
      license='GPL',
      packages=['lintools','lintools.analysis'],
      entry_points={
          'console_scripts': [
              'lintools = lintools.__main__:main'
          ]
      },
      zip_safe=False)
