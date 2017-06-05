
from setuptools import setup

setup(name='lintools',
      version='05.2017',
      description='Illustrates protein-ligand interactions',
      url='http://github.com/ldomic/lintools',
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
