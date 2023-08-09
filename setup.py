import setuptools

def read_requirements(fname):
    with open(fname, 'r', encoding='utf-8') as file:
        return [line.rstrip() for line in file]


setuptools.setup(
     name='pycistarget',
     version="1.0a1",
     packages=setuptools.find_packages(where='src'),
     package_dir={'': 'src'},
     install_requires=read_requirements('requirements.txt'),
     author="Carmen Bravo, Seppe de Winter",
     author_email="carmen.bravogonzalezblas@kuleuven.be, seppe.dewinter@kuleuven.be",
     description="pycistarget is a python module to perform motif enrichment analysis in sets of regions with different tools and identify high confidence TF cistromes",
     long_description=open('README.rst').read(),
     url="https://github.com/aertslab/pycistarget",
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
      entry_points={"console_scripts": ["pycistarget = pycistarget.cli.pycistarget:main"]}
 )

