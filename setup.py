from setuptools import setup, find_packages

def read_requirements(fname):
    with open(fname, 'r', encoding='utf-8') as file:
        return [line.rstrip() for line in file]


setup(
     name='pycistarget',
     use_scm_version=True,
     setup_requires=['setuptools_scm'],
     packages=find_packages(),
     include_dirs=["."],
     install_requires=read_requirements('requirements.txt'),
     author="Carmen Bravo, Seppe de Winter",
     author_email="carmen.bravogonzalezblas@kuleuven.be, seppe.dewinter@kuleuven.be",
     description="pycistarget is a python module to perform motif enrichment analysis in sets of regions with different tools and identify high confidence TF cistromes",
     long_description=open('README.md').read(),
     url="https://github.com/aertslab/pycistarget",
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )

