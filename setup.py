import os.path
from setuptools import setup, find_packages


script_path = os.path.dirname(__file__)
with open(os.path.join(script_path, 'README.md')) as f:
    readme = f.read()

with open(os.path.join(script_path, 'LICENSE')) as f:
    license = f.read()

setup(
    name='NetworkAnalysis',
    version='0.1',
    description='Package to analyze protein-protein interaction networks',
    long_description=readme,
    author='Joaquim Aguirre-Plans',
    author_email='joaquim.aguirre@upf.edu',
    url='',
    license=license,
    packages=find_packages(exclude=('NetworkAnalysis'))
)