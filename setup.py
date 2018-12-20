"""
ppfit: Interatomic Potential Fitting Tools
"""

from setuptools import setup, find_packages
from ppfit import __version__ as VERSION

readme = 'README.md'
long_description = open( readme ).read()

scripts = [ 'wannier2dipoles' ]

setup(
    name='ppfit',
    version=VERSION,
    description='Interatomic potential fitting tools',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Benjamin J. Morgan',
    author_email='bjm42@bath.ac.uk',
    url='https://github.com/bjmorgan/ppfit', 
    download_url='https://github.com/bjmorgan/ppfit/archive/{}.tar.gz'.format( VERSION ),
    keywords=['pimaim'], # keywords
    packages=find_packages( exclude=['docs', 'tests*'] ),
    #package_data={ 'vasppy': ['data/*.yaml'] },
    entry_points={ 'console_scripts': [
                       '{} = ppfit.scripts.{}:main'.format( s, s ) for s in scripts ] },
    license='MIT',
    install_requires=open( 'requirements.txt' ).read(),
    python_requires='>=3.5'
    )
