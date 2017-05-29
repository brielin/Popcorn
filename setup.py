import os
from setuptools import setup

HERE = os.path.dirname(__file__)
version_file = open(os.path.join(HERE,'popcorn','VERSION'))
version = version_file.read().strip()

long_description = """
popcorn is a tool for estimating genetic correlation across multiple
populations from summary GWAS data. Please see Brown et. al.
'Transethnic genetic correlation estimates from summary statistics'
AJHG 2016
"""

with open(os.path.join(HERE, 'requirements.txt'), 'r') as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
    name='popcorn',
    version=version,
    install_requires=install_requires,
    requires=['python (>=2.7, <3.0)'],
    packages=['popcorn'],
    author='Brielin Brown',
    description='A tool for estimating transethnic genetic correlation',
    long_description=long_description,
    url='github.com/brielin/popcorn',
    package_dir={'popcorn': 'popcorn'},
    package_data={'popcorn': []},
    zip_safe=False,
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'popcorn=popcorn.__main__:main',
        ],
    },
    author_email="brielin.brown@gmail.com"
)
