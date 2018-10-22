
from os import path

from setuptools import find_packages, setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

# Get package metadata from 'crc/__about__.py' file
about = {}
with open(path.join(here, 'crc', '__about__.py')) as f:
    exec(f.read(), about)

setup(
    name=about['__title__'],

    version=about['__version__'],

    description=about['__summary__'],
    long_description=long_description,

    url=about['__url__'],
    download_url=about['__download_url__'],

    maintainer=about['__maintainer__'],
    maintainer_email=about['__email__'],

    license=about['__license__'],

    packages=find_packages(exclude=['tests', 'tests.*', '*.tests', '*.tests.*']),

    install_requires=['numpy>=1.14.5', 'networkx>=1.8.1'],
    extras_require={
        'package': [
            'twine',
            'wheel',
        ],
        'test': [
            'check-manifest',
            'docutils',
            'pytest-cov',
            'flake8',
            'coverage>=4.2',
        ]
    },

    entry_points={
        'console_scripts': [
            'crc=crc.bin.crc3:main',
            ]
        },

    include_package_data=True,

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords=['bioinformatics', 'crc', 'core', 'transcriptional', 'regulatory', 'circuitry'],
)

