
from os import path

from setuptools import find_packages, setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(
        name='crc',
        version='1.0.0',
        description='Core transcriptional regulatory circuitry analysis',
        long_description=long_description,
        url='https://github.com/linlabcode/CRC',
        download_url='https://github.com/linlabcode/CRC/tarball/v1.0.0',

        classifiers=[],

        keywords=['bioinformatics', 'crc', 'core', 'transcriptional', 'regulatory', 'circuitry'],

        packages=find_packages(),

        install_requires=['numpy>=1.14.5', 'networkx>=1.8.1'],
        extras_require={},

        entry_points={
            'console_scripts': [
                'crc=crc.bin.crc3:main',
                ]
            },
        include_package_data=True,
)
