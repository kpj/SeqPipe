import os
from setuptools import setup, find_packages


own_path = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(own_path, 'README.md'), encoding='utf-8') as fd:
    long_description = fd.read()

setup(
    name='seqpipe',
    version='0.0.3',

    description='Sequencing pipeline',
    long_description=long_description,

    url='https://github.com/kpj/SeqPipe',

    author='kpj',
    author_email='kpjkpjkpjkpjkpjkpj@gmail.com',

    license='MIT',

    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.6',
    ],

    keywords='bioinformatics sequencing pipeline alignment mapping',
    packages=find_packages(exclude=['tests']),

    install_requires=[
        'numpy', 'pandas', 'seaborn', 'matplotlib', 'colorama',
        'tqdm', 'biopython', 'pysam', 'joblib', 'click', 'sh'],
    extras_require={
        'test': ['tox', 'pytest', 'coverage'],
    },

    entry_points={
        'console_scripts': [
            'seqpipe=seqpipe:run',
        ],
    },
)
