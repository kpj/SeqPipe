from setuptools import setup, find_packages


setup(
    name='seqpipe',
    version='0.0.4',

    description='Sequencing pipeline',
    long_description_markdown_filename='README.md',

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

    setup_requires=['setuptools-markdown'],
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
