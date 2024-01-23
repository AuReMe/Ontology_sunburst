from setuptools import setup

setup(
    name='ontosunburst',
    packages=['ontosunburst'],
    author="Pauline Hamon-Giraud",
    author_email="gem-aureme@inria.fr",
    description="Sunbursts for ontologies",
    long_description='Ontology_sunburst allows to represent a set of element (compounds, reactions, pathways, etc..) '
                     'throught a sunburst diagram showing proportion of ontology classes. Works for Metacyc, Chebi '
                     'roles and EC ontologies.',
    install_requires=['numpy', 'padmet', 'plotly', 'scipy', 'SPARQLWrapper'],
    url='https://github.com/AuReMe/Ontology_sunburst.git',
    license='GPLv3+',
    python_requires='>=3.9',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        # Audience
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        # Environnement, OS, languages
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)
