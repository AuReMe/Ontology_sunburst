# Ontosunburst

Sunburst visualisation of an ontology representing classes of sets
of metabolic objects

## Requirements

### Mandatory
Python 3.10 recommended

Requirements from `requirements.txt`

- numpy>=1.22.0
- padmet>=5.0.1
- plotly>=5.17.0
- scipy>=1.8.1
- SPARQLWrapper>=2.0.0

### Optional

Need *Apache Jena Fuseki* SPARQL server for ChEBI and GO requests and their OWL files.

- Download *Apache Jena Fuseki* : https://jena.apache.org/download/index.cgi 
- Download ChEBI ontology : https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/
  (chebi.owl or chebi_lite.owl)
- Download GO ontology : https://geneontology.org/docs/download-ontology/ (go.owl)

## Installation

```commandline
pip install -r requirements.txt
pip install -e .
```

### Set up Jena SPARQL server (optional : for ChEBI and GO)

Execute followed bash script to launch server.

#### ChEBI

```bash
#!/bin/bash

FUSEKI_PATH=/path/to/apache-jena-fuseki-x.x.x
CHEBI_PATH=/path/to//chebi_lite.owl

${FUSEKI_PATH}/fuseki-server --file=${CHEBI_PATH} /chebi
```

#### GO

```bash
#!/bin/bash

FUSEKI_PATH=/path/to/apache-jena-fuseki-x.x.x
GO_PATH=/path/to/go.owl

${FUSEKI_PATH}/fuseki-server --file=${GO_PATH} /go
```

## Utilisation

### Availabilities

#### 5 **Ontologies :**

With local files :
- MetaCyc (compounds, reactions, pathways)
- EC (EC-numbers)
- KEGG Ontology (modules, pathways, ko, ko_transporter, metabolite, metabolite_lipid)

With SPARQL server :
- ChEBI (chebi roles)
- Gene Ontology (works well only with small set)

#### 2 **Analysis :**

- Topology (**1 set** + 1 optional reference set) : displays proportion 
(number of occurrences) representation of all classes
- Enrichment (**1 set** + **1 reference set**) :  displays enrichment 
analysis significance of a set according to a reference set of metabolic 
objects

### Inputs

#### Objects sets

A python collection (list, set, ...) of metabolic objects IDs 

Metabolic objects from some tools outputs can be extracted easily with 
functions from `obj_extraction`

#### Ontology

#### 1. Files (MetaCyc, EC, KEGG)

Intern tool files (in `Inputs` folder) by default

Specified user local files
- MetaCyc (Classes ontology json file)
- EC (Classes ontology json file + Names json file)
- KEGG (Classes ontology json file)

#### 2. SPARQL server url (ChEBI, GO)

See **Set up Jena SPARQL server (optional : for ChEBI and GO)** for 
install guide.



## Main workflow function : `ontosunburst.ontosunburst`

### Parameters

---------- REQUIRED ----------

`ontology` : `str`
- Ontology to use, must be in : [metacyc, ec, kegg, chebi, go]

`metabolic_objects` : `Collection[str]`
- Set of metabolic objects to classify

---------- OPTIONAL ----------

`reference_set` : `Collection[str]`  `(optional, default=None)`
- Set of reference metabolic objects

`analysis` : `str` `(optional, default=topology)`
- Analysis mode, must be in : `[topology, enrichment]

`output` : `str` `(optional, default=None)`
- Path to output to save figure

`class_file` : `str` `(optional, default=None)`
- Path to class ontology file

`names_file` : `str` `(optional, default=None)`
- Path to EC_ID - EC_NAME association json file

`endpoint_url` : `str` `(optional, default=None)`
- URL of ChEBI ontology for SPARQL requests

`test` : `str` `(optional, default=binomial)`
- Type of test if analysis=enrichment, must be in : [binomial, hypergeometric]
    
`full` : `bool` `(optional, default=True)`
- True to duplicate labels if +1 parents (False to take exactly 1 random parent)
    
`total` : `bool` `(optional, default=True)`
- True to have branch values proportional of the total parent (may not work in some cases)
    
`root_cut` : `str` `(optional, default=ROOT_CUT)`
- mode for root cutting (uncut, cut, total)
    
`ref_base` : `bool` `(optional, default=False)`
- True to have the base classes representation of the reference set in the figure.
    
`show_leaves` : `bool` `(optional, default=False)`
- True to show input metabolic objets at sunburst leaves


### MetaCyc

#### Example

```python
from ontosunburst.ontosunburst import ontosunburst

MET_SET = {'CPD-24674', 'CPD-24687', 'CPD-24688'}
REF_MET = {'CPD-24674', 'CPD-24687', 'CPD-24688',
           'CPD-12782', 'CPD-12784', 'CPD-12787',
           'CPD-12788', 'CPD-12789', 'CPD-12796',
           'CPD-12797', 'CPD-12798', 'CPD-12805',
           'CPD-12806', 'CPD-12812', 'CPD-12816',
           'CPD-1282', 'CPD-12824', 'CPD-1283'}

# TOPOLOGY
fig = ontosunburst(ontology='metacyc', 
                   metabolic_objects=MET_SET, 
                   reference_set=REF_MET,
                   output='test_mc_cpd_prop', 
                   ref_base=True)

# ENRICHMENT
fig = ontosunburst(ontology='metacyc', 
                   metabolic_objects=MET_SET, 
                   reference_set=REF_MET,
                   analysis='enrichment', 
                   output='test_mc_cpd_comp', 
                   ref_base=True)
```

![MetaCyc Proportion figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_mc_cpd_prop.png)

![MetaCyc Comparison figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_mc_cpd_comp.png)


### EC

#### Example

```python
from  ontosunburst.ontosunburst import ontosunburst

EC_SET = {'2.6.1.45', '1.1.1.25', '1.1.1.140'}
REF_EC = {'2.6.1.45', '1.1.1.25', '1.1.1.140',
          '1.14.14.52', '2.7.1.137', '7.1.1.8',
          '1.17.4.5', '2.3.1.165', '3.2.1.53',
          '3.2.1.91', '6.3.4.2', '5.4.99.8'}

# TOPOLOGY
fig = ontosunburst(ontology='ec', 
                   metabolic_objects=EC_SET, 
                   reference_set=REF_EC,
                   output='test_ec_prop', 
                   ref_base=True, 
                   show_leaves=True)

# ENRICHMENT
fig = ontosunburst(ontology='ec', 
                   metabolic_objects=EC_SET, 
                   reference_set=REF_EC,
                   output='test_ec_comp', 
                   analysis='enrichment', 
                   ref_base=True)

```

![EC Proportion figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_ec_prop.png)

![EC Comparison figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_ec_comp.png)

### ChEBI

#### Example

```python
from  ontosunburst.ontosunburst import ontosunburst

URL = 'http://localhost:3030/chebi/'
CH_SET = {'38028', '28604', '85146'}
REF_CH = {'38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639'}

# TOPOLOGY
fig = ontosunburst(ontology='chebi', 
                   metabolic_objects=CH_SET, 
                   reference_set=REF_CH,
                   endpoint_url=URL, 
                   output='test_chebi_prop', 
                   ref_base=True)

# ENRICHMENT
fig = ontosunburst(ontology='chebi', 
                   metabolic_objects=CH_SET, 
                   reference_set=REF_CH,
                   endpoint_url=URL, 
                   output='test_chebi_comp', 
                   analysis='enrichment', 
                   ref_base=True)
```

![ChEBI Proportion figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_chebi_prop.png)

![ChEBI Comparison figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_chebi_comp.png)
