# Ontosunburst

Sunburst visualisation of an ontology representing classes of sets
of metabolic objects

## Requirements

Requirements from `requirements.txt`

- numpy>=1.22.0
- padmet>=5.0.1
- plotly>=5.17.0
- scipy>=1.8.1
- SPARQLWrapper>=2.0.0

### Installation

```commandline
pip install -r requirements.txt
pip install -e .
```
## Utilisation

### Availabilities

#### 1. **Ontologies :**

- MetaCyc (compounds, reactions and pathways)
- ChEBI (chebi roles)
- EC (EC-numbers)

#### 2. **Analysis :**

- Proportion (1 set) : displays proportion representation of all classes
- Comparison (1 set + 1 reference set) :  displays enrichment analysis
significance of a set according to a reference set of metabolic objects

### Inputs

#### Objects sets

A python collection (list, set, ...) of metabolic objects IDs 

Metabolic objects from some tools outputs can be extracted easily with 
functions from `obj_extraction`

#### Ontology

#### 1. Files (MetaCyc, EC)

Intern tool files (in `Inputs` folder) by default

Specified user local files
- MetaCyc (Data Base Padmet Reference (metacyc_x.xx.padmet) + 
Classes ontology json file)
- EC (Classes ontology json file + Names json file)

#### 2. SPARQL server url (ChEBI)


## Run

Codes to run workflows (proportion or comparison) to create sunburst in
`ontosunburst.py`


### MetaCyc `ontosunburst.ontosunburst.metacyc_ontosunburst`

#### Parameters

---------- REQUIRED ----------

- `metabolic_objects` : `Collection[str]`
    - Set of metabolic objects to classify

---------- OPTIONAL ----------

- `reference_set` : `Collection[str]` (optional, `default=None`)
  - Set of reference metabolic objects
- `output` : `str` (optional, `default=None`)
  - Path to output to save figure
- `class_file` : `str` (optional, `default=CLASS_FILE`)
  - Path to class ontology file
- `padmet_ref` : `str` (optional, `default=METACYC_FILE`)
  - Path to metacyc padmet ref file
- `test` : `str` (optional, `default=BINOMIAL_TEST`)
  - Type of test for enrichment analysis if reference_set is 
not None
- `full` : `bool` (optional, `default=True`)
  - True to duplicate labels if +1 parents (False to take exactly 
1 random parent)
- `total` : `bool` (optional, `default=True`)
  - True to have branch values proportional of the total parent 
(may not work in some cases)

#### Example

```python
from ontosunburst.ontosunburst import metacyc_ontosunburst

MET_SET = {'CPD-24674', 'CPD-24687', 'CPD-24688'}
REF_MET = {'CPD-24674', 'CPD-24687', 'CPD-24688',
           'CPD-12782', 'CPD-12784', 'CPD-12787',
           'CPD-12788', 'CPD-12789', 'CPD-12796',
           'CPD-12797', 'CPD-12798', 'CPD-12805',
           'CPD-12806', 'CPD-12812', 'CPD-12816',
           'CPD-1282', 'CPD-12824', 'CPD-1283'}
```


```python
# PROPORTION
metacyc_ontosunburst(metabolic_objects=REF_MET, 
                     output='test')
```
![MetaCyc Proportion figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_mc_cpd_prop.png)


```python
# COMPARISON
metacyc_ontosunburst(metabolic_objects=MET_SET, 
                     reference_set=REF_MET, 
                     output='test')
```
![MetaCyc Comparison figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_mc_cpd_comp.png)


### EC `ontosunburst.ontosunburst.ec_ontosunburst`

#### Parameters

---------- REQUIRED ----------

- `ec_set` : `Collection[str]`
  - Set of EC numbers objects to classify (format "x.x.x.x" or "x.x.x.-")

---------- OPTIONAL ----------

- `reference_set` : `Collection[str]` (optional, `default=None`)
  - Set of reference metabolic objects
- `output` : `str` (optional, `default=None`)
  - Path to output to save figure
- `class_file` : `str` (optional, `default=ENZYME_ONTO_FILE`)
  - Path to class ontology file
- `names_file` : `str` (optional, `default=NAMES_FILE`)
  - Path to EC_ID - EC_NAME association json file
- `test` : `str` (optional, `default=BINOMIAL_TEST`)
  - Type of test for enrichment analysis if reference_set is 
not None
- `full` : `bool` (optional, `default=True`)
  - True to duplicate labels if +1 parents (False to take exactly 
1 random parent)
- `total` : `bool` (optional, `default=True`)
  - True to have branch values proportional of the total parent 
(may not work in some cases)

#### Example

```python
from  ontosunburst.ontosunburst import ec_ontosunburst

EC_SET = {'2.6.1.45', '1.1.1.25', '1.1.1.140'}
REF_EC = {'2.6.1.45', '1.1.1.25', '1.1.1.140',
          '1.14.14.52', '2.7.1.137', '7.1.1.8',
          '1.17.4.5', '2.3.1.165', '3.2.1.53',
          '3.2.1.91', '6.3.4.2', '5.4.99.8'}
```


```python
# PROPORTION
ec_ontosunburst(ec_set=REF_EC, 
                output='test')
```
![EC Proportion figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_ec_prop.png)


```python
# COMPARISON
ec_ontosunburst(ec_set=EC_SET, 
                reference_set=REF_EC, 
                output='test')
```
![EC Comparison figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_ec_comp.png)

### ChEBI ` ontosunburst.ontosunburst.chebi_ontosunburst`

#### Parameters

---------- REQUIRED ----------

- `chebi_ids` : `Collection[str]`
  - Set of ChEBI IDs to classify
- `endpoint_url` : `str`
  - URL of ChEBI ontology for SPARQL requests

---------- OPTIONAL ----------

- `reference_set` : `Collection[str]` (optional, `default=None`)
  - Set of reference metabolic objects
- `output` : `str` (optional, `default=None`)
  - Path to output to save figure
- `test` : `str` (optional, `default=BINOMIAL_TEST`)
  - Type of test for enrichment analysis if reference_set is 
not None
- `full` : `bool` (optional, `default=True`)
  - True to duplicate labels if +1 parents (False to take exactly 
1 random parent)
- `total` : `bool` (optional, `default=True`)
  - True to have branch values proportional of the total parent 
(may not work in some cases)

#### Example

```python
from  ontosunburst.ontosunburst import chebi_ontosunburst

URL = 'http://localhost:3030/chebi/'
CH_SET = {'38028', '28604', '85146'}
REF_CH = {'38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639'}
```


```python
# PROPORTION
chebi_ontosunburst(chebi_ids=REF_CH, 
                   endpoint_url=URL, 
                   output='test_chebi_prop')
```
![ChEBI Proportion figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_chebi_prop.png)


```python
# COMPARISON
chebi_ontosunburst(chebi_ids=CH_SET, 
                   reference_set=REF_CH, 
                   endpoint_url=URL,
                   output='test_chebi_comp')
```
![ChEBI Comparison figure](https://github.com/PaulineGHG/Ontology_sunburst/blob/main/tests/expected_figures/png/test_chebi_comp.png)
