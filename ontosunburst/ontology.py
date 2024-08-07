from typing import List, Set, Dict, Tuple

import numpy
from SPARQLWrapper import SPARQLWrapper, JSON


# CONSTANTS ========================================================================================


METACYC = 'metacyc'
EC = 'ec'
CHEBI = 'chebi'
GO = 'go'
KEGG = 'kegg'

ROOTS = {METACYC: 'FRAMES',
         CHEBI: 'role',
         EC: 'Enzyme',
         GO: 'GO',
         KEGG: 'kegg'}

GO_ROOTS = ['cellular_component', 'biological_process', 'molecular_function']


# ==================================================================================================
# CLASSES EXTRACTION
# ==================================================================================================

# Main class extraction function
# --------------------------------------------------------------------------------------------------
def extract_classes(ontology: str, concepts: List[str], root: str,
                    d_classes_ontology: Dict[str, List[str]] = None, endpoint_url: str = None,
                    names: Dict[str, str] = None)\
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]], Dict[str, str] or None]:
    """ Extract all parent classes (until root) from a list of concepts.

    Parameters
    ----------
    ontology: str
        Name of the ontology considered
    concepts: List[str]
        List of concepts to classify
    root: str
        Root item of the ontology
    d_classes_ontology: Dict[str, List[str]], optional (default=None)
        Dictionary of the classes ontology associating for each concept its +1 parent classes.
    endpoint_url: str, optional (default=None)
        URL for the SPARQL server (for GO and ChEBI ontologies)
    names: Dict[str, str] (default=None)

    Returns
    -------
    Dict[str, Set[str]]
        Dictionary associating to each concept all its parent classes (until the root)
    Dict[str, List[str]]
        Dictionary of the classes ontology associating for each concept its +1 parent classes.
    Dict[str, str]  or None
        Dictionary of labels
    """
    if ontology == METACYC or ontology == KEGG or ontology is None:
        leaf_classes = extract_met_classes(concepts, d_classes_ontology)
        return get_all_classes(leaf_classes, d_classes_ontology, root), d_classes_ontology, names
    if ontology == EC:
        leaf_classes, d_classes_ontology = extract_ec_classes(concepts, d_classes_ontology)
        return get_all_classes(leaf_classes, d_classes_ontology, root), d_classes_ontology, names
    if ontology == CHEBI:
        return extract_chebi_roles(concepts, endpoint_url)
    if ontology == GO:
        return extract_go_classes(concepts, endpoint_url)


# For MetaCyc and Kegg Ontology
# --------------------------------------------------------------------------------------------------
def extract_met_classes(concepts: List[str], d_classes_ontology: Dict[str, List[str]]) \
        -> Dict[str, List[str]]:
    """ Extract +1 parent classes for each concept considered.

    Parameters
    ----------
    concepts: List[str]
        List of concepts considered
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each concept its +1 parent classes.

    Returns
    -------
    Dict[str, List[str]]
        Dictionary associating for each concept, the list of +1 parent classes it belongs to.
    """
    d_obj_classes = dict()
    print(f'{len(concepts)} metabolic objects to classify')
    for obj in concepts:
        try:
            d_obj_classes[obj] = d_classes_ontology[obj]
            classified = True
        except KeyError:
            classified = False
        if not classified:
            print(f'{obj} not classified.')
    print(f'{len(d_obj_classes)}/{len(concepts)} metabolic objects classified')
    return d_obj_classes


# For EC
# --------------------------------------------------------------------------------------------------
def extract_ec_classes(ec_list: List[str], d_classes_ontology: Dict[str, List[str]]) \
        -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """ Extract +1 parent classes for each EC number.

    Parameters
    ----------
    ec_list: List[str]
        List of EC numbers
    d_classes_ontology: Dict[str, List[str]]
        EC Ontology classes dictionary : Dict[object, List[parents]]

    Returns
    -------
    Dict[str, List[str]]
        Dictionary associating for each metabolic object, the list of +1 parent classes it belongs
        to.
    Dict[str, List[str]]
        EC Ontology classes dictionary : Dict[object, List[parents]]
    """
    print(f'{len(ec_list)} EC numbers to classify')
    ec_classes = dict()
    for ec in ec_list:
        dec_ec = ec.split('.')
        while len(dec_ec) < 4:
            dec_ec.append('-')
        if dec_ec.count('-') == 3 and ec in d_classes_ontology:
            parent = ROOTS[EC]
        else:
            for i in range(3):
                if dec_ec.count('-') == i:
                    dec_ec[3-i] = '-'
                    break
            parent = '.'.join(dec_ec)
        if parent in d_classes_ontology or parent == ROOTS[EC]:
            d_classes_ontology[ec] = [parent]
            ec_classes[ec] = [parent]
        else:
            print(f'{ec} not classified')
    print(f'{len(ec_classes)}/{len(ec_list)} EC numbers classified')
    return ec_classes, d_classes_ontology


# For ChEBI Ontology
# --------------------------------------------------------------------------------------------------
def extract_chebi_roles(chebi_ids: List[str], endpoint_url: str) \
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]], Dict[str, str]]:
    """ Extract all parent classes for each chebi ID + Generate ontology dictionary.

    Parameters
    ----------
    chebi_ids: List[str]
        List of ChEBI IDs to extract roles associated
    endpoint_url: str
        URL endpoint string

    Returns
    -------
    Tuple[Dict[str, Set[str]], Dict[str, List[str]], Dict[str, str]]
        Dictionary of all roles associated for each ChEBI ID
        Dictionary of roles ontology, associating for each role, its parent roles
        Dictionary of labels
    """
    root_id = '50906'
    d_roles_ontology = dict()
    all_roles = dict()
    d_labels = dict()
    chebi_ok = 0
    total_nb = len(chebi_ids)
    for chebi_id in chebi_ids:
        roles = set()
        sparql = SPARQLWrapper(endpoint_url)
        sparql.setQuery(f"""
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        PREFIX dc: <http://purl.org/dc/elements/1.1/>
        PREFIX dcterms: <http://purl.org/dc/terms/>
        PREFIX chebi: <http://purl.obolibrary.org/obo/chebi/>
        PREFIX chebidb: <http://purl.obolibrary.org/obo/CHEBI_>
        PREFIX chebirel: <http://purl.obolibrary.org/obo/chebi#>
        PREFIX obo: <http://purl.obolibrary.org/obo/>
        PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
        PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
    
        SELECT DISTINCT ?molecule ?moleculeLabel ?roleId ?roleLabel ?parentRoleLabel ?parentRoleId
        WHERE {{
            VALUES ?molecule{{ chebidb:{chebi_id}}}                                        
            
            ?molecule rdfs:label ?moleculeLabel .
            
            ?molecule rdfs:subClassOf+ ?restriction .
            ?restriction rdf:type owl:Restriction .
            ?restriction owl:onProperty obo:RO_0000087 .
            ?restriction owl:someValuesFrom/(rdfs:subClassOf*) ?role .
            
            ?role oboInOwl:id ?roleId .
            ?role rdfs:subClassOf ?parentRole .
            
            ?parentRole oboInOwl:id ?parentRoleId .
            ?parentRole rdfs:label ?parentRoleLabel .
            ?role rdfs:label ?roleLabel .
        }}
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        parent_roles = set()
        for result in results["results"]["bindings"]:
            molecule_label = result['moleculeLabel']['value']
            role_label = result['roleLabel']['value']
            role_id = result['roleId']['value'].split(':')[1]
            parent_role_label = result['parentRoleLabel']['value']
            parent_role_id = result['parentRoleId']['value'].split(':')[1]
            if role_id == root_id:
                role_id = ROOTS[CHEBI]
            if parent_role_id == root_id:
                parent_role_id = ROOTS[CHEBI]
            d_labels[role_id] = role_label
            d_labels[parent_role_id] = parent_role_label
            d_labels[chebi_id] = molecule_label
            parent_roles.add(parent_role_id)
            roles.add(role_id)
            if role_id not in d_roles_ontology:
                d_roles_ontology[role_id] = set()
            d_roles_ontology[role_id].add(parent_role_id)
        if roles:
            chebi_ok += 1
            d_roles_ontology[chebi_id] = list(roles.difference(parent_roles))
            roles.add(ROOTS[CHEBI])
            all_roles[chebi_id] = roles
        else:
            print(f'No ChEBI role found for : {chebi_id}')

    for c, p in d_roles_ontology.items():
        d_roles_ontology[c] = list(p)
    print(f'{chebi_ok}/{total_nb} chebi id with roles associated.')
    return all_roles, d_roles_ontology, d_labels


# For GO Ontology
# --------------------------------------------------------------------------------------------------
def extract_go_classes(go_ids: List[str], endpoint_url: str) \
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]], Dict[str, str]]:
    """ Extract all parent classes for each GO ID + Generate ontology dictionary.

    Parameters
    ----------
    go_ids: List[str]
        List of GO IDs
    endpoint_url: str
        URL endpoint string

    Returns
    -------
    Tuple[Dict[str, Set[str]], Dict[str, List[str]], Dict[str, str]]
        Dictionary of all roles associated for each ChEBI ID
        Dictionary of roles ontology, associating for each role, its parent roles
        Dictionary of labels
    """
    d_classes_ontology = dict()
    all_classes = dict()
    d_labels = dict()
    for go in go_ids:
        go_classes = set()
        sparql = SPARQLWrapper(endpoint_url)
        sparql.setQuery(f"""
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        PREFIX dc: <http://purl.org/dc/elements/1.1/>
        PREFIX dcterms: <http://purl.org/dc/terms/>
        PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX uniprot: <http://purl.uniprot.org/uniprot/>
        PREFIX up:<http://purl.uniprot.org/core/>
        PREFIX go: <http://purl.obolibrary.org/obo/GO_>
        PREFIX goavoc: <http://bio2rdf.org/goa_vocabulary:>
        
        SELECT ?goLabel ?parentGoLabel ?goId ?parentGoId
        WHERE {{
           {go} rdfs:subClassOf* ?go .
           ?go oboInOwl:id ?goId .
           ?go rdfs:label ?goLabel .
           ?go rdf:type owl:Class .
           ?go rdfs:subClassOf ?parentGo .
           ?parentGo rdfs:label ?parentGoLabel .
           ?parentGo oboInOwl:id ?parentGoId .
           ?parentGo rdf:type owl:Class .
        }}
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            go_id = result['goId']['value'].lower()
            go_label = result['goLabel']['value']
            parent_id = result['parentGoId']['value'].lower()
            parent_label = result['parentGoLabel']['value']
            d_labels[go_id] = go_label
            d_labels[parent_id] = parent_label
            go_classes.add(parent_id)
            if parent_label in GO_ROOTS:
                d_classes_ontology[parent_id] = [ROOTS[GO]]
            if go_id not in d_classes_ontology:
                d_classes_ontology[go_id] = []
            d_classes_ontology[go_id].append(parent_id)
        if go_classes:
            go_classes.add(ROOTS[GO])
            all_classes[go] = go_classes
        else:
            print(f'No GO class found for : {go}')

    for c, p in d_classes_ontology.items():
        d_classes_ontology[c] = list(set(p))
    return all_classes, d_classes_ontology, d_labels


# Recursive class extraction function
# --------------------------------------------------------------------------------------------------
def get_all_classes(obj_classes: Dict[str, List[str]], d_classes_ontology: Dict[str, List[str]],
                    root_item: str) -> Dict[str, Set[str]]:
    """ Extract all parent classes for each metabolite.

    Parameters
    ----------
    obj_classes: Dict[str, List[str]] (Dict[metabolite, List[class]])
        Dictionary associating for each object the list of +1 parent classes it belongs to.
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each class its +1 parent classes.
    root_item: str
        Name of the root item of the ontology.

    Returns
    -------
    Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.
    """
    all_classes_met = dict()
    for met, classes in obj_classes.items():
        all_classes = set(classes)
        for c in classes:
            if c != root_item:
                m_classes = get_parents(c, set(d_classes_ontology[c]), d_classes_ontology, root_item)
                all_classes = all_classes.union(m_classes)
        all_classes_met[met] = all_classes
    return all_classes_met


def get_parents(child: str, parent_set: Set[str], d_classes_ontology: Dict[str, List[str]],
                root_item) -> Set[str]:
    """ Get recursively from a child class, all its parents classes found in ontology.

    Parameters
    ----------
    child: str
        Child class
    parent_set: Set[str]
        Set of all parents from previous classes
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology of MetaCyc associating for each class its parent classes.
    root_item: str
        Name of the root item of the ontology

    Returns
    -------
    Set[str]
        Set of the union of the set  of child parent classes and the set of all previous parents.
    """
    parents = d_classes_ontology[child]
    for p in parents:
        parent_set.add(p)
    for p in parents:
        if p != root_item:
            parent_set = get_parents(p, parent_set, d_classes_ontology, root_item)
    return parent_set


# ==================================================================================================
# ABUNDANCES CALCULATION
# ==================================================================================================

def get_abundance_dict(abundances: List[float] or None, metabolic_objects: List[str], ref: bool)\
        -> Dict[str, float]:
    """ Generate abundances dictionary.

    Parameters
    ----------
    abundances: List[float] (size N) or None
        List of metabolic objects abundances (or None if no abundances associated --> will associate
        an abundance of 1 for each object)
    metabolic_objects: List[str] (size N)
        List of metabolic objects ID.
    ref: bool
        True if metabolic objects are the reference list, False otherwise (subset / study case).

    Returns
    -------
    Dict[str, float]
        Dictionary associating to each object its abundance.
    """
    if abundances is None:
        abundances = len(metabolic_objects) * [1]
    if len(metabolic_objects) == len(abundances):
        abundances_dict = {}
        for i in range(len(metabolic_objects)):
            abundances_dict[metabolic_objects[i]] = abundances[i]
    else:
        if ref:
            raise AttributeError(f'Length of "reference_set" parameter must be equal to '
                                 f'"ref_abundances" parameter length : {len(metabolic_objects)} '
                                 f'!= {len(abundances)}')
        else:
            raise AttributeError(f'Length of "metabolic_objects" parameter must be equal to '
                                 f'"abundances" parameter length : {len(metabolic_objects)} '
                                 f'!= {len(abundances)}')
    return abundances_dict


def get_classes_abundance(all_classes: Dict[str, Set[str]], abundances_dict: Dict[str, float],
                          show_leaves: bool) -> Dict[str, float]:
    """ Indicate for each class the number of base object found belonging to the class

    Parameters
    ----------
    all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each concept the list of all parent classes it belongs to.
    abundances_dict: Dict[str, float]
        Dictionary associating for each concept, its abundance value
    show_leaves: bool
        True to show input metabolic objets at sunburst leaves

    Returns
    -------
    Dict[str, float]
        Dictionary associating for each class the weight of concepts found belonging to the class.
    """
    classes_abondance = dict()
    for met, classes in all_classes.items():
        if show_leaves:
            if met not in classes_abondance.keys():
                classes_abondance[met] = abundances_dict[met]
            else:
                classes_abondance[met] += abundances_dict[met]
        for c in classes:
            if c not in classes_abondance.keys():
                classes_abondance[c] = abundances_dict[met]
            else:
                classes_abondance[c] += abundances_dict[met]
    return dict(reversed(sorted(classes_abondance.items(), key=lambda item: item[1])))


def get_classes_scores(all_classes, scores_dict, root):
    classes_scores = dict()
    for met, classes in all_classes.items():
        if met in scores_dict.keys():
            classes_scores[met] = scores_dict[met]
        else:
            classes_scores[met] = numpy.nan
    classes_scores[root] = numpy.nan
    return classes_scores


# ==================================================================================================
# UTILS
# ==================================================================================================

def reduce_d_ontology(d_classes_ontology: Dict[str, List[str]],
                      classes_abundance: Dict[str, float]) -> Dict[str, List[str]]:
    """ Extract the sub-graph of the d_classes_ontology dictionary conserving only nodes implicated
    with the concepts studied.

    Parameters
    ----------
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the ontology complete graph
    classes_abundance: Dict[str, float]
        Dictionary of abundances (keys are all nodes implicated to be conserved)

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of the ontology sub-graph conserving only nodes implicated with the concepts
        studied.
    """
    reduced_d_ontology = dict()
    for k, v in d_classes_ontology.items():
        if k in classes_abundance:
            reduced_d_ontology[k] = v
    return reduced_d_ontology
