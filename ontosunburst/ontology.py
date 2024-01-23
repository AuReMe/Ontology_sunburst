import padmet.classes
from typing import List, Set, Dict, Tuple, Collection
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


# For MetaCyc Ontology
# --------------------------------------------------------------------------------------------------
def get_dict_ontology(padmet_ref: padmet.classes.PadmetRef):
    # TODO: fix, understand why it doesn't work
    dict_ontology = dict()
    for e, node in padmet_ref.dicOfNode.items():
        if node.type == 'class' or node.type == 'compound' and e != 'FRAMES':
            try:
                rel = padmet_ref.dicOfRelationIn[e]
                dict_ontology[e] = list()
                for r in rel:
                    if r.type == 'is_a_class':
                        dict_ontology[e].append(r.id_out)
            except KeyError:
                pass
    return dict_ontology


def extract_metacyc_classes(metabolic_objects: Collection[str],
                            d_classes_ontology: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """ Extract +1 parent classes for each metabolite.

    Parameters
    ----------
    metabolic_objects: Collection[str]
        Collection of metabolic objects.
    d_classes_ontology: Dict[str, List[str]]
        MetaCyc classes dictionary

    Returns
    -------
    Dict[str, List[str]] (Dict[metabolite, List[class]])
        Dictionary associating for each metabolic object, the list of +1 parent classes it belongs
        to.
    """
    d_obj_classes = dict()
    print(f'{len(metabolic_objects)} metabolic objects to classify')
    for obj in metabolic_objects:
        try:
            d_obj_classes[obj] = d_classes_ontology[obj]
            classified = True
        except KeyError:
            classified = False
        if not classified:
            print(f'{obj} not classified.')
    print(f'{len(d_obj_classes)}/{len(metabolic_objects)} metabolic objects classified')
    return d_obj_classes


# For EC
# --------------------------------------------------------------------------------------------------
def extract_ec_classes(ec_set: Collection[str], d_classes_ontology):
    ec_classes = dict()
    for ec in ec_set:
        parent = ec.split('.')
        parent[-1] = '-'
        parent = '.'.join(parent)
        d_classes_ontology[ec] = [parent]
        ec_classes[ec] = [parent]
    return ec_classes, d_classes_ontology


# For ChEBI Ontology
# --------------------------------------------------------------------------------------------------
def extract_chebi_roles(chebi_ids: Collection[str], endpoint_url: str) \
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]]]:
    """

    Parameters
    ----------
    chebi_ids: Collection[str]
        Collection of ChEBI IDs to extract roles associated
    endpoint_url: str
        URL endpoint string

    Returns
    -------
    Tuple[Dict[str, Set[str]], Dict[str, List[str]]]
        Dictionary of all roles associated for each ChEBI ID
        Dictionary of roles ontology, associating for each role, its parent roles
    """
    d_roles_ontology = dict()
    all_roles = dict()
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
    
        SELECT DISTINCT ?molecule ?moleculeLabel ?role ?roleLabel ?parentRoleLabel
        WHERE {{
            VALUES ?molecule{{ chebidb:{chebi_id }}}                                        
            
            ?molecule rdfs:label ?moleculeLabel .
            
            ?molecule rdfs:subClassOf+ ?restriction .
            ?restriction rdf:type owl:Restriction .
            ?restriction owl:onProperty obo:RO_0000087 .
            ?restriction owl:someValuesFrom/(rdfs:subClassOf*) ?role .
            
            ?role rdfs:subClassOf ?parentRole .
            
            ?parentRole rdfs:label ?parentRoleLabel .
            ?role rdfs:label ?roleLabel .
        }}
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        parent_roles = set()
        for result in results["results"]["bindings"]:
            role = result['roleLabel']['value']
            parent_role = result['parentRoleLabel']['value']
            parent_roles.add(parent_role)
            roles.add(result['roleLabel']['value'])
            if role not in d_roles_ontology:
                d_roles_ontology[role] = set()
            d_roles_ontology[role].add(parent_role)
        if roles:
            d_roles_ontology[chebi_id] = list(roles.difference(parent_roles))
            roles.add(ROOTS[CHEBI])
            all_roles[chebi_id] = roles
        else:
            print(f'No ChEBI role found for : {chebi_id}')

    for c, p in d_roles_ontology.items():
        d_roles_ontology[c] = list(p)

    return all_roles, d_roles_ontology


def extract_go_classes(go_ids: Collection[str], endpoint_url: str) \
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]]]:
    """

    Parameters
    ----------
    go_ids: Collection[str]
        Collection of GO IDs
    endpoint_url: str
        URL endpoint string

    Returns
    -------
    Tuple[Dict[str, Set[str]], Dict[str, List[str]]]
        Dictionary of all roles associated for each ChEBI ID
        Dictionary of roles ontology, associating for each role, its parent roles
    """
    d_classes_ontology = dict()
    all_classes = dict()
    for go in go_ids:
        go = go.lower()
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
        
        SELECT ?goLabel ?parentGoLabel
        WHERE {{
           {go} rdfs:subClassOf* ?go .
           ?go rdfs:label ?goLabel .
           ?go rdf:type owl:Class .
           ?go rdfs:subClassOf ?parentGo .
           ?parentGo rdfs:label ?parentGoLabel .
           ?parentGo rdf:type owl:Class .
        }}
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        parent_classes = set()
        for result in results["results"]["bindings"]:
            go_class = result['goLabel']['value']
            parent_class = result['parentGoLabel']['value']
            parent_classes.add(parent_class)
            go_classes.add(go_class)
            go_classes.add(parent_class)
            if parent_class in GO_ROOTS:
                d_classes_ontology[parent_class] = [ROOTS[GO]]
            if go_class not in d_classes_ontology:
                d_classes_ontology[go_class] = set()
            d_classes_ontology[go_class].add(parent_class)
        if go_classes:
            d_classes_ontology[go] = list(go_classes.difference(parent_classes))
            go_classes.add(ROOTS[GO])
            all_classes[go] = go_classes
        else:
            print(f'No GO class found for : {go}')

    for c, p in d_classes_ontology.items():
        d_classes_ontology[c] = list(p)
    return all_classes, d_classes_ontology


def extract_classes(ontology, metabolic_objects, root, d_classes_ontology=None, endpoint_url=None):
    if ontology == METACYC:
        leaf_classes = extract_metacyc_classes(metabolic_objects, d_classes_ontology)
        return get_all_classes(leaf_classes, d_classes_ontology, root), d_classes_ontology
    if ontology == EC:
        leaf_classes, d_classes_ontology = extract_ec_classes(metabolic_objects, d_classes_ontology)
        return get_all_classes(leaf_classes, d_classes_ontology, root), d_classes_ontology
    if ontology == CHEBI:
        return extract_chebi_roles(metabolic_objects, endpoint_url)
    if ontology == GO:
        return extract_go_classes(metabolic_objects, endpoint_url)
    if ontology == KEGG:
        leaf_classes = extract_metacyc_classes(metabolic_objects, d_classes_ontology)
        return get_all_classes(leaf_classes, d_classes_ontology, root), d_classes_ontology


# For all Ontology - Utils
# --------------------------------------------------------------------------------------------------

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

    if parents != [root_item]:
        for p in parents:
            parent_set = get_parents(p, parent_set, d_classes_ontology, root_item)

    return parent_set


def get_all_classes(met_classes: Dict[str, List[str]], d_classes_ontology: Dict[str, List[str]],
                    root_item: str) -> Dict[str, Set[str]]:
    """ Extract all parent classes for each metabolite.

    Parameters
    ----------
    met_classes: Dict[str, List[str]] (Dict[metabolite, List[class]])
        Dictionary associating for each metabolite the list of +1 parent classes it belongs to.
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each class its parent classes.
    root_item: str
        Name of the root item of the ontology.

    Returns
    -------
    Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.
    """
    all_classes_met = dict()
    for met, classes in met_classes.items():

        all_classes = set(classes)
        for c in classes:
            m_classes = get_parents(c, set(d_classes_ontology[c]), d_classes_ontology, root_item)
            all_classes = all_classes.union(m_classes)

        all_classes_met[met] = all_classes
    return all_classes_met


def get_classes_abondance(all_classes: Dict[str, Set[str]], show_leaves: bool) -> Dict[str, int]:
    """ Indicate for each class the number of metabolites found belonging to the class

    Parameters
    ----------
    all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.
    show_leaves: bool
        True to show input metabolic objets at sunburst leaves

    Returns
    -------
    Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    """
    classes_abondance = dict()
    for met, classes in all_classes.items():
        if show_leaves:
            classes_abondance[met] = 1
        for c in classes:
            if c not in classes_abondance.keys():
                classes_abondance[c] = 1
            else:
                classes_abondance[c] += 1
    return dict(reversed(sorted(classes_abondance.items(), key=lambda item: item[1])))


def get_children_dict(parent_dict: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """ Create the children dictionary from the parents dictionary.

    Parameters
    ----------
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes

    Returns
    -------
    Dict[str, List[str]]
        Dictionary associating for each class, its children classes
    """
    children_dict = dict()
    for c, ps in parent_dict.items():
        for p in ps:
            if p not in children_dict.keys():
                children_dict[p] = list()
            children_dict[p].append(c)
    return children_dict


def reduce_d_ontology(d_classes_ontology, classes_abundance):
    reduced_d_ontology = dict()
    for k, v in d_classes_ontology.items():
        if k in classes_abundance:
            reduced_d_ontology[k] = v
    return reduced_d_ontology
