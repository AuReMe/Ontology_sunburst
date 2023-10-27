import padmet.classes
from typing import List, Set, Dict, Tuple, Collection
from SPARQLWrapper import SPARQLWrapper, JSON


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
                            padmet_ref: padmet.classes.PadmetRef) -> Dict[str, List[str]]:
    """ Extract +1 parent classes for each metabolite.

    Parameters
    ----------
    metabolic_objects: Collection[str]
        Collection of metabolic objects.
    padmet_ref: padmet.classes.PadmetRef
        The PadmetRef instance

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
            rel = padmet_ref.dicOfRelationIn[obj]
        except KeyError:
            rel = None
        classified = False

        if rel is not None:
            for r in rel:
                if r.type == 'is_a_class':
                    classified = True
                    if obj not in d_obj_classes.keys():
                        d_obj_classes[obj] = list()
                    d_obj_classes[obj].append(r.id_out)
        if not classified:
            print(f'{obj} not classified.')
    print(f'{len(d_obj_classes)}/{len(metabolic_objects)} metabolic objects classified')
    return d_obj_classes


# For EC
# --------------------------------------------------------------------------------------------------
def extract_ec_classes(ec_set: Collection[str]):
    ec_classes = dict()
    for ec in ec_set:
        parent = ec.split('.')
        parent[-1] = '-'
        parent = '.'.join(parent)
        ec_classes[ec] = [parent]
    return ec_classes


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
        molecule = chebi_id
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

        for result in results["results"]["bindings"]:
            role = result['roleLabel']['value']
            parent_role = result['parentRoleLabel']['value']
            roles.add(result['roleLabel']['value'])
            if role not in d_roles_ontology:
                d_roles_ontology[role] = set()
            d_roles_ontology[role].add(parent_role)
        roles.add('role')
        all_roles[chebi_id] = roles

    for c, p in d_roles_ontology.items():
        d_roles_ontology[c] = list(p)

    return all_roles, d_roles_ontology


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


def get_classes_abondance(all_classes: Dict[str, Set[str]]) -> Dict[str, int]:
    """ Indicate for each class the number of metabolites found belonging to the class

    Parameters
    ----------
    all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.

    Returns
    -------
    Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    """
    classes_abondance = dict()
    for met, classes in all_classes.items():
        # classes_abondance[met] = 1
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


def select_root(count_input: Dict[str, int], parent: Dict[str, List[str]], act_root) \
        -> Tuple[Set[str], str, Dict[str, List[str]]]:
    # TODO: Fix, understand why it doesn't works everytime
    """ Select the root to use and return tax IDs to exclude from the figure. Keep as root the tax
    ID shared by all taxa having the lowest taxonomic rank.

    Parameters
    ----------
    count_input: Dict[str, int]
        Dictionary associating to each class, the count.
    parent: Dict[str, List[str]]
        Dictionary associating to each class, its parent classes.

    Returns
    -------
    Tuple[Set[str], str, Dict[str, List[str]]]
        Set of tax IDs to exclude and Parent dictionary modified for the root ID.
    """
    max_count = max(count_input.values())
    roots = {x for x, y in count_input.items() if y == max_count}
    root = 1
    for met, ch_id in parent.items():
        if ch_id != [act_root]:
            for ch in ch_id:
                if parent[ch][0] in roots and ch not in roots:
                    root = parent[ch][0]
    parent[root] = ['']
    return set(roots).difference({root}), root, parent
