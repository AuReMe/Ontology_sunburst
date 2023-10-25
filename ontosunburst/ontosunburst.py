import os
import json
from typing import Collection
import plotly.graph_objects as go

from padmet.classes import PadmetRef

from ontosunburst.ontology import get_all_classes, get_classes_abondance, get_children_dict, \
    extract_chebi_roles, extract_metacyc_classes, extract_ec_classes

from ontosunburst.sunburst_fig import get_fig_parameters, get_data_prop_diff, get_data_proportion, \
    generate_sunburst_fig, BINOMIAL_TEST


# CONSTANTS ========================================================================================

# DEFAULT FILES
# -------------
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# For MetaCyc
CLASS_FILE = os.path.join(CURRENT_DIR, 'Inputs/classes.json')
METACYC_FILE = os.path.join(CURRENT_DIR, 'Inputs/metacyc_26.0_prot70.padmet')
# For EC numbers
ENZYME_ONTO_FILE = os.path.join(CURRENT_DIR, 'Inputs/enzymes_ontology.json')
NAMES_FILE = os.path.join(CURRENT_DIR, 'Inputs/enzymes_class_names.json')

# Other
# -----
COMPARISON_METHOD = 'comparison'
PROPORTION_METHOD = 'proportion'

# ONTOLOGIES
# ----------
METACYC_ROOT = 'FRAMES'
CHEBI_ROLE_ROOT = 'role'
EC_ROOT = 'Enzyme'


# WORKFLOW =========================================================================================
def metacyc_ontosunburst(metabolic_objects: Collection[str], reference_set: Collection[str] = None,
                         output: str = None, class_file: str = CLASS_FILE,
                         padmet_ref: str = METACYC_FILE, test: str = BINOMIAL_TEST,
                         full: bool = True) -> go.Figure:
    """ Classify and plot a sunburst from a list of metabolic objects with MetaCyc ontology Ids

    Parameters
    ----------
    metabolic_objects: Collection[str]
        Set of metabolic objects to classify
    reference_set: Collection[str] (optional, default=None)
        Set of reference metabolic objects
    output: str (optional, default=None)
        Path to output to save figure
    class_file: str (optional, default=CLASS_FILE)
        Path to class ontology file
    padmet_ref: str (optional, default=METACYC_FILE)
        Path to metacyc padmet ref file
    test: str (optional, default=BINOMIAL_TEST)
        Type of test for enrichment analysis if reference_set is not None
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    # Load files
    padmet_ref = PadmetRef(padmet_ref)
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)

    # Extract set information
    obj_leaf_classes = extract_metacyc_classes(metabolic_objects, padmet_ref)
    obj_all_classes = get_all_classes(obj_leaf_classes, d_classes_ontology, METACYC_ROOT)
    classes_abundance = get_classes_abondance(obj_all_classes)

    # Comparison figure
    if reference_set is not None:
        ref_leaf_classes = extract_metacyc_classes(reference_set, padmet_ref)
        ref_all_classes = get_all_classes(ref_leaf_classes, d_classes_ontology, METACYC_ROOT)
        ref_classes_abundance = get_classes_abondance(ref_all_classes)
        data = get_fig_parameters(classes_abundance, d_classes_ontology,
                                  get_children_dict(d_classes_ontology), METACYC_ROOT, full)
        data = get_data_prop_diff(data, ref_classes_abundance)
        if output is not None:
            write_met_classes(ref_all_classes, output, padmet_ref)
        return generate_sunburst_fig(data, output, COMPARISON_METHOD, ref_classes_abundance, test)

    # Proportion figure
    else:
        data = get_fig_parameters(classes_abundance, d_classes_ontology,
                                  get_children_dict(d_classes_ontology), METACYC_ROOT, full)
        data = get_data_proportion(data)
        if output is not None:
            write_met_classes(obj_all_classes, output, padmet_ref)
        return generate_sunburst_fig(data, output, PROPORTION_METHOD)


def chebi_ontosunburst(chebi_ids: Collection[str], endpoint_url: str,
                       reference_set: Collection[str] = None, output: str = None,
                       test: str = BINOMIAL_TEST, full: bool = True) -> go.Figure:
    """ Classify and plot a sunburst from a list of ChEBI IDs with ChEBI roles ontology

    Parameters
    ----------
    chebi_ids: Collection[str]
        Set of ChEBI IDs to classify
    endpoint_url: str
        URL of ChEBI ontology for SPARQL requests
    reference_set: Collection[str] (optional, default=None)
        Set of reference ChEBI IDs
    output: str (optional, default=None)
        Path to output to save figure
    test: str (optional, default=BINOMIAL_TEST)
        Type of test for enrichment analysis if reference_set is not None
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    # Extract set information
    all_classes, d_roles_ontology = extract_chebi_roles(chebi_ids, endpoint_url)
    classes_abondance = get_classes_abondance(all_classes)

    # Comparison figure
    if reference_set is not None:
        ref_all_roles, d_roles_ontology = extract_chebi_roles(reference_set, endpoint_url)
        ref_roles_abundance = get_classes_abondance(ref_all_roles)
        data = get_fig_parameters(classes_abondance, d_roles_ontology,
                                  get_children_dict(d_roles_ontology), CHEBI_ROLE_ROOT, full)
        data = get_data_prop_diff(data, ref_roles_abundance)
        return generate_sunburst_fig(data, output, COMPARISON_METHOD, ref_roles_abundance, test)

    # Proportion figure
    else:
        data = get_fig_parameters(classes_abondance, d_roles_ontology,
                                  get_children_dict(d_roles_ontology), CHEBI_ROLE_ROOT, full)
        data = get_data_proportion(data)
        return generate_sunburst_fig(data, output, PROPORTION_METHOD)


def ec_ontosunburst(ec_set: Collection[str], reference_set: Collection[str] = None,
                    output: str = None, class_file: str = ENZYME_ONTO_FILE,
                    names_file: str = NAMES_FILE, test: str = BINOMIAL_TEST,
                    full: bool = True) -> go.Figure:
    """ Classify and plot a sunburst from a list of EC numbers with EC ontology Ids

    Parameters
    ----------
    ec_set: Collection[str]
        Set of EC numbers objects to classify (format "x.x.x.x" or "x.x.x.-")
    reference_set: Collection[str] (optional, default=None)
        Set of reference EC numbers
    output: str (optional, default=None)
        Path to output to save figure
    class_file: str (optional, default=ENZYME_ONTO_FILE)
        Path to class ontology file
    names_file: str (optional, default=NAMES_FILE)
        Path to EC_ID - EC_NAME association json file
    test: str (optional, default=BINOMIAL_TEST)
        Type of test for enrichment analysis if reference_set is not None
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    # Load files
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    with open(names_file, 'r') as f:
        names = json.load(f)

    # Extract set information
    ec_classes = extract_ec_classes(ec_set)
    all_classes = get_all_classes(ec_classes, d_classes_ontology, EC_ROOT)
    classes_abundance = get_classes_abondance(all_classes)

    # Comparison figure
    if reference_set is not None:
        ref_leaf_classes = extract_ec_classes(reference_set)
        ref_all_classes = get_all_classes(ref_leaf_classes, d_classes_ontology, EC_ROOT)
        ref_classes_abundance = get_classes_abondance(ref_all_classes)
        data = get_fig_parameters(classes_abundance, d_classes_ontology,
                                  get_children_dict(d_classes_ontology), EC_ROOT, full, names)
        print(data)
        data = get_data_prop_diff(data, ref_classes_abundance)
        return generate_sunburst_fig(data, output, COMPARISON_METHOD, ref_classes_abundance, test,
                                     True)

    # Proportion figure
    else:
        data = get_fig_parameters(classes_abundance, d_classes_ontology,
                                  get_children_dict(d_classes_ontology), EC_ROOT, full, names)
        data = get_data_proportion(data)
        return generate_sunburst_fig(data, output, PROPORTION_METHOD)


# EXTRAS ===========================================================================================

def write_met_classes(all_classes, output, pref):
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['Compound', 'Classes', 'Common names', 'MetaCyc link']) + '\n')
        for met, classes, in all_classes.items():
            try:
                name = ' / '.join(pref.dicOfNode[met].misc['COMMON-NAME'])
            except KeyError:
                name = ''
            link = f'https://metacyc.org/compound?orgid=META&id={met}'
            f.write('\t'.join([met, ', '.join(classes), name, link]) + '\n')
