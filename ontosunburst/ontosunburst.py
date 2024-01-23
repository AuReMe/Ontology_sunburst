import os
import json
from typing import Collection
import plotly.graph_objects as go

from ontosunburst.ontology import get_classes_abondance, get_children_dict, extract_classes, \
    reduce_d_ontology, METACYC, CHEBI, EC, GO, KEGG, ROOTS

from ontosunburst.sunburst_fig import get_fig_parameters, get_data_proportion, \
    generate_sunburst_fig, BINOMIAL_TEST, TOPOLOGY_A, ENRICHMENT_A, ROOT_CUT

# ==================================================================================================
#                                           CONSTANTS
# ==================================================================================================

# DEFAULT FILES
# -------------
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# For MetaCyc
METACYC_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'MetaCyc26_0_classes.json')
# For EC numbers
EC_ONTO_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_ontology.json')
EC_NAMES_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_class_names.json')
# For KEGG
KEGG_ONTO_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'kegg_onto.json')
# For ChEBI
CHEBI_URL = 'http://localhost:3030/chebi/'
# For GO
GO_URL = 'http://localhost:3030/go/'


# ==================================================================================================
#                                            WORKFLOW
# ==================================================================================================

def ontosunburst(ontology: str,
                 metabolic_objects: Collection[str],
                 reference_set: Collection[str] = None,
                 analysis: str = TOPOLOGY_A,
                 output: str = None,
                 class_file: str = None,
                 names_file: str = None,
                 endpoint_url: str = None,
                 test: str = BINOMIAL_TEST,
                 full: bool = True,
                 total: bool = True,
                 root_cut: str = ROOT_CUT,
                 ref_base: bool = False,
                 show_leaves: bool = False) -> go.Figure:
    """

    Parameters
    ----------
    ontology: str
        Ontology to use, must be in : [metacyc, ec, chebi]
    metabolic_objects: Collection[str]
        Set of metabolic objects to classify
    reference_set: Collection[str] (optional, default=None)
        Set of reference metabolic objects
    analysis: str (optional, default=topology)
        Analysis mode, must be in : [topology, enrichment]
    output: str (optional, default=None)
        Path to output to save figure
    class_file: str (optional, default=None)
        Path to class ontology file
    names_file: str (optional, default=None)
        Path to EC_ID - EC_NAME association json file
    endpoint_url: str (optional, default=None)
        URL of ChEBI ontology for SPARQL requests
    test: str (optional, default=binomial)
        Type of test if analysis=enrichment, must be in : [binomial, hypergeometric]
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)
    total: bool (optional, default=True)
        True to have branch values proportional of the total parent (may not work in some cases)
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)
    ref_base: bool (optional, default=False)
        True to have the base classes representation of the reference set in the figure.
    show_leaves: bool (optional, default=False)
        True to show input metabolic objets at sunburst leaves

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """

    # METACYC --------------------------------------------------------------------------------------
    if ontology == METACYC:
        if class_file is None:
            class_file = METACYC_FILE
        names = None
        with open(class_file, 'r') as f:
            d_classes_ontology = json.load(f)
    # EC -------------------------------------------------------------------------------------------
    elif ontology == EC:
        if class_file is None:
            class_file = EC_ONTO_FILE
        if names_file is None:
            names_file = EC_NAMES_FILE
        with open(class_file, 'r') as f:
            d_classes_ontology = json.load(f)
        with open(names_file, 'r') as f:
            names = json.load(f)
    # KEGG -----------------------------------------------------------------------------------------
    elif ontology == KEGG:
        if class_file is None:
            class_file = KEGG_ONTO_FILE
        names = None
        with open(class_file, 'r') as f:
            d_classes_ontology = json.load(f)
    # CHEBI ----------------------------------------------------------------------------------------
    elif ontology == CHEBI:
        endpoint_url = CHEBI_URL
        d_classes_ontology = None
        names = None
    # GO -------------------------------------------------------------------------------------------
    elif ontology == GO:
        endpoint_url = GO_URL
        d_classes_ontology = None
        names = None
    # ELSE -----------------------------------------------------------------------------------------
    else:
        raise ValueError(f'ontology parameter must be in : {[METACYC, EC, KEGG, CHEBI, GO]}')
    # WORKFLOW -------------------------------------------------------------------------------------
    return _global_analysis(ontology=ontology, analysis=analysis,
                            metabolic_objects=metabolic_objects, reference_set=reference_set,
                            d_classes_ontology=d_classes_ontology, endpoint_url=endpoint_url,
                            output=output, full=full, names=names, total=total, test=test,
                            root=ROOTS[ontology], root_cut=root_cut, ref_base=ref_base,
                            show_leaves=show_leaves)


# ==================================================================================================
#                                             FUNCTIONS
# ==================================================================================================
def _global_analysis(ontology, analysis, metabolic_objects, reference_set, d_classes_ontology,
                     endpoint_url, output, full, names, total, test, root, root_cut, ref_base,
                     show_leaves):
    obj_all_classes, d_classes_ontology = extract_classes(ontology, metabolic_objects, root,
                                                          d_classes_ontology=d_classes_ontology,
                                                          endpoint_url=endpoint_url)
    classes_abundance = get_classes_abondance(obj_all_classes, show_leaves)

    if output is not None and ontology == METACYC:
        write_met_classes(obj_all_classes, output)

    if reference_set is not None:
        ref_all_classes, d_classes_ontology = extract_classes(ontology, reference_set, root,
                                                              d_classes_ontology=d_classes_ontology,
                                                              endpoint_url=endpoint_url)
    else:
        ref_all_classes = None

    # Enrichment figure
    if analysis == ENRICHMENT_A:
        if ref_all_classes is not None:
            return _enrichment_analysis(ref_all_classes=ref_all_classes,
                                        classes_abundance=classes_abundance,
                                        d_classes_ontology=d_classes_ontology,
                                        output=output, full=full, names=names, total=total,
                                        test=test, root=root, root_cut=root_cut, ref_base=ref_base,
                                        show_leaves=show_leaves)
        else:
            raise AttributeError('Missing reference set parameter')

    # Proportion figure
    elif analysis == TOPOLOGY_A:
        return _topology_analysis(ref_all_classes=ref_all_classes,
                                  classes_abundance=classes_abundance,
                                  d_classes_ontology=d_classes_ontology,
                                  output=output, full=full, names=names, total=total,
                                  root=root, root_cut=root_cut, ref_base=ref_base,
                                  show_leaves=show_leaves)

    else:
        raise ValueError(f'Value of analysis parameter must be in : {[TOPOLOGY_A, ENRICHMENT_A]}')


def _topology_analysis(ref_all_classes, classes_abundance, d_classes_ontology, output, full,
                       names, total, root, root_cut, ref_base, show_leaves) -> go.Figure:
    """ Performs the proportion analysis

    Parameters
    ----------
    ref_all_classes
    classes_abundance
    d_classes_ontology
    output
    full
    names
    total
    root
    root_cut
    ref_base

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    if ref_all_classes is not None:
        ref_classes_abundance = get_classes_abondance(ref_all_classes, show_leaves)
        d_classes_ontology = reduce_d_ontology(d_classes_ontology, ref_classes_abundance)

        if ref_base:
            data = get_fig_parameters(classes_abondance=ref_classes_abundance,
                                      parent_dict=d_classes_ontology,
                                      children_dict=get_children_dict(d_classes_ontology),
                                      root_item=root, subset_abundance=classes_abundance,
                                      full=full, names=names)

        else:
            data = get_fig_parameters(classes_abondance=classes_abundance,
                                      parent_dict=d_classes_ontology,
                                      children_dict=get_children_dict(d_classes_ontology),
                                      root_item=root, full=full, names=names)

        data = get_data_proportion(data, total)
        names = names is not None
        return generate_sunburst_fig(data=data, output=output, analysis=TOPOLOGY_A,
                                     ref_classes_abundance=ref_classes_abundance, names=names,
                                     total=total, root_cut=root_cut, ref_base=ref_base)

    else:
        d_classes_ontology = reduce_d_ontology(d_classes_ontology, classes_abundance)
        data = get_fig_parameters(classes_abondance=classes_abundance,
                                  parent_dict=d_classes_ontology,
                                  children_dict=get_children_dict(d_classes_ontology),
                                  root_item=root, full=full, names=names)
        data = get_data_proportion(data, total)
        names = names is not None
        return generate_sunburst_fig(data=data, output=output, analysis=TOPOLOGY_A,
                                     names=names, total=total, root_cut=root_cut, ref_base=ref_base)


def _enrichment_analysis(ref_all_classes, classes_abundance, d_classes_ontology, output, full,
                         names, total, test, root, root_cut, ref_base, show_leaves) -> go.Figure:
    """ Performs the comparison analysis

    Parameters
    ----------
    ref_all_classes
    classes_abundance
    d_classes_ontology
    output
    full
    names
    total
    test
    root
    root_cut
    ref_base

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    ref_classes_abundance = get_classes_abondance(ref_all_classes, show_leaves)

    if ref_base:
        data = get_fig_parameters(classes_abondance=ref_classes_abundance,
                                  parent_dict=d_classes_ontology,
                                  children_dict=get_children_dict(d_classes_ontology),
                                  root_item=root, subset_abundance=classes_abundance,
                                  full=full, names=names)
    else:
        data = get_fig_parameters(classes_abondance=classes_abundance,
                                  parent_dict=d_classes_ontology,
                                  children_dict=get_children_dict(d_classes_ontology),
                                  root_item=root, full=full, names=names)

    data = get_data_proportion(data, total)
    names = names is not None
    return generate_sunburst_fig(data=data, output=output, analysis=ENRICHMENT_A,
                                 ref_classes_abundance=ref_classes_abundance, test=test,
                                 names=names, total=total, root_cut=root_cut)


def write_met_classes(all_classes, output):
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['Compound', 'Classes', 'Common names', 'MetaCyc link']) + '\n')
        for met, classes, in all_classes.items():
            try:
                # name = ' / '.join(pref.dicOfNode[met].misc['COMMON-NAME'])
                name = ''
            except KeyError:
                name = ''
            link = f'https://metacyc.org/compound?orgid=META&id={met}'
            f.write('\t'.join([met, ', '.join(classes), name, link]) + '\n')
