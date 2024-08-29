import os
import json
from typing import List, Dict
from time import time
import plotly.graph_objects as go

from ontosunburst.ontology import get_abundance_dict, get_classes_abundance, get_classes_scores, \
    extract_classes, reduce_d_ontology, METACYC, CHEBI, EC, GO, KEGG, ROOTS

from ontosunburst.data_table_tree import DataTable, get_name, BINOMIAL_TEST, ROOT_CUT, PATH_UNCUT
from ontosunburst.sunburst_fig import generate_sunburst_fig, TOPOLOGY_A, ENRICHMENT_A

# ==================================================================================================
#                                           CONSTANTS
# ==================================================================================================

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# Dictionary json files
# ---------------------
DEFAULT_FILE = {METACYC: os.path.join(CURRENT_DIR, 'Inputs', 'MetaCyc26_0_classes.json'),
                EC: os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_ontology.json'),
                KEGG: os.path.join(CURRENT_DIR, 'Inputs', 'kegg_onto.json')}
# Names json files
# ----------------
DEFAULT_NAMES = {EC: os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_class_names.json'),
                 METACYC: None, KEGG: None, CHEBI: None, GO: None}
DEFAULT = 'default'
# Sparql URL
# ----------
DEFAULT_URL = {CHEBI: 'http://localhost:3030/chebi/',
               GO: 'http://localhost:3030/go/'}


# ==================================================================================================
#                                            WORKFLOW
# ==================================================================================================

def ontosunburst(interest_set: List[str],
                 ontology: str = None,
                 root: str = None,
                 abundances: List[float] = None,
                 scores: Dict[str, float] = None,
                 reference_set: List[str] = None,
                 ref_abundances: List[float] = None,
                 analysis: str = TOPOLOGY_A,
                 output: str = 'sunburst',
                 write_output: bool = True,
                 class_ontology: str or Dict[str, str] = None,
                 labels: str or Dict[str, str] = DEFAULT,
                 endpoint_url: str = None,
                 test: str = BINOMIAL_TEST,
                 root_cut: str = ROOT_CUT,
                 path_cut: str = PATH_UNCUT,
                 ref_base: bool = False,
                 show_leaves: bool = False,
                 **kwargs) -> go.Figure:
    """ Main function to be called generating the sunburst figure

    Parameters
    ----------
    interest_set: List[str]
        Interest set of concepts to classify
    ontology: str (optional, default=None, values in ['metacyc', 'ec', 'chebi', 'kegg', 'go', None])
        Ontology to use
    root: str (optional, default=None)
        Root item of the ontology.
    abundances: List[str] (optional, default=None)
        Abundance values associated to interest_set list parameter
    scores: Dict[str, float] (optional, default=None)
        Dictionary associating for each ontology ID, its enrichment score. If None enrichment will
        be calculated.
    reference_set: List[str] (optional, default=None)
        Reference set of concepts
    ref_abundances: List[str] (optional, default=None)
        Abundance values associated to reference_set list parameter
    analysis: str (optional, default='topology', values in ['topology', 'enrichment'])
        Analysis mode : topology or enrichment.
    output: str (optional, default='sunburst')
        Path of the output to save figure, if None, outputs will be sunburst.html and sunburst.tsv
        files
    write_output: bool (optional, default=True)
        True to write the html figure and tsv class files, False to only return plotly sunburst
        figure
    class_ontology: str or Dict[str, str] (optional, default=None)
        Class ontology dictionary or json file.
    labels: str or Dict[str, str] (optional, default='default')
        Path to ID-LABELS association json file or ID-LABELS association dictionary or 'default'
        to use default files. If None ontology IDs will be used as labels.
    endpoint_url: str (optional, default=None)
        URL of ChEBI or GO ontology for SPARQL requests
    test: str (optional, default='binomial', values in ['binomial', 'hypergeometric'])
        Type of test if analysis=enrichment, binomial or hypergeometric test.
    root_cut: str (optional, default='cut', values in ['uncut', 'cut', 'total'])
        mode for root cutting (uncut, cut or total)
    path_cut: str (optional, default='uncut', values in ['uncut', 'deeper', 'higher', 'bound'])
        mode for nested path cutting (uncut, deeper, higher or bound)
    ref_base: bool (optional, default=False)
        True to have the base classes representation of the reference set in the figure.
    show_leaves: bool (optional, default=False)
        True to show input metabolic objets at sunburst leaves
    **kwargs

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    start_time = time()
    # LOAD NAMES -----------------------------------------------------------------------------------
    if labels == DEFAULT:
        if ontology is not None:
            labels = DEFAULT_NAMES[ontology]
        else:
            labels = None
    if labels is not None:
        if type(labels) == str:
            with open(labels, 'r') as f:
                names = json.load(f)
        else:
            names = labels
    else:
        names = None
    # DICTIONARY / JSON INPUT ----------------------------------------------------------------------
    if ontology == METACYC or ontology == EC or ontology == KEGG or ontology is None:
        if ontology is None:
            if class_ontology is None:
                raise ValueError('If no default ontology, must fill class_ontology parameter')
            if root is None:
                raise ValueError('If no default ontology, must fill root parameter')
        else:
            if class_ontology is None:
                class_ontology = DEFAULT_FILE[ontology]
        if type(class_ontology) == str:
            with open(class_ontology, 'r') as f:
                class_ontology = json.load(f)
    # SPARQL URL INPUT -----------------------------------------------------------------------------
    elif ontology == CHEBI or ontology == GO:
        if endpoint_url is None:
            endpoint_url = DEFAULT_URL[ontology]
    # ELSE -----------------------------------------------------------------------------------------
    else:
        raise ValueError(f'ontology parameter must be in {[METACYC, EC, KEGG, GO, CHEBI, None]}')
    # GET ROOT -------------------------------------------------------------------------------------
    if ontology is not None:
        root = ROOTS[ontology]
    # WORKFLOW -------------------------------------------------------------------------------------
    fig = _global_analysis(ontology=ontology, analysis=analysis,
                           metabolic_objects=interest_set, abundances=abundances,
                           scores=scores,
                           reference_set=reference_set, ref_abundances=ref_abundances,
                           d_classes_ontology=class_ontology, endpoint_url=endpoint_url,
                           output=output, write_output=write_output, names=names,
                           test=test, root=root, root_cut=root_cut, path_cut=path_cut,
                           ref_base=ref_base, show_leaves=show_leaves, **kwargs)
    end_time = time()
    print(f'Execution time : {end_time - start_time} seconds')
    return fig


# ==================================================================================================
#                                             FUNCTIONS
# ==================================================================================================
def _global_analysis(ontology, analysis, metabolic_objects, abundances, scores, reference_set,
                     ref_abundances, d_classes_ontology, endpoint_url, output, write_output, names,
                     test, root, root_cut, path_cut, ref_base, show_leaves, **kwargs):
    """

    Parameters
    ----------
    ontology
    analysis
    metabolic_objects
    abundances
    reference_set
    ref_abundances
    d_classes_ontology
    endpoint_url
    output
    write_output
    names
    test
    root
    root_cut
    path_cut
    ref_base
    show_leaves
    kwargs

    Returns
    -------

    """
    # EXTRACT CLASSES
    # ----------------------------------------------------------------------------------------------
    obj_all_classes, d_classes_ontology, names = extract_classes(ontology, metabolic_objects, root,
                                                                 d_classes_ontology, endpoint_url,
                                                                 names)
    if not obj_all_classes:
        print('No object classified, passing.')
    if write_output:
        write_met_classes(ontology, obj_all_classes, output, names)

    # ABUNDANCES
    # ----------------------------------------------------------------------------------------------
    # Interest
    abundances_dict = get_abundance_dict(abundances=abundances,
                                         metabolic_objects=metabolic_objects,
                                         ref=False)
    classes_abundance = get_classes_abundance(obj_all_classes, abundances_dict, show_leaves)
    # Reference
    if reference_set is not None:
        ref_abundances_dict = get_abundance_dict(abundances=ref_abundances,
                                                 metabolic_objects=reference_set,
                                                 ref=True)
        ref_all_classes, d_classes_ontology, names = extract_classes(ontology, reference_set, root,
                                                                     d_classes_ontology,
                                                                     endpoint_url, names)
        ref_classes_abundance = get_classes_abundance(ref_all_classes, ref_abundances_dict,
                                                      show_leaves)
    else:
        ref_classes_abundance = None
    # Scores
    classes_scores = None
    if scores is not None:
        classes_scores = get_classes_scores(classes_abundance, scores, root)

    # DATA TABLE
    # ----------------------------------------------------------------------------------------------
    data = DataTable()
    if ref_classes_abundance is not None:
        d_classes_ontology = reduce_d_ontology(d_classes_ontology,
                                               {**ref_classes_abundance, **classes_abundance})
        ref_set = True
        data.fill_parameters(set_abundance=classes_abundance, ref_abundance=ref_classes_abundance,
                             parent_dict=d_classes_ontology, root_item=root, names=names,
                             ref_base=ref_base)
    else:
        ref_set = False
        d_classes_ontology = reduce_d_ontology(d_classes_ontology, classes_abundance)
        data.fill_parameters(set_abundance=classes_abundance, ref_abundance=classes_abundance,
                             parent_dict=d_classes_ontology,  root_item=root, names=names,
                             ref_base=ref_base)
    data.calculate_proportions(ref_base)
    significant = None
    if analysis == ENRICHMENT_A:
        significant = data.make_enrichment_analysis(test, classes_scores)
    data.cut_root(root_cut)
    data.cut_nested_path(path_cut, ref_base)

    # FIGURE
    # ----------------------------------------------------------------------------------------------
    return generate_sunburst_fig(data=data, output=output, analysis=analysis, test=test,
                                 significant=significant, ref_set=ref_set,
                                 write_fig=write_output, **kwargs)


def write_met_classes(ontology: str, all_classes: Dict[str, List[str]], output: str,
                      names: Dict[str, str]):
    """ Writes, for each input class, all its ancestors in a .tsv file.

    Parameters
    ----------
    ontology
    all_classes
    output
    names
    """
    if ontology is None:
        ontology = ''
    links_dict = {METACYC: 'https://metacyc.org/compound?orgid=META&id=',
                  CHEBI: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  EC: 'https://enzyme.expasy.org/EC/',
                  KEGG: 'https://www.genome.jp/entry/',
                  GO: 'https://amigo.geneontology.org/amigo/term/',
                  '': ''}
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['ID', 'Label', 'Classes ID', 'Classes Label', 'Link']) + '\n')
        for met_id, classes_id, in all_classes.items():
            link = links_dict[ontology] + met_id
            met_lab = get_name(met_id, names)
            classes_lab = [get_name(cl, names) for cl in classes_id]
            f.write('\t'.join([met_id, met_lab, ', '.join(classes_id), ', '.join(classes_lab),
                               link]) + '\n')
