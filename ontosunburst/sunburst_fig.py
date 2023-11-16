from typing import List, Dict, Tuple
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import scipy.stats as stats

# CONSTANTS
# ==================================================================================================
# Comparison tests
BINOMIAL_TEST = 'Binomial'
HYPERGEO_TEST = 'Hypergeometric'

# Analysis method
COMPARISON_METHOD = 'comparison'
PROPORTION_METHOD = 'proportion'

MAX_RELATIVE_NB = 1000000

# Root cut
ROOT_CUT = 'cut'
ROOT_TOTAL_CUT = 'total'
ROOT_UNCUT = 'uncut'

# Keys
# ----
IDS = 'ID'
PARENT = 'Parent'
LABEL = 'Label'
COUNT = 'Count'
PROP = 'Proportion'
R_PROP = 'Relative_proportion'
PROP_DIF = 'Proportion_difference'
PVAL = 'Pvalue'


# FUNCTIONS
# ==================================================================================================
def get_fig_parameters(classes_abondance: Dict[str, int], parent_dict: Dict[str, List[str]],
                       children_dict: Dict[str, List[str]], root_item, full: bool = True,
                       names: Dict[str, str] = None) -> Dict[str, List]:
    """ Returns a dictionary of parameters to create the sunburst figure.

    Parameters
    ----------
    classes_abondance: Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes
    children_dict: Dict[str, List[str]]
        Dictionary associating for each class, its children classes
    full: bool (default=True)
        True for a full figure with class duplication if a class has +1 parents.
        False to have a reduced vew with all label appearing only once (1 random parent chosen)
    root_item: str
        Name of the root item of the ontology
    names: Dict[str, str]
        Dictionary associating metabolic object ID to its Name

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    """
    data = {IDS: list(),
            PARENT: list(),
            LABEL: list(),
            COUNT: list()}
    for c_label, c_abundance in classes_abondance.items():
        if c_label != root_item:
            c_parents = parent_dict[c_label]
            if names is None:
                data = add_value_data(data, c_label, c_label, c_abundance, c_parents[0])
            else:
                data = add_value_data(data, c_label, names[c_label], c_abundance, c_parents[0])
            if len(c_parents) > 1 and full:
                for p in c_parents[1:]:
                    suffix = '__' + p
                    c_id = c_label + suffix
                    data = add_value_data(data, c_id, c_label, c_abundance, p)
                    if c_label in children_dict.keys():
                        c_children = children_dict[c_label]
                        for c in c_children:
                            data = add_children(data, suffix, c, c_id, classes_abondance,
                                                children_dict)
        else:
            data = add_value_data(data, c_label, c_label, c_abundance, '')

    return data


def add_value_data(data: Dict[str, List], m_id: str, label: str, value: int, parent: str) -> \
        Dict[str, List]:
    """ Fill the data dictionary for a metabolite class.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    m_id: str
        ID of the metabolite class to add
    label: str
        Label (name) of the metabolite class to add
    value: int
        Abundance value of the metabolite class to add
    parent: str
        Parent metabolite class of the metabolite class to add

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    """
    data[IDS].append(m_id)
    data[LABEL].append(label)
    data[PARENT].append(parent)
    data[COUNT].append(value)
    return data


def add_children(data: Dict[str, List], origin: str, child: str, parent: str,
                 classes_abondance: Dict[str, int], children_dict: Dict[str, List[str]]) \
        -> Dict[str, List]:
    """ Add recursively all children of a given class to the data dictionary.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    origin: str
        Origin of propagation : parent class of parent
    child: str
        Child metabolite class
    parent: str
        Parent metabolite class
    classes_abondance: Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    children_dict: Dict[str, List[str]]
        Dictionary associating for each class, its children classes

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    """
    if child in classes_abondance.keys():
        data[IDS].append(child + origin)
        data[LABEL].append(child)
        data[PARENT].append(parent)
        data[COUNT].append(classes_abondance[child])
        if child in children_dict.keys():
            origin_2 = origin + '__' + child
            cs = children_dict[child]
            for c in cs:
                add_children(data, origin_2, c, child + origin, classes_abondance, children_dict)
    return data


def get_data_proportion(data: Dict[str, List], total: bool) -> Dict[str, List]:
    """ Add a proportion value for color. If total add relative proportion to +1 parent for branch
    value.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    total: bool
        True to have branch values proportional of the total parent

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    """
    # Get total proportion
    max_abondance = np.max(data[COUNT])
    data[PROP] = [x / max_abondance for x in data[COUNT]]
    # Get proportion relative to +1 parent proportion for total branch value
    if total:
        data[R_PROP] = [x for x in data[PROP]]
        p = ''
        data = get_relative_prop(data, p)
        # TODO : IDK WHY IT WORKS ???
        missed = [data[IDS][i] for i in range(len(data[IDS])) if data[R_PROP][i] < 1]
        if missed:
            parents = {data[PARENT][data[IDS].index(m)] for m in missed}
            for p in parents:
                data = get_relative_prop(data, p)
            missed = [data[IDS][i] for i in range(len(data[IDS])) if data[R_PROP][i] < 1]
    return data


def get_relative_prop(data: Dict[str, List], p_id: str):
    """ Get recursively relative proportion of a parent children to itself. Add id to data
    Relative_proportion.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    p_id: str
        ID of the parent

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
            - branch proportion : Relative_prop --> + actual children values
    """
    if p_id == '':
        prop_p = MAX_RELATIVE_NB
        count_p = max(data[COUNT])
    else:
        prop_p = data[R_PROP][data[IDS].index(p_id)]
        count_p = data[COUNT][data[IDS].index(p_id)]
    index_p = [i for i, v in enumerate(data[PARENT]) if v == p_id]
    p_children = [data[IDS][i] for i in index_p]
    count_p_children = [data[COUNT][i] for i in index_p]
    if sum(count_p_children) > count_p:
        total = sum(count_p_children)
    else:
        total = count_p
    for i, c in enumerate(p_children):
        prop_c = int((count_p_children[i] / total) * prop_p)
        data[R_PROP][data[IDS].index(c)] = prop_c
    for c in p_children:
        if c in data[PARENT]:
            data = get_relative_prop(data, c)
    return data


def get_data_prop_diff(data: Dict[str, List], ref_classes_abundance: Dict[str, int]):
    i_max_abondance = np.max(data[COUNT])
    i_prop = [x / i_max_abondance for x in data[COUNT]]
    b_max_abondance = np.max(list(ref_classes_abundance.values()))
    b_prop = [ref_classes_abundance[x] / b_max_abondance if x in ref_classes_abundance.keys() else 0 for
              x in data[LABEL]]
    diff = [i_prop[i] - b_prop[i] for i in range(len(i_prop))]
    data[PROP_DIF] = diff
    return data


def get_data_enrichment_analysis(data: Dict[str, List], ref_classes_abundance: Dict[str, int],
                                 test: str, names: bool) -> Tuple[Dict[str, List], Dict[str, float]]:
    """ Performs statistical tests for enrichment analysis.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    ref_classes_abundance: Dict[str, int]
        Abundances of reference set classes
    test: str
        Type of test Binomial or Hypergeometric
    names: bool
        True if names associated with labels, False otherwise

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
            - branch proportion : Relative_prop
            - p-value : P-values of enrichment analysis
    Dict[str, float]
        Dictionary of significant metabolic object label associated with their p-value
    """
    M = np.max(list(ref_classes_abundance.values()))
    if names:
        m_list = [ref_classes_abundance[x] if x in ref_classes_abundance.keys() else 0 for x in
                  data[IDS]]
    else:
        m_list = [ref_classes_abundance[x] if x in ref_classes_abundance.keys() else 0 for x in
                  data[LABEL]]
    N = np.max(data[COUNT])
    n_list = data[COUNT]
    data[PVAL] = list()
    nb_classes = len(set(data[IDS]))
    significant_representation = dict()
    for i in range(len(m_list)):
        # Binomial Test
        if test == BINOMIAL_TEST:
            p_val = stats.binomtest(n_list[i], N, m_list[i] / M, alternative='two-sided').pvalue
        # Hypergeometric Test
        elif test == HYPERGEO_TEST:
            p_val = stats.hypergeom.sf(n_list[i] - 1, M, m_list[i], N)
        else:
            raise ValueError(f'test parameter must be in : {[BINOMIAL_TEST, HYPERGEO_TEST]}')
        if ((n_list[i] / N) - (m_list[i] / M)) > 0:
            data[PVAL].append(-np.log10(p_val))
        else:
            data[PVAL].append(np.log10(p_val))
        if p_val < 0.05 / nb_classes:
            significant_representation[data[LABEL][i]] = p_val.round(10)
    significant_representation = dict(
        sorted(significant_representation.items(), key=lambda item: item[1]))
    return data, significant_representation


def data_cut_root(data: Dict[str, List], mode: str) -> Dict[str, List]:
    """ Filter data to cut (or not) the root to remove not necessary 100% represented classes.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary of figure parameters
    mode: str
        Mode of root cutting
        - uncut: doesn't cut and keep all nodes from ontology root
        - cut: keep only the lowest level 100% shared node
        - total: remove all 100% shared nodes (produces a pie at center)

    Returns
    -------
    data: Dict[str, List]
        Dictionary of figure parameters with root cut applied
    """
    if mode not in {ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT}:
        raise ValueError(f'Root cutting mode {mode} unknown, '
                         f'must be in {[ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT]}')
    if mode == ROOT_UNCUT:
        return data
    else:
        roots_ind = [i for i in range(len(data[IDS])) if data[R_PROP][i] == MAX_RELATIVE_NB]
        roots = [data[IDS][i] for i in roots_ind]
        for root_id in roots:
            root_ind = data[IDS].index(root_id)
            for v in data.values():
                del v[root_ind]

        if mode == ROOT_TOTAL_CUT:
            data[PARENT] = ['' if x in roots else x for x in data[PARENT]]

    return data


def generate_sunburst_fig(data: Dict[str, List[str or int or float]], output: str = None,
                          sb_type: str = PROPORTION_METHOD, ref_classes_abundance=None,
                          test=BINOMIAL_TEST, names: bool = False, total: bool = True,
                          root_cut: str = ROOT_CUT) -> go.Figure:
    """ Generate a Sunburst figure and save it to output path.

    Parameters
    ----------
    data: Dict[str, List[str or int or float]]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
    output: str (optional, default=None)
        Path to output to save the figure without extension
    sb_type: str (optional, default=Proportion)
        Type of sunburst : Proportion or Comparison
    ref_classes_abundance: Dict[str, int] (optional, default=None)
        Abundances of reference set classes
    test: str (optional, default=Binomial)
        Type of test for enrichment analysis : Binomial or Hypergeometric
    names: bool (optional, default=False)
        True if names associated with labels, False otherwise
    total: bool (optional, default=True)
        True to have branch values proportional of the total parent
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)

    Returns
    -------
    go.Figure
    """
    data = data_cut_root(data, root_cut)
    if total:
        branch_values = 'total'
        values = data[R_PROP]
    else:
        branch_values = 'remainder'
        values = data[COUNT]
    if sb_type == PROPORTION_METHOD:
        fig = go.Figure(go.Sunburst(labels=data[LABEL], parents=data[PARENT], values=values,
                                    ids=data[IDS],
                                    hoverinfo='label+text', maxdepth=7,
                                    branchvalues=branch_values,
                                    hovertext=[f'Count: {data[COUNT][i]}<br>'
                                               f'Proportion: {round(data[PROP][i]*100, 2)}%<br>'
                                               f'ID: {data[IDS][i]}'
                                               for i in range(len(data[PROP]))],
                                    marker=dict(colors=data[COUNT],
                                                colorscale=px.colors.sequential.Viridis,
                                                cmid=0.5 * max(data[COUNT]), showscale=True,
                                                colorbar=dict(title=dict(text='Count')))))
        fig.update_layout(title=dict(text='Proportion of classes', x=0.5, xanchor='center'))
    elif sb_type == COMPARISON_METHOD:
        data, signif = get_data_enrichment_analysis(data, ref_classes_abundance, test, names)
        # m = np.mean(data[PROP_DIF])
        fig = make_subplots(rows=1, cols=2,
                            column_widths=[0.3, 0.7],
                            vertical_spacing=0.03,
                            subplot_titles=('Significant p-values',
                                            'Classes enrichment representation'),
                            specs=[[{'type': 'table'}, {'type': 'sunburst'}]])

        fig.add_trace(go.Sunburst(labels=data[LABEL], parents=data[PARENT],
                                  values=values, ids=data[IDS],
                                  hovertext=[f'P value: {10 ** (-data[PVAL][i])}<br>'
                                             f'Count: {data[COUNT][i]}<br>'
                                             f'Proportion: {round(data[PROP][i]*100, 2)}%<br>'
                                             f'ID: {data[IDS][i]}'
                                             if data[PVAL][i] > 0 else
                                             f'P value: {10 ** data[PVAL][i]}<br>'
                                             f'Count: {data[COUNT][i]}<br>'
                                             f'Proportion: {round(data[PROP][i]*100, 2)}%<br>'
                                             f'ID: {data[IDS][i]}'
                                             for i in range(len(data[PVAL]))],
                                  hoverinfo='label+text', maxdepth=7,
                                  branchvalues=branch_values,
                                  marker=dict(colors=data[PVAL],
                                              colorscale=px.colors.diverging.RdBu,
                                              cmid=0, cmax=10.0, cmin=-10.0, showscale=True,
                                              colorbar=dict(title=dict(text='Log10(p-value)')))),
                      row=1, col=2)

        fig.add_trace(go.Table(header=dict(values=['Metabolite', f'{test} test P-value'],
                                           fill=dict(color='#666666'), height=40,
                                           font=dict(size=20)),
                               cells=dict(values=[list(signif.keys()), list(signif.values())],
                                          fill=dict(color='#777777'), height=35,
                                          font=dict(size=16))),
                      row=1, col=1)
    else:
        raise ValueError('Wrong type input')
    fig.update_layout(paper_bgcolor="#888888", font_color='#111111', font_size=20)
    fig.update_annotations(font_size=28)
    if output is not None:
        fig.write_html(f'{output}.html')
    return fig
