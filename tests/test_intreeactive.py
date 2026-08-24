from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo

from src.intreeactive import intreeactive

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# Constants and fixtures:
node_list = ['A', 'F', 'B', 'E', 'C', 'D']
id_column = 'ID'

root = Path(__file__).parents[1]
path_to_metadata = Path(root / 'tests/test_snp_metadata.csv')


@pytest.fixture
def test_tree():
    """
    :return: tree, Phylo object - 4 leaves, 6 nodes, rooted and ladderized.
    """
    handle = StringIO('(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F')
    tree = Phylo.read(handle, 'newick')
    tree.root_with_outgroup({'name': 'A'})
    tree.ladderize(reverse=False)
    # Phylo.draw_ascii(tree)
    """ Y on the left, X corresponds to the branch lengths.
    1     , A
         _|
    2     |      __0.2__ B
          |_0.1_|F
    3           |              ___0.3___ C
                |_____0.5_____|E
    4                         |____0.4____ D
    """
    return tree


@pytest.fixture
def x_y_coords(test_tree) -> tuple[dict, dict]:
    """
    Get xcoords and ycoords dictionaries, using the class instances as keys and coords as vals.
    """
    xcoords = test_tree.depths()

    nodes = list(test_tree.find_clades())
    raw_coords = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 3.5, 'F': 2.75}
    ycoords = {}
    for clade in nodes:
        if clade.name is not None:
            ycoords[clade] = raw_coords[str(clade.name)]
        else:
            ycoords[clade] = 1.875

    return xcoords, ycoords


@pytest.fixture
def snp_dist_matrix() -> pd.DataFrame:
    """
    Create a snp distance matrix for 5 samples - noting 'O' is not in the tree.
    """
    snp_dists = {
        'A': [0, 3, 9, 10, 12],
        'B': [3, 0, 6, 7, 12],
        'C': [9, 6, 0, 1, 2],
        'D': [10, 7, 1, 0, 1],
        'O': [12, 12, 2, 1, 0],
    }
    snp_dists_matrix = pd.DataFrame.from_dict(snp_dists)
    snp_dists_matrix.index = ['A', 'B', 'C', 'D', 'O']
    return snp_dists_matrix


@pytest.fixture
def metadata_with_neighbours() -> pd.DataFrame:
    """
    Create metadata to accompany the tree. Note that 'O' does not exist in the tree, and internal nodes F and E have no
    corresponding metadata.
    """
    metadata = {
        'A': ['A', 'Ashley', 'Anon', '2024-01-01', True, 'B=3'],
        'B': ['B', 'Barbara', 'Bettername', '2024-01-02', True, 'A=3'],
        'C': ['C', 'Charlie', 'Choo', '2024-01-03', False, 'D=1'],
        'D': ['D', 'Deborah', 'Doe', '2024-01-04', False, 'C=1<br>O=1'],
        'O': ['O', 'Otto', 'Otherington', '2024-01-05', True, 'D=1'],
    }
    metadata_df = pd.DataFrame.from_dict(
        metadata, orient='index', columns=['ID', 'name', 'last', 'date', 'somebool', 'Nearest_neighbour']
    ).reset_index(drop=True)
    return metadata_df


# Helper functions
def get_coords(coords_dict: dict, clade_name: str) -> float:
    """
    Nicely return the value from the dictionary where keys are the clade class instance.
    """
    key_to_get = [clade for clade in coords_dict if clade.name == clade_name]
    assert len(key_to_get) == 1, 'More than one node with the same clade name.'
    return coords_dict[key_to_get[0]]


# Tests
def test_exit_if_sample_in_tree_but_not_in_snpdist_matrix(test_tree, metadata_with_neighbours, snp_dist_matrix):
    snp_dist_matrix = snp_dist_matrix.drop('A')
    snp_dist_matrix = snp_dist_matrix.drop('A', axis=1)
    print('\n  If sample is not in the snp distance matrix but is in the tree, it should exit.')
    with pytest.raises(SystemExit):
        intreeactive.check_ids(
            tree=test_tree, metadata=metadata_with_neighbours, id_column=id_column, snp_dists=snp_dist_matrix
        )


def test_exit_if_sample_in_tree_but_not_in_metadata(test_tree, metadata_with_neighbours, snp_dist_matrix):
    metadata_with_neighbours = metadata_with_neighbours[metadata_with_neighbours['ID'] != 'A']
    print('\n  If sample is not in the metadata but is in the tree, it should exit.')
    with pytest.raises(SystemExit):
        intreeactive.check_ids(
            tree=test_tree, metadata=metadata_with_neighbours, id_column=id_column, snp_dists=snp_dist_matrix
        )


def test_check_ids_occur_where_needed(test_tree, metadata_with_neighbours, snp_dist_matrix):
    outcome = intreeactive.check_ids(
        tree=test_tree, metadata=metadata_with_neighbours, id_column=id_column, snp_dists=snp_dist_matrix
    )
    print(
        f'\n  Sample "O" is in the metadata but not in the tree, so should be dropped from the metadata: \n'
        f'{metadata_with_neighbours} \n'
        f'{outcome}\n'
    )
    assert outcome.size == 24
    assert 'O' not in outcome['ID'].to_list()


def test_get_x_coordinates(test_tree):
    """
    Test function that gets x-coordinates from the phylo tree object.
    Expect dict with clade class instance as key and coordinate (floating point) as value.
    """
    xcoords = intreeactive.get_x_coordinates(test_tree)
    # dict format:
    # key = Clade(branch_length, name), val = x coord
    print(
        f"\n   Expect: Clade(branch_length=0.0, name='A'), val = 0.1)\n"
        f'  got: {[(clade, val) for clade, val in xcoords.items() if clade.name == "A"]}'
    )
    assert len(xcoords.keys()) == 7
    assert round(get_coords(xcoords, 'B'), 1) == 0.3
    assert round(get_coords(xcoords, 'D'), 1) == 1.0
    assert round(get_coords(xcoords, 'A'), 1) == 0.0


def test_get_y_coordinates(test_tree):
    """
    Test function that gets y-coordinates from the phylo tree object.
    Expect dict with clade class instance as key and coordinate (floating point) as value. In this case, y-coordinates
    are multiplied by the constant 1 for easier calculation.
    """
    ycoords = intreeactive.get_y_coordinates(test_tree, dist=1)
    print(
        f"\n   Expect: Clade(branch_length=0.0, name='A'): 1 \n"
        f'   got: {[(clade, val) for clade, val in ycoords.items() if clade.name == "A"]}'
    )
    assert len(ycoords.keys()) == 7  # There are 6 nodes, plus the tree object is wrapped in a Clade()
    assert len(test_tree.get_terminals()) == 4  # Only A, B, C, D are leaves, E and F are internal nodes.
    assert get_coords(ycoords, 'B') == 2
    assert get_coords(ycoords, 'D') == 4
    # Check the calc_row function is working:
    assert get_coords(ycoords, 'F') == 2.75
    assert get_coords(ycoords, 'E') == 3.5


def test_get_y_coordinates_scaling(test_tree):
    """
    Test function that gets y-coordinates from the phylo tree object.
    Expect dict with clade class instance as key and coordinate (floating point) as value. In this case, y-coordinates
    are multiplied by the constant 1 for easier calculation.
    """
    ycoords = intreeactive.get_y_coordinates(test_tree, dist=10)
    print(
        f"\n   Expect: Clade(branch_length=0.0, name='A'): 10 \n"
        f'   got: {[(clade, val) for clade, val in ycoords.items() if clade.name == "A"]}'
    )
    assert len(ycoords.keys()) == 7  # There are 6 nodes, plus the tree object is wrapped in a Clade()
    assert len(test_tree.get_terminals()) == 4  # Only A, B, C, D are leaves, E and F are internal nodes.
    assert get_coords(ycoords, 'B') == 20
    assert get_coords(ycoords, 'D') == 40
    # Check the calc_row function is working:
    assert get_coords(ycoords, 'F') == 27.5
    assert get_coords(ycoords, 'E') == 35


def test_get_clade_lines(test_tree, x_y_coords):
    """
    Test function that draws the all the horizontal and vertical lines in the tree.
    Tree looks roughly like this, with Y axis on the left, X corresponds to the branch lengths.
    1     , A
         _|
    2     |      __0.2__ B
          |_0.1_|F
    3           |              ___0.3___ C
                |_____0.5_____|E
    4                         |____0.4____ D
    """
    xcoords, ycoords = x_y_coords
    test_tree_line_shapes = []
    intreeactive.draw_clade(
        clade=test_tree.root,
        x_coords=xcoords,
        y_coords=ycoords,
        line_shapes=test_tree_line_shapes,
        line_colour='rgb(25,25,25)',
        line_width=1,
    )

    # Horizontal lines:
    horizontal_line_parent_to_f = {'x0': 0, 'x1': 0.1, 'y0': 2.75, 'y1': 2.75}
    print(
        """\n   Expect horizontal line (*) between parental node to F internal node:
           , A
          _|
           |      _______ B
           |__*__|F
                 |
          to be plotted as:""",
        horizontal_line_parent_to_f,
    )
    assert any(
        horizontal_line_parent_to_f.items() <= line_shape_dict.items() for line_shape_dict in test_tree_line_shapes
    )

    horizontal_line_parent_to_d = {'x0': 0.6, 'x1': 1.0, 'y0': 4, 'y1': 4}
    print(
        """\n   Expect horizontal line (*) between parental node to leaf D to be plotted as:
           _________ C
          |
          |_____*_____ D
          """,
        horizontal_line_parent_to_d,
    )
    assert any(
        horizontal_line_parent_to_d.items() <= line_shape_dict.items() for line_shape_dict in test_tree_line_shapes
    )

    # Vertical lines:
    between_c_and_d = {'x0': 0.6, 'x1': 0.6, 'y0': 4, 'y1': 3}
    print(
        """\n   Expect vertical line (*) plotted here between leaves C and D:'
           _________ C
          *
          |___________ D
          to have these coordinates:""",
        between_c_and_d,
    )
    assert any(between_c_and_d.items() <= line_shape_dict.items() for line_shape_dict in test_tree_line_shapes)

    between_a_and_f = {'x0': 0, 'x1': 0, 'y0': 2.75, 'y1': 1}
    print(
        """\n   Expected vertical line (*) plotted here between leaf A and internal node F:
           , A
          _|
           *      _______ B
           |_____|F
                 |
          to have these coordinates:""",
        between_a_and_f,
    )
    assert any(between_a_and_f.items() <= line_shape_dict.items() for line_shape_dict in test_tree_line_shapes)


@pytest.mark.parametrize('leaf,expect', [('A', ['B=3']), ('C', ['D=1']), ('D', ['C=1', 'O=1'])])
def test_get_nearest_neighbours(snp_dist_matrix, leaf, expect):
    """
    Test function that gets the nearest neighbours using the snp distance matrix - ignoring itself, find which other
    sample(s) has the lowest SNP distance.
    """
    actual = intreeactive.get_nearest_neighbours(snp_dist_matrix, leaf, do_join=False)
    print(f'\n   Return nearest neighbour and SNP distance (that is not itself). \nExpect: {expect}.actual: {actual}')
    assert actual == expect


def test_missing_snp_dist_entry(snp_dist_matrix):
    """
    what if a sample is present in the metadata but not in the snp distance matrix?
    Noting that if a sample is in the tree but not in the metadata or not in the snp-distance matrix, the whole thing
    falls over earlier than this. If there are samples in the metadata not in the tree, these should be dropped so this
    should never trigger anyway.
    """
    actual = intreeactive.get_nearest_neighbours(snp_dist_matrix, "doesn't exist", do_join=False)
    expect = ''
    print(
        f'\n If ID is in the metadata but is not in the snp distance matrix: \n'
        f" Expect: empty string '{expect}',\n got: {actual, type(actual)}"
    )
    assert actual == expect


def test_make_hover_text(metadata_with_neighbours):
    """
    Test the function that makes hover text. Hover text is given to plotly as a list, where the index in the list
    corresponds to the nodes in a particular order.
    A list of nodes is given and a list of info is returned, where the info includes all the metadata available for that
    sample, separated by <br>.
    """
    hover_text = intreeactive.make_hover_text(metadata_with_neighbours, id_column, node_list)

    a_hover_text = 'ID: A<br>name: Ashley'
    b_hover_text = 'Nearest_neighbour: A=3'
    c_hover_text = 'date: 2024-01-03<br>somebool: False'
    print(f"\n   Expect list of hover text, with metadata separated by <br>, \nFor example, 'A' got: {hover_text['A']}")

    assert a_hover_text in str(hover_text['A'])
    assert b_hover_text in hover_text['B']
    assert c_hover_text in hover_text['C']
    assert hover_text['F'] == 'F'


def test_get_colourings(metadata_with_neighbours):
    """
    Test that list of colours are created for leaves, with internal nodes coloured black rgb(0,0,0,1.0).
    """
    internal_node_colour = 'rgb(0,0,0)'
    colourings = intreeactive.get_colourings(
        metadata_with_neighbours,
        id_column=id_column,
        category=id_column,
        node_list=node_list,
        intermediate_node_colour=internal_node_colour,
    )
    print(colourings)

    assert colourings['F'] == internal_node_colour
    assert colourings['E'] == internal_node_colour
    assert colourings['A'] != internal_node_colour


def test_get_colourings_missing_data(metadata_with_neighbours):
    """
    Test that list of colours are created for leaves even if data is missing, NoneType or float (NaN), with internal
    nodes coloured black rgb(0,0,0).
    """
    internal_node_colour = 'rgb(0,0,0)'
    metadata_with_neighbours.at[1, 'name'] = np.NaN
    metadata_with_neighbours.at[2, 'name'] = None
    metadata_with_neighbours.at[3, 'name'] = ''

    
    colourings = intreeactive.get_colourings(
        metadata_with_neighbours,
        id_column=id_column,
        category='name',
        node_list=node_list,
        intermediate_node_colour=internal_node_colour,
    )

    assert (len(set(colourings.values()))) == 5, 'Expecting 5 different colours'


def test_get_colourings_numerical_data(metadata_with_neighbours):
    """
    Test that list of colours are created for numerical data, with internal nodes coloured black rgb(0,0,0).
    """
    internal_node_colour = 'rgb(0,0,0)'
    metadata_with_neighbours['nums'] = pd.Series([300, 2, 22, 2, 100])
    colourings = intreeactive.get_colourings(
        metadata_with_neighbours,
        id_column=id_column,
        category='nums',
        node_list=node_list,
        intermediate_node_colour=internal_node_colour,
    )
    assert (len(set(colourings.values()))) == 4, 'Expecting 4 different colours'
    assert colourings['E'] == internal_node_colour and colourings['F'] == internal_node_colour


def test_get_colourings_mixed_data(metadata_with_neighbours):
    """
    Test that list of colours are created for mixed data, with internal nodes coloured black rgb(0,0,0).
    """
    internal_node_colour = 'rgb(0,0,0)'
    metadata_with_neighbours['mixed'] = pd.Series([None, 2.2, 22, '4', ''])
    colourings = intreeactive.get_colourings(
        metadata_with_neighbours,
        id_column=id_column,
        category='mixed',
        node_list=node_list,
        intermediate_node_colour=internal_node_colour,
    )
    
    assert (len(set(colourings.values()))) == 5, f'Expecting 5 different colours, got {len(set(colourings))}'


def test_continuous_colours(metadata_with_neighbours):
    """
    This tests that dates get a continuous gradient.
    The gradient scale runs from 0-1, royal blue (rgb(0, 0, 131)) - maroon red (rgb(128, 0, 0)).
    The newest date (O in metadata) should be blue, the oldest date should be red (A in metadata).
    """
    node_list.append('O')
    internal_node_colour = 'rgb(0, 0, 0)'
    colours = intreeactive.get_continuous_colourings(
        metadata_with_neighbours,
        id_column=id_column,
        date_category='date',
        node_list=node_list,
        intermediate_node_colour=internal_node_colour,
    )
    print(colours)
    print(
        f'\n Colours are assigned as a gradient for the dates, such that:\n'
        f'    - no dates are given rgb(0, 0, 0),\n'
        f'    - oldest is given royal blue rgb(128, 0, 0),\n'
        f'    - newest is given maroon red rgb(0, 0, 131).\n'
        f'For the metadata, the colour assigned for the date is: \n'
        f'{dict(zip(node_list, colours, strict=False))}'
    )
    assert colours['A'] == 'rgb(128, 0, 0)'
    assert colours['O'] == 'rgb(0, 0, 131)'
    assert colours['E'] == 'rgb(0, 0, 0)' and colours['F'] == 'rgb(0, 0, 0)'


def test_all_same_date():
    """
    Test what happens when all the metadata have the same date.
    All the samples should be assigned maroon nodes (rgb(0, 0, 131), with grey internal nodes.
    """
    metadata_same_date = {
        'A': ['A', '01-01-2021'],
        'B': ['B', '01-01-2021'],
        'C': ['C', '01-01-2021'],
        'D': ['D', '01-01-2021'],
    }
    metadata_same_date = pd.DataFrame.from_dict(metadata_same_date, orient='index', columns=['ID', 'date']).reset_index(
        drop=True
    )
    nodes_test_all_same_date = ['internal1', 'A', 'B', 'C', 'D', 'internal2']
    internal_node_colour = 'rgb(128, 128, 128)'
    colours = intreeactive.get_continuous_colourings(
        metadata_same_date,
        id_column='ID',
        date_category='date',
        node_list=nodes_test_all_same_date,
        intermediate_node_colour=internal_node_colour,
    )
    assert colours['internal1'] and colours['internal2'] == 'rgb(128, 128, 128)'  # these are internal nodes
    assert len([name for name in colours if colours[name] == 'rgb(0, 0, 131)']) == 4
    


def test_gradient_colours_with_missing_dates():
    """
    Check that entries without dates (empty string, NoneType or NaN) are assigned black.
    Check that entries with incomplete dates are assigned to the 1st of the year
    Check that entries with the same dates are assigned the same colour.
    """
    metadata_missing_dates = {
        'A': ['A', '2021-01-01'],
        'B': ['B', None],
        'C': ['C', '2021-01-02'],
        'D': ['D', '2021-01'],
        'E': ['E', '2021'],
        'F': ['F', np.NaN],
        'G': ['G', ''],
    }
    metadata_missing_dates = pd.DataFrame.from_dict(
        metadata_missing_dates, orient='index', columns=['ID', 'date']
    ).reset_index(drop=True)
    nodes_test_missing_dates = ['internal1', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'internal2']
    internal_node_colour = 'rgb(128, 128, 128)'
    colours = intreeactive.get_continuous_colourings(
        metadata_missing_dates,
        id_column='ID',
        date_category='date',
        node_list=nodes_test_missing_dates,
        intermediate_node_colour=internal_node_colour,
    )
    
    print(
        f'\n Colours are assigned as a gradient for the dates, such that:\n'
        f'    - no dates are given rgb(0, 0, 0),\n'
        f'    - oldest is given royal blue rgb(128, 0, 0),\n'
        f'    - newest is given maroon red rgb(0, 0, 131).\n'
        f'For the metadata, the colour assigned for the date is: \n'
        f'{metadata_missing_dates})'
    )
    assert colours['internal1'] and colours['internal2'] == 'rgb(128, 128, 128)'  # these are internal nodes
    assert colours['B'] == 'rgb(0, 0, 0)'  # these are the empty dates
    assert colours['F'] == 'rgb(0, 0, 0)'
    assert colours['G'] == 'rgb(0, 0, 0)'


def test_gradient_colours_with_all_same_or_missing():
    """
    If dates are either the same one date or empty, those with a date are maroon (rgb(0, 0, 131)) and those without
    are black (rgb(0, 0, 0)).
    """
    metadata_mixed_dates = {
        'A': ['A', '2021-01-01'],
        'B': ['B', '2021-01-01'],
        'C': ['C', '2021-01-01'],
        'D': ['D', '2021-01-01'],
        'E': ['E', ''],
        'F': ['F', None],
        'G': ['G', '2021-01-01'],
    }
    metadata_mixed_dates_df = pd.DataFrame.from_dict(
        metadata_mixed_dates, orient='index', columns=['ID', 'date']
    ).reset_index(drop=True)

    nodes_test_mixed_dates = ['internal1', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'internal2']
    internal_node_colour = 'rgb(128, 128, 128)'
    colours = intreeactive.get_continuous_colourings(
        metadata_mixed_dates_df,
        id_column='ID',
        date_category='date',
        node_list=nodes_test_mixed_dates,
        intermediate_node_colour=internal_node_colour,
    )
    print(
        f'\n If date is either the same one date it is maroon (rgb(0, 0, 131)) or empty it is black (rgb(0, 0, 0)).\n'
        f'{metadata_mixed_dates_df})'
    )
    assert colours['internal1'] and colours['internal2'] == 'rgb(128, 128, 128)'
    assert colours['E'] == 'rgb(0, 0, 0)' and colours['F'] == 'rgb(0, 0, 0)'
    assert len([name for name in colours if colours[name] == 'rgb(0, 0, 131)']) == 5


def test_read_in_metadata():
    print('Check metadata is parsed correctly when not specifying id column.')
    df, id_col = intreeactive.read_in_metadata(path_to_metadata)
    assert id_col == 'ID'
    assert 'ID' in df.columns.values


def test_read_in_metadata_id_col():
    print('Check metadata is parsed correctly when specifying id column.')
    df, id_col = intreeactive.read_in_metadata(path_to_metadata, 'ID_col')
    assert id_col == 'ID'
    assert 'ID' in df.columns.values


def test_read_in_metadata_id_wrong():
    print('Check metadata is parsed correctly when specifying incorrect id column.')
    df, id_col = intreeactive.read_in_metadata(path_to_metadata, 'id_col')
    assert id_col == 'ID'
    assert 'ID' in df.columns.values


def test_create_tree(test_tree, metadata_with_neighbours, snp_dist_matrix, tmp_path):
    intreeactive.write_interactive_tree(
        tree=test_tree,
        output_name=Path(tmp_path / 'test_tree.html'),
        metadata=metadata_with_neighbours,
        id_column='ID',
        snp_distance_matrix=snp_dist_matrix,
        exclude_internal_nodes=False,
        title='Test Tree',
    )

def test_create_tree_no_internal_nodes(test_tree, metadata_with_neighbours, snp_dist_matrix, tmp_path):
    intreeactive.write_interactive_tree(
        tree=test_tree,
        output_name=Path(tmp_path / 'test_tree.html'),
        metadata=metadata_with_neighbours,
        id_column='ID',
        snp_distance_matrix=snp_dist_matrix,
        exclude_internal_nodes=True,
        title='Test Tree',
    )


# def test_make_12snp_clusters():
#     from scipy.cluster.hierarchy import linkage, fcluster
#     from scipy.spatial.distance import squareform
#
#     # Sample SNP distance matrix (replace with your actual data)
#     dist_matrix = pd.read_csv('../example/tb_in_middle_earth_snpdists_matrix.txt', sep=None, engine='python',
#                               header=0, index_col=0)
#     # Assume `dist_matrix` is your symmetric pandas DataFrame with IDs as index and columns
#     ids = dist_matrix.index.tolist()
#     # Convert the distance matrix to condensed form
#     condensed_matrix = squareform(dist_matrix.values)
#     # Perform single linkage clustering
#     linkage_matrix = linkage(condensed_matrix, method='single')
#     # Set a threshold for clustering - 12 SNPs
#     threshold = 12
#     # Get cluster assignments based on the threshold
#     cluster_labels = fcluster(linkage_matrix, t=threshold, criterion='distance')
#     clusters = pd.DataFrame({'ID': ids, 'Cluster': cluster_labels})
#     # Print the cluster assignments
#     cluster_1 = clusters[clusters['Cluster'] == 1]['ID'].values
#     # Get a reduced filtered matrix for cluster 1
#     filtered_dist_matrix = dist_matrix.loc[cluster_1, cluster_1]
#     print(filtered_dist_matrix)
