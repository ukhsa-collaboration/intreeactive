import os.path
import json
import datetime
import base64
import sys
import plotly.colors

from importlib.resources import files, as_file
from pathlib import Path
from Bio import Phylo
from typing import List, Optional, Literal
import pandas as pd
from bs4 import BeautifulSoup as bs

import plotly.express as px
import plotly.graph_objects as go

# Set up html_res path:
html_res = files('html_res')


def read_in_tree(*,
                 path_to_tree: str | os.PathLike,
                 outgroup: str = None,
                 tree_format: Literal["newick", "nexus", "nexml", "phyloxml", "cdao"] = None):
    """
    Parse tree using Bio.Phylo.read.
    :param path_to_tree: path to tree file.
    :param outgroup: outgroup to root the tree on (str, default = None).
    :param tree_format: format of tree, should be newick, nexus, phyloxm, nexml, cdao (str, default = newick).
    :return: tree - class instance of Bio.Phylo.<tree_format>.Tree. For the default Newick tree, this returns
    Bio.Phylo.Newick.Tree object.
    """
    if not tree_format:
        tree_format = 'newick'
    tree = Phylo.read(path_to_tree, tree_format)
    if outgroup:
        tree.root_with_outgroup({'name': outgroup})
    tree.ladderize(reverse=True)
    return tree


def read_in_metadata(path_to_metadata: str | os.PathLike,
                     id_column: str = None) -> tuple[pd.DataFrame, str]:
    """
    Read in metadata and read with Pandas - all cells should be strings.
    Change the id_column to "ID" and make sure there are no other ID columns.
    """
    metadata_df = pd.read_csv(path_to_metadata, sep=None, engine='python', encoding='windows-1252')
    metadata_df = metadata_df.astype(str)
    # If id_column supplied and present in the dataframe, use that, else use the first column.
    if not id_column or id_column not in metadata_df.columns.values:
        id_column = metadata_df.columns.values[0]
    # Check if there are any other columns called "ID":
    if "ID" in metadata_df.columns.values and metadata_df.columns.values[0] != "ID":
        metadata_df.rename(columns={"ID": "other_id_x"}, inplace=True)
    new_id_column = "ID"
    metadata_df.rename(columns={id_column: new_id_column}, inplace=True)
    return metadata_df, new_id_column


def read_in_snp_dist_matrix(path_to_snp_dists: str | os.PathLike) -> pd.DataFrame:
    """
    Read in snp distance matrix using pandas. The row names and column names should match, and should be parsed as
    strings.
    :param path_to_snp_dists: string or path to the snp distance matrix.
    """
    snpdist_matrix = pd.read_csv(path_to_snp_dists, sep=None, engine='python', index_col=0)
    # Assert that the order of the columns is identical to the order of the rows that way we know we only need to store
    # the index once
    assert list(snpdist_matrix) == list(snpdist_matrix.columns)
    return snpdist_matrix


def check_ids(*,
              tree,
              metadata: pd.DataFrame,
              id_column: str = 'ID',
              snp_dists: pd.DataFrame,
              ignore_ids: Optional[list | str] = None) -> pd.DataFrame:
    """
    If a sample ID occurs in the tree, it must have metadata and snp distances (unless it is an ignored_id
    like the outgroup)
    If a sample ID occurs in the metadata and not in the tree, irrelevant of existence in the snp-dist matrix,
    it can be dropped (for performance later).
    If a sample ID occurs in the snp_dists and not in the tree, it can be ignored.

    :param tree: Bio.Phylo.Newick.Tree object
    :param metadata: a dataframe with metadata - all cells should be strings.
    :param snp_dists: a dataframe matrix of the snp distances. The row names and column names should match, and should
    be parsed as strings.
    :param id_column: string, the name of the ID column used to link up the SNP distances.
    :param ignore_ids: optional; a list of strings or a string of IDs that should be ignored when checking presence of
    an ID in the tree, metadata and snp distance matrix, for example a reference or outgroup.

    :returns: metadata dataframe if changes made, else empty dataframe.
    """
    samples_in_tree = [leaf.name for leaf in tree.get_terminals()]
    samples_in_metadata = metadata[id_column].to_list()
    samples_in_snpdists = list(snp_dists)
    # Check all samples in tree are in metadata and snp dists, if not, exit.
    for sample in samples_in_tree:
        if ignore_ids:
            # Handle whether ignore ids is a list of samples or one sample:
            if isinstance(ignore_ids, list):
                if sample in ignore_ids:
                    print(f"Allowing {sample}")
                    continue
            elif isinstance(ignore_ids, str):
                if sample == ignore_ids:
                    print(f"Allowing {sample}")
                    continue
        # If sample is in tree but not in snpdists - this is critical failure - exit:
        if sample not in samples_in_snpdists:
            sys.exit(f"Error: \n "
                     f"Sample ID {sample} occurs in tree but not in snp distance matrix. \n"
                     f"Exiting...")
        # If sample is in tree but not in metadata - this is critical failure - exit
        elif sample not in samples_in_metadata:
            sys.exit(f"Error: \n "
                     f"Sample ID {sample} occurs in tree but not in the metadata. \n"
                     f"Exiting...")
    # If a sample is in the metadata but not in the tree, drop sample:
    for sample in samples_in_metadata:
        # Again, check for ignore ids:
        if ignore_ids:
            if isinstance(ignore_ids, list):
                if sample in ignore_ids:
                    print(f"Allowing {sample}")
                    continue
            elif isinstance(ignore_ids, str):
                if sample == ignore_ids:
                    print(f"Allowing {sample}")
                    continue
        if sample not in samples_in_tree:
            metadata = metadata[metadata[id_column] != sample]
            return metadata
    return pd.DataFrame()


# https://github.com/empet/Phylogenetic-trees - source for functions to parse the tree and format into plotly-ready
# structures. These functions are based on the same named `get_x_coordinates()`, `get_y_coordinates()` functions in
# Biopython: https://github.com/biopython/biopython/blob/master/Bio/Phylo/_utils.py.
# (They assign cartesian coordinates to the tree nodes.)

def get_x_coordinates(input_tree) -> dict:
    """
    Create dictionary with clades (node in tree) as key and their x-coordinate as value. If no branch lengths in tree,
    use unit branch lengths.
    :param input_tree: phylo object
    :return: dict of clade: x-coord.
    """
    xcoords = input_tree.depths()

    # If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = input_tree.depths(unit_branch_lengths=True)

    return xcoords


def get_y_coordinates(input_tree, dist=0.1) -> dict:
    """
    Create dictionary with clades (node in tree) as key and their y-coordinate as value. Y-coordinates are multiple of
    dist - multiplying by a constant to give constant vertical distance between two connected leaves.
    :param input_tree: phylo object.
    :param dist: constant, value to ensure leaves are constant width apart.
    :return: a dict of clade: y-coord.
    """
    _max_height = input_tree.count_terminals()  # Counts the number of tree leaves.
    ycoords = dict((leaf, (_max_height - ind) * dist) for ind, leaf in enumerate(reversed(input_tree.get_terminals())))

    def calc_row(clade):
        """ Recursively add all parental nodes to terminal leaves."""
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] +
                          ycoords[clade.clades[-1]]) / 2

    if input_tree.root.clades:
        calc_row(input_tree.root)

    return ycoords


def get_clade_lines(orientation: str = 'horizontal', y_curr: str = "0", x_start: str = "0", x_curr: str = "0",
                    y_bot: str = "0", y_top: str = "0", line_colour: str = 'rgb(25,25,25)',
                    line_width: float = 0.5) -> dict:
    """
    Define a Plotly shape of type 'line', for each branch
    :param orientation: horizontal or vertical.
    :param y_curr: current y coord
    :param x_start: starting x coord
    :param x_curr: current x coord
    :param y_bot: bottom y coord
    :param y_top: top y coord
    :param line_colour: line colour, default is black.
    :param line_width: line width, default is 0.5.
    :return: dict of information for line.
    """
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_colour,
                                 width=line_width)
                       )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")

    return branch_line


def draw_clade(clade, x_coords: dict, y_coords: dict, line_shapes: List[dict], x_start: str = "0",
               line_colour: str = 'rgb(15,15,15)', line_width: int = 1):
    """
    Recursively define the tree lines (branches), starting from the argument clade.
    :param y_coords: dict of x-coordinates.
    :param x_coords: dict of y-coordinates.
    :param clade: name of the clade to start from, most likely tree root.
    :param x_start: x coord to start on, default 0
    :param line_shapes: list of line_shapes already created.
    :param line_colour: default line colour is black.
    :param line_width: width of lines.
    :return:
    """
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

    # Draw a horizontal line
    branch_line = get_clade_lines(orientation='horizontal',
                                  y_curr=y_curr,
                                  x_start=x_start,
                                  x_curr=x_curr,
                                  line_colour=line_colour,
                                  line_width=line_width)

    line_shapes.append(branch_line)

    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]

        line_shapes.append(get_clade_lines(orientation='vertical',
                                           x_curr=x_curr,
                                           y_bot=y_bot,
                                           y_top=y_top,
                                           line_colour=line_colour,
                                           line_width=line_width))

        # Draw descendants
        for child in clade:
            draw_clade(child, x_coords, y_coords, line_shapes, line_colour=line_colour, x_start=x_curr)


def get_nearest_neighbours(snp_distances: pd.DataFrame, node: str, *, do_join: bool = True) -> str | list:
    """
    Get a dictionary of the nearest neighbours and their snp distance to given node (sample or ID in df).
    :param snp_distances: Pandas dataframe, snp distance matrix with empty top left corner.
    :param node: the sample or ID to get the nearest neighbours for.
    :param do_join: if True, join nearest neighbours list into HTML string, if False return the raw list (default: True)
    :return: list of dicts, nearest neighbours and their distance
    """
    try:
        snp_distances = snp_distances.drop(node)  # drop row that contains itself, as minimum will always be itself
    except KeyError:
        # If the ID is in the metadata but is not in the snp_dist matrix:
        return ""
    nearest = list(
        snp_distances[
            snp_distances[node] == snp_distances[node].min()
            ]
        .index
    )
    nearest_neighbours = []
    for neighbour in nearest:
        try:
            neighbour_distance = int(snp_distances[node].loc[neighbour])
        except TypeError:
            print(f"Issues with snp distance matrix - check sample {node} only occurs once in the matrix: "
                  f"{snp_distances[node].loc[neighbour]}")
            sys.exit()
        nearest_neighbours.append(f'{neighbour}={neighbour_distance}')  # format nicely to give distance in SNPs

    if do_join:
        return "<br>".join(nearest_neighbours)  # for html
    else:
        return nearest_neighbours


def make_hover_text(metadata_df: pd.DataFrame,
                    id_column: str,
                    node_list: list) -> list:
    """
    Make the hover text for the plotly interactive plot. This will create a pop-up box on hover over a node.
    :param metadata_df: the Pandas dataframe of metadata, must have default numerical index names.
    :param id_column: name of column in metadata that corresponds to taxa in tree (str).
    :param node_list: list of nodes in order of appearance in the tree.
    :return: list of hover text.
    """
    # For each row in the metadata field with the index, get the index in text list
    hover_text = node_list.copy()
    for df_index, sample_id in enumerate(metadata_df[id_column]):
        for index, node_name in enumerate(hover_text):
            if node_name is not None and node_name.startswith(sample_id):
                node_index = index
                # Only keep entries in the metadata that match in the tree
                metadata_list = [
                    f'{column_name}: {metadata_df.loc[df_index, column_name]}<br>' for column_name in
                    metadata_df.columns.values
                ]
                hover_text[node_index] = hover_text[node_index] + "<br>" + str("".join(metadata_list))
    return hover_text


def get_colourings(metadata_df: pd.DataFrame,
                   id_column: str,
                   category: str,
                   number_of_nodes: int,
                   node_list: list,
                   intermediate_node_colour: str = 'rgb(100,100,100)') -> list:
    """
    For a given column in the metadata dataframe provided, create a list of colours (strings) such that each unique
    item in that column (category) is given its own colour. Categories (columns) with more than 48 unique items are not
    given colours and are not made available to colour the nodes (reasoning: >48 colours is difficult to visualise).

    :param metadata_df: the Pandas dataframe of metadata.
    :param id_column: the id column, automatically detected in the metadata dataframe.
    :param category: the column name to be coloured from the Pandas dataframe.
    :param number_of_nodes: number of nodes
    :param node_list: list of nodes
    :param intermediate_node_colour: colour for intermediate nodes (default = grey).
    :return: list of colours for each node in the order that the nodes occur.
    """
    # Create list of colours - use plotly built in, for total of 48 unique colours.
    available_colours = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    # Create a dict and assign each value in the category a different colour (sort list to keep colours consistent).
    colouring_dict = dict(zip(sorted(list(set(metadata_df[category]))), available_colours))
    # Create the list of colours - the list is as long as the number of nodes in the tree.
    list_of_colours = [intermediate_node_colour] * int(number_of_nodes)
    # For each row in the metadata field with the index, get the index in text list
    for df_index, sample_id in enumerate(metadata_df[id_column]):
        for index, node_name in enumerate(node_list):
            if node_name is not None and node_name.startswith(sample_id):
                # Only colour entries in the metadata that match in the tree
                list_of_colours[index] = (colouring_dict[metadata_df.loc[df_index, category]])
    return list_of_colours


def get_continuous_colourings(metadata_df: pd.DataFrame,
                              id_column: str,
                              date_category: str,
                              number_of_nodes: int,
                              node_list: list,
                              intermediate_node_colour: str = 'rgb(100,100,100)') -> list:
    """
    For a date column, create a list of colours (in rgb) such that dates are coloured in a gradient. Newest date is
    coloured in royal blue, through to the oldest date in maroon red.

    :param metadata_df: the Pandas dataframe of metadata.
    :param id_column: the id column, automatically detected in the metadata dataframe.
    :param date_category: the name of the column that contained dates to be coloured in a gradient from the metadata.
    :param number_of_nodes: number of nodes
    :param node_list: list of nodes
    :param intermediate_node_colour: colour for intermediate nodes (default = grey).
    :return: list of colours for each node in the order that the nodes occur.
    """
    # Copy metadata:
    metadata_copy = metadata_df.copy()
    # Make date column actually dates, set anything not a date to not a time, drop later.
    metadata_copy[date_category] = pd.to_datetime(metadata_copy[date_category],
                                                  errors='coerce',
                                                  format='mixed',
                                                  yearfirst=True,
                                                  dayfirst=True)
    # If all dates are empty, assign all nodes to black (rgb(0, 0, 0):
    if len(metadata_copy[date_category].unique()) == 1 \
            and any([pd.isna(i) for i in metadata_copy[date_category].unique()]):
        metadata_copy['date_delta'] = metadata_copy[date_category].apply(lambda x: '0')
        gradient_colouring_dict = dict([('0', 'rgb(0, 0, 0)')])
    # Dates will be given a colour ranging from maroon (most recent) to royal blue (oldest).
    # If there is only one date all samples will be maroon (rgb(0, 0, 151)). Any empty dates get black (rgb(0, 0, 0)).
    else:
        # Get most recent
        most_recent_date = metadata_copy[date_category].max()
        # Add date delta column
        metadata_copy['date_delta'] = metadata_copy[date_category].apply(lambda x: most_recent_date - x)
        # Convert to days (ints)
        metadata_copy['date_delta'] = metadata_copy['date_delta'].dt.days
        # Get max date delta to scale all date deltas between 0 and 1:
        max_delta = metadata_copy['date_delta'].max()
        max_delta = max_delta if max_delta != 0 else 1  # Date delta can be 0 if dates are the same empty
        metadata_copy['date_delta'] = metadata_copy['date_delta'].apply(lambda x: x / max_delta)
        # Make dict of colours based on the date delta:
        date_deltas = sorted(metadata_copy['date_delta'].dropna().to_list())
        colour_gradient = plotly.colors.sample_colorscale('Jet', date_deltas)
        gradient_colouring_dict = dict(zip(date_deltas, colour_gradient))
        # Replace NaN with 'no_date' and add to dict to return black:
        metadata_copy['date_delta'] = metadata_copy['date_delta'].fillna('no_date')
        gradient_colouring_dict['no_date'] = 'rgb(0, 0, 0)'
    # Create the list of colours - the list is as long as the number of nodes in the tree.
    list_of_gradient_colours = [intermediate_node_colour] * int(number_of_nodes)
    # For each row in the metadata field with the index, get the index in text list
    for df_index, sample_id in enumerate(metadata_copy[id_column]):
        for index, node_name in enumerate(node_list):
            if node_name is not None and node_name.startswith(sample_id):
                # Only colour entries in the metadata that match in the tree
                list_of_gradient_colours[index] = (gradient_colouring_dict[metadata_copy.loc[df_index, 'date_delta']])
    return list_of_gradient_colours


def inline_html_images(html_res_path: os.PathLike | str, input_html: str) -> str:
    input_html_path = Path(html_res_path, input_html)
    html_str = input_html_path.read_text()
    soup = bs(html_str, features="lxml")
    image_tags = soup.findAll("img")

    for tag in image_tags:
        # get the relative location of the image file,
        # read as bytes and encode as base64
        b64_img = base64.b64encode(Path(html_res_path).joinpath(Path(tag.get("src"))).read_bytes()).decode("utf-8")
        # change the src attribute to the b64 encoded str
        tag["src"] = f"data:image/png;base64,{b64_img}"

    return str(soup)


def generate_html(html_res_path: os.PathLike | str, input_fig: go.Figure) -> str:
    """
    Takes in a plotly figure, reads some resources from the html_res directory and generates a static HTML string as
    output.
    :param html_res_path: path to the html_res directory.
    :param input_fig: plotly go figure object.
    :return: string of html.
    """
    return_html = f"""<!DOCTYPE html>
<html lang="en">
    <head>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta charset="UTF-8"/>
        <title>Intreeactive report (generated: {datetime.datetime.now()})</title>
        <link rel="icon" href="data:image/png;base64,{base64.encodebytes(Path(html_res_path, 'favicon.png').read_bytes()).decode('utf-8')}" />
        <style>
            {Path(html_res_path, 'main.css').read_text()}
        </style>
        <script>
            {Path(html_res_path, 'tempmetadata.js').read_text()}
        </script>
    </head>
    <body>
        <div id="leftPanel">
                {
    input_fig.to_html(
        full_html=False,
        default_height=900,
    )
    }
        </div>
        <div id="rightPanel">
            <div id="ashleyDiv" class="fancyDiv">
                <table>
                    <tr>
                        <td>Highlight by ID</td>
                        <td>
                            <input type="text" id="highlightInput" oninput="doHighlight(this.value)"></input>
                            <button onclick="restoreOriginalColours(); highlightInputElm.value = ''">Clear</button>
                        </td>
                    </tr>
                    <tr>
                        <td>Label by</td>
                        <td id="labelDropdownCell">
                        </td>
                    </tr>
                    <tr>
                        <td>Labels</td>
                        <td>
                            <button onclick="showLabels()">Show labels</button>
                            <button onclick="hideLabels()">Hide labels</button>
                        </td>
                    </tr>
                    <tr>
                        <td>SNP threshold</td>
                        <td>
                            <input checked type="radio" name="snpThresholdRadio" value="Use min">
                            <label for="Use min">Use min</label><br>
                            <input type="radio" name="snpThresholdRadio" value="Set manually">
                            <label for="Set manually">Set manually</label>
                            <input type="number" value=1 min=1 id="snpThresholdSpinner" onchange="snpThresholdSpinnerChange()">
                        </td>
                    </tr>
                    <tr>
                        <td>Scale</td>
                        <td>
                            <button class="zoomBtn zoomBtnX" onclick="doZoomX(0)">X+</button>
                            <button class="zoomBtn zoomBtnX" onclick="doZoomX(1)">X-</button>
                            <button class="zoomBtn zoomBtnY" onclick="doZoomY(0)">Y+</button>
                            <button class="zoomBtn zoomBtnY" onclick="doZoomY(1)">Y-</button>
                            <button class="zoomBtn zoomBtnReset" onclick="resetScale()">Reset</button>
                        </td>
                    </tr>
                </table>
            </div>

            <div id="metadataDivContainer" class="fancyDiv">
                <h2>Search and compare metadata</h2>
                <div id="metadataDivControls">
                </div>
                <div id="metadataDiv">
                </div>
            </div>
        </div>
        <div id="helpModal" class="hidden">
            {inline_html_images(html_res_path, 'help.html')}
        </div>
    </body>
    <script>
        {Path(html_res_path, 'main.js').read_text()}
    </script>
</html>"""

    return return_html


def write_interactive_tree(*,
                           tree,
                           output_name: str | os.PathLike,
                           metadata: pd.DataFrame,
                           id_column: str = 'ID',
                           snp_distance_matrix: pd.DataFrame,
                           title: str = None) -> None:
    """
    Create an interactive phylogeny (html) file for a given phylogeny file.
    :param tree: Bio Phylo Tree object.
    :param output_name: A file name for or path and file name for the output (str). Do not add file format suffix.
    :param metadata: Pandas dataframe of metadata. Each column will be read in as hover text for the tree
        (pandas df)
    :param id_column: the column in the metadata that corresponds to the taxa in the tree. (str, default = 'ID')
    :param snp_distance_matrix: Pandas dataframe with all against all SNP distances. Column order must match row
        order (pandas df)
    :param title: string, title of the plot.
    :return: None, but html files are created
    """
    ##################
    # Set up:
    print(f'Creating Tree \'{title}\'... \n Number of leaves: {len(list(tree.get_terminals()))}')

    # Add nearest neighbours to the metadata dataframe
    metadata["Nearest_neighbour"] = metadata[id_column].apply(
        lambda x: get_nearest_neighbours(snp_distance_matrix, x))

    ###########
    # 1. Get x and y coordinated for the tree
    tree_x_coords = get_x_coordinates(tree)
    tree_y_coords = get_y_coordinates(tree)

    ###########
    # 2. Draw tree line shapes by iteratively moving through the branches starting at the tree root.
    # Set the line shapes for each node and leaf in the tree, ready for Plotly to plot:
    tree_line_shapes = []

    draw_clade(clade=tree.root,
               x_coords=tree_x_coords,
               y_coords=tree_y_coords,
               line_shapes=tree_line_shapes,
               line_colour='rgb(25,25,25)',
               line_width=1)

    # Combine all the shapes into one SVG path
    svg_path = []
    for record_dict in tree_line_shapes:
        svg_path.append(
            f"M {record_dict['x0']} {record_dict['y0']} L {record_dict['x1']} {record_dict['y1']}")
    svg_path = " ".join(svg_path)

    ###########
    # 3. Create the text for the hover text
    # Get the node coordinates:
    my_tree_clades = tree_x_coords.keys()
    x_nodes = []  # list of nodes x-coordinates
    y_nodes = []  # list of nodes y-coordinates
    node_list = []  # list of nodes as they appear in order in the tree

    for clade in my_tree_clades:
        x_nodes.append(tree_x_coords[clade])
        y_nodes.append(tree_y_coords[clade])
        node_list.append(clade.name)

    hover_text = make_hover_text(metadata, id_column, node_list)

    ###########
    # 4. Set colours for the nodes - select a suitable column from the metadata:
    # get the default category to colour on first:
    count = 1
    default_category = ""
    while count <= len(metadata.columns.values):
        if len(set(metadata[metadata.columns.values[count]])) <= 48:
            default_category = metadata.columns.values[count]
            break
        else:
            count += 1

    colourings = get_colourings(metadata_df=metadata,
                                id_column=id_column,
                                category=default_category,
                                number_of_nodes=len(x_nodes),
                                node_list=node_list)

    ###########
    # 5. Create traces for plotly plot - These are the nodes.
    trace = go.Scattergl(x=x_nodes,
                         y=y_nodes,
                         mode='markers',
                         marker=dict(color=colourings,
                                     size=10),
                         opacity=1.0,
                         text=hover_text,
                         hoverinfo='text')

    ###########
    # 6. Create the drop-down functionality
    # Add a fancy drop_down list to change the node colours in the interactive tree:
    drop_down_update = list([dict(
        buttons=[],
        direction='down',
        pad={'r': 10, 't': 10},
        showactive=True,
        x=1.01,
        xanchor='right',
        y=1.075,
        yanchor='top'
    )
    ])

    # Create a colour list for every column in the metadata dataframe if >=48 things to colour:
    for category_to_colour in metadata.columns.values:
        if "date" in category_to_colour.lower():
            drop_down_update_dict = dict(label=category_to_colour,
                                         method='update',
                                         args=[{'marker.color': [
                                             get_continuous_colourings(
                                                 metadata_df=metadata,
                                                 id_column=id_column,
                                                 date_category=category_to_colour,
                                                 number_of_nodes=len(x_nodes),
                                                 node_list=node_list
                                             )]}]
                                         )
            drop_down_update[0]['buttons'].append(drop_down_update_dict)
        else:
            list_of_categories = sorted(list(set(metadata[category_to_colour])))
            if len(list_of_categories) <= 48:
                drop_down_update_dict = dict(label=category_to_colour,
                                             method='update',
                                             args=[{'marker.color': [
                                                 get_colourings(
                                                     metadata_df=metadata,
                                                     id_column=id_column,
                                                     category=category_to_colour,
                                                     number_of_nodes=len(x_nodes),
                                                     node_list=node_list
                                                 )]}]
                                             )
                drop_down_update[0]['buttons'].append(drop_down_update_dict)

    ###########
    # 7. Add layout to plotly plot - add the lines between nodes.
    # Prep a title:
    if title:
        graph_title = f'{title}; (n={len(list(tree.get_terminals()))})'
    else:
        graph_title = f"Interactive Tree - {datetime.datetime.now()}"

    # The branches are already defined and stored as Plotly shapes that are included in the plot layout below:
    layout = go.Layout(title=dict(text=graph_title, yanchor='top', y=0.95),
                       font=dict(family='Arial', size=14),
                       showlegend=False,
                       autosize=True,
                       xaxis=dict(showline=True,
                                  zeroline=False,
                                  showgrid=False,
                                  ticklen=4,
                                  showticklabels=True,
                                  title='Branch Length'),
                       yaxis=dict(visible=False),
                       hovermode='closest',
                       plot_bgcolor='rgb(250,250,250)',
                       margin=dict(l=10, t=150),
                       # shapes=tree_line_shapes,  # lines for tree branches
                       updatemenus=drop_down_update  # This adds the drop-down menu to change the node colours.
                       )

    fig = go.Figure(data=[trace], layout=layout)

    # Add in the SVG path
    fig.update_layout(
        shapes=[
            dict(
                type="path",
                path=svg_path,
                line_color="rgb(25,25,25)",
                layer="below"
            )
        ]
    )

    # Add text annotations as a separate trace
    fig.add_trace(
        go.Scattergl(
            x=x_nodes,
            y=y_nodes,
            mode="text",
            text=[f"\t\t\t{x.split('<br>')[0]}" if x else "" for x in hover_text],
            textposition="middle right",
            visible=False
        )
    )

    ###########
    # 8. Jsonify data for javascript shenanigans
    metadata.index = metadata[id_column]
    # split NN data back out into an array
    metadata["Nearest_neighbour"] = metadata["Nearest_neighbour"].apply(lambda x: x.split("<br>"))

    # basket case to make sure we have normal jsonable types instead of weird numpy types
    matrix_json = json.loads(snp_distance_matrix.to_json(orient="split"))
    del matrix_json["columns"]

    # Technically we only need ~half the data because snp_dist(A,B) == snp_dist(B,A) but we're keeping all of it for now
    matrix_json = json.dumps(matrix_json)

    ###########
    # 9. Write outputs
    # use the html_res context manager to open the correct temp path to get access to the html_res folder needed for
    # outputs in step 9:

    # Write out the datasets to a javascript file -
    with as_file(files('html_res')) as html_res_path:
        with open(html_res_path / "tempmetadata.js", "w", encoding="utf-8") as outfile:
            # write out the metadata
            outfile.write("const inputMetadata = ")
            outfile.write(metadata.to_json(
                indent=4,
                orient="index"
            ))
            outfile.write("\n\n")

            # write out the distance matrix
            outfile.write("const inputSnpMatrix = ")
            outfile.write(matrix_json)
            outfile.write("\n")

        # Write the html file
        with open(f'{output_name}', "w", encoding="utf-8") as outfile:
            outfile.write(generate_html(html_res_path, fig))
