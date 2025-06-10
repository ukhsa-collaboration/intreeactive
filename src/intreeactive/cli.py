import argparse
import os
import glob
import sys
import datetime

from intreeactive import intreeactive


#####################
# Setup CL arguments
def get_args():
    parser = argparse.ArgumentParser(
        prog='intreeactive',
        epilog='''
        Create an interactive tree. 
        It has interactive features such as move, zoom and pan, on-hover sample metadata, node colouring, custom node 
        label creation and identifying neighbours using custom SNP distance thresholds.
        All information and code needed is included in one html file. 
        '''
    )
    parser.add_argument(
        '--tree',
        '-t',
        dest='tree_path',
        type=str,
        required=True,
        help='Required: supply the path to the tree file. If file is not Newick format, specify the tree type using '
             '--tree-format/-tf'
    )
    parser.add_argument(
        '--metadata',
        '-m',
        dest='metadata_path',
        type=str,
        required=True,
        help='Required: supply the path to the metadata file. The first column will be used as the sample ID unless '
             'specified with --id-column/-id. This sample ID will be used to match samples in the tree, metadata and '
             'SNP distance matrix. If any samples are not present in the metadata but are in the tree or snp distance '
             'matrix, use the --ignore-ids/-x to ignore these.'
    )
    parser.add_argument(
        '--snp-distance-matrix',
        '-s',
        dest='snp_distance_matrix_path',
        type=str,
        required=True,
        help='Required: Supply path to the SNP distance matrix. Can use any seperator, the order of the columns must be '
             'identical to the order of the rows.'
    )
    parser.add_argument(
        '--tree-format',
        '-T',
        dest='tree_format',
        type=str,
        required=False,
        default='newick',
        help='Optional: if the tree file is not Newick (.new or .newick), supply the tree file format. '
             'Default: "Newick"'
    )
    parser.add_argument(
        '--outgroup',
        '-O',
        dest='outgroup',
        type=str,
        required=False,
        default=None,
        help='Optional: supply the name of the ID for the outgroup. If the outgroup is supplied, the tree will be '
             'rooted here.'
    )
    parser.add_argument(
        '--id-column',
        '-I',
        dest='id_column',
        type=str,
        required=False,
        default='ID',
        help='Optional: supply the name of the column that contains the ID to match samples in the metadata to the tree'
             'leaves and the SNP distance matrix. Default="ID"'
    )
    parser.add_argument(
        '--ignore',
        '-x',
        dest='ignore_ids',
        type=str,
        action='append',
        required=False,
        default=None,
        help='Optional: supply the name(s) of the ID(s) to be ignored - these are IDs that are present in the tree but '
             'not in the SNP distance matrix or the metadata, for example the outgroup/reference. Default=None.'
    )
    parser.add_argument(
        '--output',
        '-o',
        dest='output_file',
        type=str,
        required=False,
        default='interactive_tree',
        help='Optional: Filename or path with filename to be used as the output. Do not include the suffix, .html will '
             'be added. Intreeactive will not overwrite files with the same name. Use --force to overwrite a file with '
             'the provided file name of path. Default="interactive_tree".'
    )
    parser.add_argument(
        '--output-dir',
        '-d',
        dest='output_dir',
        type=str,
        required=False,
        default=None,
        help='Optional: Name of directory or path with directory to be used to save the output into. '
             'If it does not already exist, it will be created. Default=current working directory.'
    )
    parser.add_argument(
        '--title',
        '-y',
        dest='title',
        type=str,
        required=False,
        default=None,
        help='Optional: Title to be added to the interactive tree. Default="Interactive Tree - <date today>".'
    )
    parser.add_argument(
        '--force',
        dest='force',
        help='Overwrite the output directory if it already exists',
        action='store_true')
    return parser.parse_args()


def _handle_outdir(outdir: str | os.PathLike = None) -> os.PathLike:
    """
    Handle the output directory, using either the default or commandline arg. Create dir
    if needed.
    :param outdir: the outdir argument from commandline.
    :return: path to main output directory.
    """
    if outdir is not None:
        output_dir_path = os.path.abspath(outdir)
        if not os.path.exists(output_dir_path):
            print(f'Output directory {outdir} does not exist, creating...')
            os.mkdir(output_dir_path)
    else:
        output_dir_path = os.getcwd()
    return output_dir_path


def _check_output_exists(output_file_path: os.PathLike | str, force: bool) -> None:
    """
    Check if output file already exists, if so, exit unless --force is used.
    :param force: bool, arg for force
    :param output_file_path: Name of the output file. Default is 'current_working_directory/interactive_tree.html'.
    :return: None
    """
    if glob.glob(str(output_file_path)):
        if not force:
            message = f'File "{output_file_path}" exists. Use --force to force overwrite. Exiting.'
            print(message, file=sys.stderr)
            sys.exit()
        else:
            print(f'Files in directory {output_file_path} found, --force passed, overwriting...')


def setup_and_check_files(args: argparse.Namespace) -> os.PathLike | str:
    """
    Check if output directory exists, if not make it.
    :param args: arg parse object
    :return: Output path with filename.
    """
    # Set up main output dir:
    outdir_path: os.PathLike | str = _handle_outdir(args.output_dir)
    # Check output already exists:
    output_file_path = os.path.join(outdir_path, (args.output_file + '.html'))
    _check_output_exists(output_file_path, args.force)
    return output_file_path


#####################
def main():
    args = get_args()
    output_path = setup_and_check_files(args)

    ### Set up files:
    # Tree: read tree file (specify format) and parse as Bio.Phylo tree object, specifying an outgroup will root the tree.
    tree = intreeactive.read_in_tree(path_to_tree=args.tree_path,
                                     tree_format=args.tree_format,
                                     outgroup=args.outgroup
                                     )

    # Metadata: the metadata is used to add information to the hover text, as well as matching up the nearest
    # neighbours. This reads in from csv into pandas df, sets content to strings, gets the id_column if not specified
    # (This is the main ID of the sample). Dates are read in using pandas to_datetime, and prefers 2024-01-01, but
    # it will try to parse other formats, preferring year first, and then day first. Any column names that contain the
    # string "date" (not case-sensitive) will be given a colour gradient of dates.
    metadata_df, id_column = intreeactive.read_in_metadata(args.metadata_path, id_column=args.id_column)

    # SNP Distance Matrix: This is used to add functionality to the tree and find nearest neighbours.
    # This reads a file with any delimiter, sets the sample names (in the first column) to the index, and assumes
    # the order of the column names must match the row names (to store the index once).
    snp_distance_matrix = intreeactive.read_in_snp_dist_matrix(args.snp_distance_matrix_path)

    # This checks if sample IDs match up in the various files - the samples in the tree must have metadata
    # and nearest neighbour information.
    checked_metadata = intreeactive.check_ids(tree=tree,
                                              metadata=metadata_df,
                                              id_column=id_column,
                                              snp_dists=snp_distance_matrix,
                                              ignore_ids=args.ignore_ids)

    # The check_ids function also reduces the metadata down to only entries needed for the tree to save on computation:
    metadata_df = checked_metadata if not checked_metadata.empty else metadata_df

    # Set up the title for the plot:
    today = datetime.date.today().strftime("%Y%m%d")
    title = args.title if args.title else f"Interactive Phylogeny, {today}"

    # This is the main function - it takes in the Phylo tree object, a path or name to the output file, metadata Pandas
    # Dataframe, the name of the ID column
    intreeactive.write_interactive_tree(tree=tree,
                                        output_name=output_path,
                                        metadata=metadata_df,
                                        id_column=id_column,
                                        snp_distance_matrix=snp_distance_matrix,
                                        title=title)
