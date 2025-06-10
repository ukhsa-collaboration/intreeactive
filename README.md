# 🌲 intreeactive 🌲

intreeactive is used to create a single simple interactive tree in HTML format. It is primarily for
datasets that can leverage SNP distance information. This was primarily developed for TB, but can be
used for any data where a distance matrix is used.

## 🌵 Features 🌵

Plotly in Python is used to plot the tree and take advantage of in-built
interactivity:

- Zoom
- Pan
- Hover text
- Choose node colour

The hover text shows all columns provided in a metadata table, and also the
nearest neighbour(s) with the SNP distance for a highlighted sample.
The node colours are chosen from a drop dox box, where each option is a metadata
column, if there are fewer than 49 unique items in that column (more than 49 is
too many to reasonable colour).

Custom Javascript and css has been added to provided a range of features such
as:

- Search for ID
- Choose leaf labels from metadata column
- Show or hide labels
- Compare metadata between two IDs
- highlight selected nearest neighbour from list.
- Set SNP threshold to show nearest neighbours (default is use minimum, user
  can set this.)

## 🌱 Info 🌱

The code is available through a Python script, which handles the custom
JavaScript and css to create one html file.

There is also the addition of a standalone html file that renders the SNP
distance matrix with interactivity.

Some key functions were used from https://github.com/empet/Phylogenetic-trees
This code utilises some functions from a notebook published in this git repo
(commit - 93e8e80).

## 🔧 Setup 🔧

This is a python project that utilises some prefabricated pngs, javascript and html files. It can be pip installed,
and you may wish to do this in a conda or virtual environment.

*NB: The instructions below use a conda environment. You may wish to use some other form of environment management -
such as venv*

1. Clone the repo using either ssh (or https), such as:  
   `git clone git@github.com:ukhsa-collaboration/intreeactive.git`

2. Move to root of repo:  
   `cd intreeactive`

3. Create and activate conda environment (has been tested with python=3.12.5):  
   `conda create --name intreeactive python=3.12.5`

   `conda activate intreeactive`

4. Install the distribution package:   
   `pip install .`

5. Check it has installed by running the below which should bring up the help text:  
   `intreeactive -h`

***Note: for some Windows installations, the executable might have .exe appended.***

## 🏁 Usage 🏁

### ⏩ Quick start ⏩

There is some example data included for a quick start guide. Simply run the command below:

```
intreeactive.exe --tree example/tb_in_middle_earth_tree.new --metadata example/tb_in_middle_earth_metadata.csv \
--snp-distance-matrix example/tb_in_middle_earth_snpdists_matrix.txt -d example -o example_interactive_tree \
-O AL123456.3 -x AL123456.3 -I 'ID_col' --title 'Example interactive phylogeny: TB in Middle Earth.' --force
```

This will run intreeactive with the example tree, metadata and SNP distance matrix. The output will be saved in the
`example` directory with the name `example_interactive_tree` (`.html` is automatically added to the end of the file
name.). The tree will be rooted at the reference 'AL123456.3', and as this reference does not appear in the metadata,
it is ignored. The ID column in the metadata (the first column) is called 'ID_col'. The title on the tree will
be `Example interactive phylogeny - TB in Middle Earth.`. The file already
exists as an example output, so force overwrite is applied.

### 📃 Inputs: 📃

| Argument                  | Required? | Description                                                                                                                                                                                                                                                                                                                                                |
|---------------------------|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --tree, -t                | Yes       | Supply the path to the tree file. If file is not Newick format, specify the tree type using '--tree-format/-T'.                                                                                                                                                                                                                                            |
| --metadata, -m            | Yes       | Supply the path to the metadata file. The first column will be used as the sample ID unless specified with --id-column/-I. This sample ID will be used to match samples in the tree, metadata and SNP distance matrix. If any samples are not present in the metadata but are in the tree or snp distance matrix, use the --ignore-ids/-x to ignore these. |
| --snp-distance-matrix, -s | Yes       | Supply path to the SNP distance matrix. Can use any seperator, the order of the columns must be identical to the order of the rows.                                                                                                                                                                                                                        |
| --tree-format, -T         | No        | Optional: if the tree file is not Newick (.new or .newick), supply the tree file format. Default: "Newick"                                                                                                                                                                                                                                                 |
| --outgroup, -O            | No        | Optional: supply the name of the ID for the outgroup. If the outgroup is supplied, the tree will be rooted here.'                                                                                                                                                                                                                                          |
| --id-column, -I           | No        | Optional: supply the name of the column that contains the ID to match samples in the metadata to the tree leaves and the SNP distance matrix. Default="ID"                                                                                                                                                                                                 |
| --ignore, -x              | No        | Optional: supply the name(s) of the ID(s) to be ignored - these are IDs that are present in the tree but not in the SNP distance matrix or the metadata, for example the outgroup/reference. Can supply one or many -x arguments. Default=None.                                                                                                            |
| --output, -o              | No        | Optional: Filename or path with filename to be used as the output. Do not include the suffix, .html will be added. Intreeactive will not overwrite files with the same name. Use --force to overwrite a file with the provided file name of path. Default="interactive_tree".                                                                              |
| --output-dir, -d          | No        | Optional: Name of directory or path with directory to be used to save the output into. If it does not already exist, it will be created. Default=current working directory.'                                                                                                                                                                               |
| --title, -y               | No        | Optional: Title to be added to the interactive tree. Default="Interactive Tree - <date today>".                                                                                                                                                                                                                                                            |
| --force, -f               | No        | Overwrite the output directory if it already exists, default = False.                                                                                                                                                                                                                                                                                      |

## 📤 Outputs: 📤

There is one output file called `<output>.html`. This can now be opened in a browser. It contains all the javascript,
html and css needed to render the tree, which uses Plotly libraries. It also has custom javascript to include the SNP
distance information and functionality to find neighbours in the panel on the right hand side.

There is an in-built help menu, accessible by clicking the question mark button on the bottom left of the screen.

## 🔨 Troubleshooting: 🔨

- There are pytests available in repository_root/tests. 
