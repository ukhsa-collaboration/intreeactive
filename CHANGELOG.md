# Changelog

# v1.2.0
Displaying internal nodes is now a command line option. 
Setting the x-axis range is an arg to write_interactive_tree function in the lib.

## Added: 
- command line flag for exclude-internal to choose to exclude internal nodes in the tree. 
- x-axisrange is an arg to write_interactive_tree, default is None which results in plotly autoscaling the x-axis. Expects None or list of ints with length 2, e.g. [-10, 10]. Raises ValueError if does not match format.

## Changed:
- node x and y coords, colours, hover_text are now stored in a dict, with the clade name (ID) as a 
key. The 'is_terminal' boolean is also stored in this dict. 
- the list of nodes, colours and hover text supplied to plotly to plot are filtered to remove those 
where is_terminal is false if --exclude-internal is used. 
- Default colour is always the ID, and the colours will cycle (will no be unique.)