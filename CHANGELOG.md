# Changelog

# v1.2.0
Displaying internal nodes is now a command line option. 

## Added: 
- command line flag for exclude-internal to choose to exclude internal nodes in the tree. 

## Changed:
- node x and y coords, colours, hover_text are now stored in a dict, with the clade name (ID) as a 
key. The 'is_terminal' boolean is also stored in this dict. 
- the list of nodes, colours and hover text supplied to plotly to plot are filtered to remove those 
where is_terminal is false if --exclude-internal is used. 
- Default colour is always the ID, and the colours will cycle (will no be unique.)