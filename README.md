# sage_classes
New class for graphs in sage which inherits the standard Graph-class but has additional features, such as local complementations and Pauli-measurements.

To use the class, start sage and load the script by writing `load("MyClasses.sage")`.

Example:
```
G=MyGraph(range(10))            #Creates a complete graph on the vertices 0,1,...,9
G*0                             #Does a local complementation (\tau) on vertex 0.
G-0                             #Deletes vertex 0
G.meas_seq([0,3],['X','Y'])     #Measures node 0 in X-basis and then node 3 in Y-basis.
G.is_LC_eq(graphs.StarGraph(9)) #Tests if G is LC-equivalent to the star graph.
G.can_meas(graphs.StarGraph(3)) #Tests if G can be measured into a star graph on 4 vertices.
```

To easily input larger graphs the script input_graph_single.py can be run. This is nothing super-fancy but can be convenient.
Run `python input_graph_single.py` (requires matplotlib and numpy) and input your graph by first placing vertices (right-click when done). Then draw edges by clicking on vertices (paths and walks can be input in one go) and finish a path/walk by right-clicking. To finish the full graph left-click one time outside the vertices and then right-click. If you're happy with your graph, click the green circle and the program outputs the adjacency-dictionary which can be copy-pasted into sage. The red circle restarts the process.
