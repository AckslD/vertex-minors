# Code for local complementations and vertex-minor

The file GraphClasses.sage contain classes for sage that can be used to study graphs under local complementations and vertex-deletions, i.e. vertex-minors.
Using the class SimpleGraph one can efficiently test whether two graphs are equivalent under local complementations and therefore if their corresponding graph states are equivalent under local Clifford operations.
Furthermore one can test if a graph is a vertex-minor of another using the method can_meas.
If the a graph is a vertex-minor of another then its corresponding graph state is a qubit-minor of the other, i.e. can be reached by local Clifford operations, local Pauli measurement and classical communication.
Furthermore there are classes for isotropic systems which describe equivalence classes of graphs under local complementations.
Also a class for double occurence words, related to circle graphs and their equivalence classes.
More details of the theory can be found in....

To use the classes, start sage and load the script by writing `load("GraphClasses.sage")`.

Example:
```
G=SimpleGraph(range(10))           #Creates a complete graph on the vertices 0,1,...,9
G*0                                #Does a local complementation (\tau) on vertex 0.
G-0                                #Deletes vertex 0
G.meas_seq([0,3],['X','Y'])        #Measures node 0 in X-basis and then node 3 in Y-basis.
G.is_LC_eq(graphs.StarGraph(9))    #Tests if G is LC-equivalent to the star graph.
G.can_meas(graphs.StarGraph(3))    #Tests if G can be measured into a star graph on 4 vertices.
```

To easily input larger graphs the script input_graph_single.py can be run. This is nothing super-fancy but can be convenient.
Run `python input_graph_single.py` (requires matplotlib and numpy) and input your graph by first placing vertices (right-click when done). Then draw edges by clicking on vertices (paths and walks can be input in one go) and finish a path/walk by right-clicking. To finish the full graph left-click one time outside the vertices and then right-click. If you're happy with your graph, click the green circle and the program outputs the adjacency-dictionary which can be copy-pasted into sage. The red circle restarts the process.

The files in MONA implements the monadic second-order formula describing whether a graph is a vertex-minor of another, which we describe in the paper ??.
