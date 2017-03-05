class MyGraph(Graph):
    def __init__(self,*args,**kwargs):
        """Can be initialized in the same way as Graph(). An additional option is to give a list, which by default will give a complete graph on the vertices provided in the list, if format='empty' is given the graph is instead a graph on the vertices in the list with no edges."""
        try:
            super(MyGraph,self).__init__(*args,**kwargs)
            try:
                self._order=list(args[0].get_order())
            except:
                self._order=self.vertices()
        except StandardError as e:
            try:
                data=list(args[0])
                n=len(data)
                if kwargs.has_key('format'):
                    if kwargs['format']=='empty':
                        A=zero_matrix(GF(2),n)
                    else:
                        A=ones_matrix(GF(2),n)+identity_matrix(GF(2),n)
                else:
                    A=ones_matrix(GF(2),n)+identity_matrix(GF(2),n)
                super(MyGraph,self).__init__(A)
                self.relabel(data)
                self._order=self.vertices()
            except:
                raise e
    def __add__(self,other):
        """Addition of graphs is given by addition of their corresponding adjacency matrices, or equivalently as the symmetric difference of their edge-sets."""
        Gself=MyGraph(self)
        Gother=MyGraph(other)
        Vself=Set(Gself.vertices())
        Vother=Set(Gother.vertices())
        Gself.add_vertices(Vother-Vself)
        Gother.add_vertices(Vself-Vother)
        vert_order=Gself.vertices()
        G=MyGraph(matrix(GF(2),Gself.adjacency_matrix(vertices=vert_order))+matrix(GF(2),Gother.adjacency_matrix(vertices=vert_order)))
        G.relabel(vert_order)
        return G
    def __sub__(self,other):
        """Subtracting something from a graph means deleting that vertex, if the vertex does not exist an error is raised."""
        G=MyGraph(self)
        G.delete_vertex(other)
        G._order.remove(other)
        return G
    def __mul__(self,m):
        """G*v gives the graph where a local complementation is done on vertex v. If v is a list or a string, a sequence of local complementations are performed, starting with the first in the list."""
        G=MyGraph(self)
        try:
            for v in m:
                G.tau(int(v),inplace=True)
        except TypeError:
            G.tau(m,inplace=True)
        return G
    def relabel(self,data):
        super(MyGraph,self).relabel(data)
        self._order=self.vertices()
    def get_order(self):
        return self._order
    def set_order(self,order):
        self.relabel(sorted(order))
        self._order=order
    def N(self,X):
        try:
            neigh=Set([])
            for v in X:
                neigh^^=self.neighbors(v)
            return neigh
        except TypeError:
            return self.neighbors(X)
    def tau(self,v,inplace=False):
        Na=self.neighbors(v)
        if inplace:
            for e in Subsets(Na,2):
                if self.has_edge(e):
                    self.delete_edge(e)
                else:
                    self.add_edge(e)
            self.name('')
        else:
            G=MyGraph(self)
            for e in Subsets(Na,2):
                if G.has_edge(e):
                    G.delete_edge(e)
                else:
                    G.add_edge(e)
            G.name('')
            return G
    def tau_seq(self,m,inplace=False):
        if inplace:
            for v in m:
                self.tau(int(v),inplace=True)
        else:
            G=MyGraph(self)
            for v in m:
                G.tau(v,inplace=True)
            return G
    def to_isotropic(self,Ain=None,Bin=None):
        return Isotropic(self,Ain=Ain,Bin=Bin)
