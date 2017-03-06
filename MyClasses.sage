load("LC.sage")
class MyGraph(Graph):
    def __init__(self,*args,**kwargs):
        """Can be initialized in the same way as Graph(). An additional option is to give a list, which by default will give a complete graph on the vertices provided in the list, if format='empty' is given the graph is instead a graph on the vertices in the list with no edges."""
        try:
            super(MyGraph,self).__init__(*args,**kwargs)
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
    def N(self,X):
        """Returns the symmetric difference of all the neighbor-sets of the vertices in X. If X is a single vertex, the neighbors are returned."""
        try:
            neigh=Set([])
            for v in X:
                neigh^^=Set(self.neighbors(v))
            return neigh
        except TypeError:
            return self.neighbors(X)
    def tau(self,v,inplace=False):
        """Performs a tau on vertex v. If inplace, this is done on the graph itself and nothing is returned, otherwise a new graph is returned."""
        if inplace:
            G=self
        else:
            G=MyGraph(self)
        Na=G.neighbors(v)
        for e in Subsets(Na,2):
            if G.has_edge(e):
                G.delete_edge(e)
            else:
                G.add_edge(e)
        G.name('')
        if not inplace:
            return G
    def tau_seq(self,m,inplace=False):
        """Performs a sequence of local complementation of the vertices in m. If inplace, this is done on the graph itself and nothing is returned, otherwise a new graph is returned."""
        if inplace:
            G=self
        else:
            G=MyGraph(self)
        for v in m:
            G.tau(int(v),inplace=True)
        G.name('')
        if not inplace:
            return G
    def meas_Z(self,a,inplace=False):
        """Measures the graph at vertex a in Z-basis. If inplace, this is done on the graph itself and nothing is returned, otherwise a new graph is returned."""
        if inplace:
            G=self
        else:
            G=MyGraph(self)
        G.delete_vertex(a)
        G.name('')
        if not inplace:
            return G
    def meas_Y(self,a,inplace=False):
        """Measures the graph at vertex a in Y-basis. If inplace, this is done on the graph itself and nothing is returned, otherwise a new graph is returned."""
        if inplace:
            G=self
        else:
            G=MyGraph(self)
        G.tau(a)
        G.delete_vertex(a)
        G.name('')
        if not inplace:
            return G
    def meas_X(self,a,b=None,inplace=False):
        """Measures the graph at vertex a in X-basis. If inplace, this is done on the graph itself and nothing is returned, otherwise a new graph is returned."""
        if inplace:
            G=self
        else:
            G=MyGraph(self)
        if len(G.neighbors(a))==0:
            G.delete_vertex(a)
        else:
            if b==None:
                bp=G.neighbors(a)[0]
            else:
                bp=b
            if not G.has_edge((a,bp)):
                raise ValueError("b in not adjacent to a")
            G.tau_seq([bp,a,bp],inplace=True)
            G.delete_vertex(a)
        G.name('')
        if not inplace:
            return G
    def meas_seq(self,inds,bases,ret_meas=False,inplace=False):
        """Meas the graph on vertices inds and in bases. For X-measuremnts the first neighbor is choosen as special if there is a neighbor. If inplace, this is done on the graph itself and nothing is returned, otherwise a new graph is returned."""
        if not len(inds)==len(bases):
            raise ValueError("List of vertices and list of bases should be of equal length.")
        if inplace:
            G=self
        else:
            G=MyGraph(self)
        if ret_meas:
            meas=""
        for (i,b) in zip(inds,bases):
            if b=="Y":
                G.meas_Y(i,inplace=True)
                if ret_meas:
                    meas+="Y"
            elif b=="Z":
                G.meas_Z(i,inplace=True)
                if ret_meas:
                    meas+="Z"
            elif b=="X":
                if ret_meas:
                    if len(G.neighbors(i))==0:
                        meas+="x"
                    else:
                        meas+=str(G.neighbors(i)[0])
                G.meas_X(i,inplace=True)
            else:
                raise ValueError("Not correct basis")
        if ret_meas:
            return (G,meas)
        return G
    def is_LC_eq(self,other,test_Vs=True,allow_disc=False):
        """Tests whether the graphs G1 and G2 are LC-eq, if allow_disc, then disconnected graphs are allowed and the components are tested seperately. If test_Vs=False, the vertex-sets are not tested to be equal."""
        if allow_disc:
            comp1=Set(map(Set,G1.connected_components()))
            comp2=Set(map(Set,G2.connected_components()))
            if comp1!=comp2:
                return False
            for s in comp1:
                if not MyGraph(self.subgraph(s)).is_LC_eq(other.subgraph(s)):
                    return False
            return True
        else:
            if not (self.is_connected() and other.is_connected()):
                raise ValueError("Graphs can't be disconnected if allow_disc=False")
            if test_Vs:
                if set(self.vertices())==set(other.vertices()):
                    return test_cond(self.adjacency_matrix(),other.adjacency_matrix())
                else:
                    return False
            return test_cond(self.adjacency_matrix(),other.adjacency_matrix())
    def can_meas(self,other,inds=None,start=True):
        """Test if G can be turned into Gp by measuring vertices in inds in some bases. Returns a (ans,meas) where ans is True of False depending on if it is possible to meas G to Gp or not and meas is a dictionary with which indices should be measured in which basis. Can take a long time for larger graphs, in particular if there is not meas-sequence."""
        if not (self.is_connected() and other.is_connected()):
            raise NotImplementedError("Both graphs should be connected")
        if start:
            if inds==None:
                inds=list(set(self.vertices())-set(other.vertices()))
            rem_inds=set(self.vertices())-set(inds)
            if not rem_inds==set(other.vertices()):
                print("V(Gp)!=V(G)\inds")
                return (False,None)
        if len(inds)==0:
            if self.is_LC_eq(other):
                return (True,{})
            return (False,None)
        noneq=[]
        for bas in ["Z","Y","X"]:
            (post,meas)=self.meas_seq(inds[:1],[bas],ret_meas=True)
            if post.is_connected():
                if not any(post.is_LC_eq(x) for x in noneq):
                    (ans,bases)=post.can_meas(other,inds=inds[1:],start=False)
                    if ans:
                        bases[inds[0]]=meas
                        return (True,bases)
                    noneq+=[post]
        return (False,None)
    def can_meas_all(self,other,inds=None):
        """Test if G can be turned into Gp by measuring vertices in inds in some bases, returns all possible sequences. Takes a long time for larger graphs."""
        if inds==None:
            inds=list(set(self.vertices())-set(other.vertices()))
        rem_inds=set(self.vertices())-set(inds)
        if not rem_inds==set(other.vertices()):
            print("V(Gp)!=V(G)\inds")
            return (False,None)
        rem_inds=list(rem_inds)
        S=[]
        for bases in Tuples("ZYX",len(inds)):
            (post,meas)=self.meas_seq(inds,bases[::-1],ret_meas=True)
            if post.is_connected():
                if post.is_LC_eq(other):
                    meas_dct={inds[i]:meas[i] for i in range(len(inds))}
                    S+=[meas_dct]
        return S
