import matplotlib.pyplot as plt
import numpy as np
import sys
def draw_sphere(c,r=0.02,n=30,col='b'):
    t=np.linspace(0,2*np.pi,num=n)
    plt.plot(r*np.cos(t)+c[0],r*np.sin(t)+c[1],col)
def draw_edge(a,b):
    plt.plot([V[a][0],V[b][0]],[V[a][1],V[b][1]],'r')
def whichvertex(p):
    for i in range(n):
        if np.linalg.norm(np.array(V[i])-np.array(p))<0.02:
            return i
    return -1
def draw_yn():
    draw_sphere(c1,r=0.03,col='g')
    draw_sphere(c2,r=0.03,col='r')
    plt.text(0.82,0.9,"Satisfied?")
def which_ans(p):
    for (i,c) in enumerate([c1,c2]):
        if np.linalg.norm(np.array(c)-np.array(p))<0.03:
            return i
    return -1
c1=(0.82,0.8)
c2=(0.91,0.8)

restart=True
while restart:
    plt.clf()
    plt.plot()
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.axes().set_aspect('equal','datalim')
    plt.draw()
    print("Please click")
    V=plt.ginput(n=0,timeout=0,mouse_stop=3,mouse_pop=2)
    n=len(V)

    plt.clf()
    plt.plot()
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.axes().set_aspect('equal','datalim')
    for p in V:
        draw_sphere(p)

    finished=False
    E=[]
    while not finished:
        plt.draw()
        pathcoord=plt.ginput(n=0,timeout=0,mouse_stop=3,mouse_pop=2)
        if len(pathcoord)>1:
            path=list(map(whichvertex,pathcoord))
            if all(map(lambda x:x!=-1,path)):
                for i in range(len(list(path))-1):
                    E+=[[path[i],path[i+1]]]
        else:
            break
        for e in E:
            draw_edge(e[0],e[1])
        plt.draw()
    A=np.zeros((n,n))
    for e in E:
        A[tuple(e)]=1
        A[tuple(e[::-1])]=1
    ind=np.arange(n)
    A={i:list(ind[A[i]==1]) for i in ind}
    draw_yn()
    plt.draw()
    no_ans=True
    while no_ans:
        ans_coord=plt.ginput(n=1,timeout=0)
        ans=which_ans(ans_coord)
        if ans==0:
            restart=False
            break
        elif ans==1:
            break
print(A)
