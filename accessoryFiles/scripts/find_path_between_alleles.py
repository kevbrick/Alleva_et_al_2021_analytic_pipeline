# program to check if there is exist a path between two vertices
# of a graph
import sys
import copy
import re
from collections import defaultdict

#This class represents a directed graph using adjacency list representation
class Graph:

  def __init__(self,allelesList):
    self.V= len(allelesList) #No. of vertices
    self.alleles = allelesList #List of alleles
    self.shortestPath = 100
    self.nshortestPath = 0
    self.graph = defaultdict(list) # default dictionary to store graph
    self.paths = list()

  # function to add an edge to graph
  def addEdge(self,u,v):
    self.graph[v].append(u)

  # function to add an edge to graph
  def printMe(self,u,v):
    print(self.graph[u])

  # Use BFS to check path between s and d
  def isReachable(self, s, d):
    # Mark all the vertices as not visited
    visited =[False]*(self.V)

    # Create a queue for BFS
    queue=[]

    # Mark the source node as visited and enqueue it
    queue.append(s)

    visited[s] = True
    steps = 0
    bestPath = 999999999999
    retpath=list()

    while queue:

      #Dequeue a vertex from queue
      n = queue.pop(0)

      # If this adjacent node is the destination node,
      # then return true
      if n == d:
        steps = steps + 1
        if len(queue) > 0 & len(queue) < bestPath:
          bestPath = len(queue)
          print("A path [" + str(bestPath) + "]--> ")
          retpath = queue
          print(retpath)

      #  Else, continue to do BFS
      for i in self.graph[n]:
        if visited[i] == False:
          queue.append(i)
          visited[i] = True

    # If BFS is complete without visited d
    print(retpath)
    return retpath

  def printAllPathsUtil(self, u, d, visited, path):

    # Mark the current node as visited and store in path
    visited[u]= True
    path.append(u)

    # If current vertex is same as destination, then print
    # current path[]

    #if (len(path) > min(self.shortestPath,int(sys.argv[1]))):
    if (len(path) <= int(sys.argv[1])):
      if u == d:
        mypath = self.alleles[path[0]]
        for i in path:
          if i != path[0]:
            mypath = mypath + "--" + self.alleles[i]

        if len(path) < self.shortestPath:
          self.shortestPath  = len(path)
          self.nshortestPath = 1
          allPaths.clear()

        if len(path) == self.shortestPath:
          allPaths.append(mypath)
          #print(mypath)
          #print("Shortest path: %i [%i]" % (self.shortestPath - 1,self.nshortestPath))
          self.nshortestPath = self.nshortestPath + 1

      else:
        # If current vertex is not destination
        # Recur for all the vertices adjacent to this vertex
        for i in self.graph[u]:
          if visited[i]== False:
            self.printAllPathsUtil(i, d, visited, path)

    # Remove current vertex from path[] and mark it as unvisited
    path.pop()
    visited[u]= False

  # Prints all paths from 's' to 'd'
  def printAllPaths(self, s, d):

    # Mark all the vertices as not visited
    visited =[False]*(self.V)

    # Create an array to store paths
    path = []

    # Call the recursive helper function to print all paths
    self.printAllPathsUtil(s, d, visited, path)

###################################################
###################################################
# Create a graph given in the above diagram
Counter = 0

alleles = list()
directParents = dict()

with open("dict.txt") as f:
  for line in f:
    (child, parent, N) = line.split()
    directParents[parent+"___"+child] = 1
    if child not in alleles:
      if child != 'PrZFA':
        alleles.append(child)

baseGraph = Graph(alleles)

with open("dict.txt") as f:
  for line in f:
    (child, parent, N) = line.split()
    if child != 'PrZFA':
      if len(sys.argv) > 2 and sys.argv[2] == 'pub':
        if not re.match('v:', child):
          baseGraph.addEdge(alleles.index(parent),alleles.index(child))
      else:
        baseGraph.addEdge(alleles.index(parent),alleles.index(child))

if len(sys.argv) > 3:
  parentAlleles = list(); parentAlleles.append(sys.argv[2])
  childAlleles = list(); childAlleles.append(sys.argv[3])
else:
  parentAlleles = alleles
  childAlleles = alleles

for parent in parentAlleles:
  for child in childAlleles:
    pcKey = parent+"___"+child

    if (parent == child) or (pcKey in directParents.keys()):
      if parent == child:
        print("%s\t%s\t%s\t%s\t" % (parent, child, "1", "1"))
      else:
        print("%s\t%s\t%s\t%s\t" % (parent, child, "1", "1"))
    else:
      g = copy.copy(baseGraph)
      allPaths=list()
      g.printAllPaths(alleles.index(child), alleles.index(parent))

      if len(allPaths) > 0:
        counter = allPaths[0].count('--')
        print("%s\t%s\t%s\t%s\t" % (parent, child, counter, len(allPaths)))
      else:
        print("%s\t%s\t%s\t%s\t" % (parent, child, "NA", "NA"))

#for p in allPaths:
#  print (p)
#
# lstUV = g.isReachable(alleles.index(u), alleles.index(v))
# lstVU = g.isReachable(alleles.index(v), alleles.index(u))
#
# if len(lstUV) > 0:
#
#   print(lstUV)
#   pPath = u
#   for step in lstUV:
#     pPath = pPath + "--" + alleles[step-1]
#   pPath = pPath + "--" + v
#
#   print("There is a path from %s to %s [path = %s ;steps = %i]" % (u,v,pPath, len(lstUV)))
# else :
#   print("There is no path from %s to %s" % (u,v))
#
# if len(lstVU) > 0:
#
#   pPath = u
#   for step in lstVU:
#     pPath = pPath + "--" + alleles[step-1]
#   pPath = pPath + "--" + v
#
#   print("There is a path from %s to %s [path = %s; steps = %i]" % (v,u,pPath, len(lstVU)))
# else :
#   print("There is no path from %s to %s" % (v,u))
