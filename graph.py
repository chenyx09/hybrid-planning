from collections import defaultdict
import pdb

class Graph:
    def __init__(self):
        self.nodes = set()
        self.edges = defaultdict(list)
        self.distances = {}

    def add_node(self, value):
        self.nodes.add(value)

    def add_edge(self, from_node, to_node, distance):
        if to_node not in self.edges[from_node]:
            self.edges[from_node].append(to_node)
        if from_node not in self.edges[to_node]:
            self.edges[to_node].append(from_node)
        self.distances[(from_node, to_node)] = distance
        self.distances[(to_node, from_node)] = distance
    def remove_edge(self, from_node, to_node):
        if to_node in self.edges.keys():
            if from_node in self.edges[to_node]:
                self.edges[to_node].remove(from_node)
        if from_node in self.edges.keys():
            if to_node in self.edges[from_node]:
                self.edges[from_node].remove(to_node)
        if (from_node, to_node) in self.distances.keys():
            del self.distances[(from_node, to_node)]
        if (to_node, from_node) in self.distances.keys():
            del self.distances[(to_node, from_node)]
    def remove_node(self, node):
        try:
            if node in self.nodes:
                self.nodes.remove(node)
            if node in self.edges.keys():
                del self.edges[node]
        except:
            pdb.set_trace()
        for node1 in self.edges.keys():
            self.remove_edge(node, node1)
    def reachable_node(self,start):
        reachable_set = [start]
        done = False
        while not done:
            done = True
            for node in reachable_set:
                try:
                    for new_node in self.edges[node]:
                        if new_node not in reachable_set:
                            done = False
                            reachable_set.append(new_node)
                except:
                    pdb.set_trace()
        return reachable_set
    def copy(self):
        g = Graph()
        g.nodes = self.nodes.copy()
        g.edges = self.edges.copy()
        g.distances = self.distances.copy()
        return g

def graph_similar(g1,g2):
    return g1.nodes==g2.nodes and g1.edges==g2.edges


def dijsktra(graph, initial,dest):
    visited = {initial: 0}
    path = {}

    nodes = set(graph.nodes)

    while nodes:
        min_node = None
        for node in nodes:
            if node in visited:
                if min_node is None:
                    min_node = node
                elif visited[node] < visited[min_node]:
                    min_node = node

        if min_node is None:
            break


        nodes.remove(min_node)
        current_weight = visited[min_node]

        for edge in graph.edges[min_node]:
            try:
                weight = current_weight + graph.distances[(min_node, edge)]
            except:
                pdb.set_trace()
            if edge not in visited or weight < visited[edge]:
                visited[edge] = weight

    if dest in visited.keys():
        path = [dest]

        while path[0]!=initial:

            for node in graph.edges[path[0]]:
                if node in visited.keys() and node!=path[0]:
                    if abs(visited[node] - visited[path[0]]+graph.distances[path[0],node])<1e-5:

                        path = [node]+path
                        break
    else:
        path = []


    return visited, path
