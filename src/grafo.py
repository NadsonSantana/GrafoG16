"""Estruturas de grafo baseadas em lista de adjacência."""

import csv


class Grafo:
    """Grafo simples ponderado, não-direcionado (ou direcionado, se dirigido=True).

    A lista de adjacência é um dicionário {vértice: {vizinho: peso}}.
    """

    def __init__(self, dirigido=False):
        self.dirigido = dirigido
        self.adj = {}

    def adicionar_vertice(self, v):
        if v not in self.adj:
            self.adj[v] = {}

    def adicionar_aresta(self, u, v, peso=1.0):
        self.adicionar_vertice(u)
        self.adicionar_vertice(v)
        self.adj[u][v] = peso
        if not self.dirigido:
            self.adj[v][u] = peso

    def remover_vertices(self, vertices):
        """Retorna uma cópia do grafo sem os vértices informados."""
        removidos = set(vertices)
        h = Grafo(self.dirigido)
        for u, viz in self.adj.items():
            if u in removidos:
                continue
            h.adj[u] = {v: p for v, p in viz.items() if v not in removidos}
        return h

    def subgrafo(self, vertices):
        manter = set(vertices)
        h = Grafo(self.dirigido)
        for u in manter:
            h.adj[u] = {v: p for v, p in self.adj[u].items() if v in manter}
        return h

    def vertices(self):
        return list(self.adj)

    def vizinhos(self, v):
        return self.adj[v]

    def grau(self, v):
        return len(self.adj[v])

    def num_vertices(self):
        return len(self.adj)

    def num_arestas(self):
        total = sum(len(viz) for viz in self.adj.values())
        return total if self.dirigido else total // 2

    def arestas(self):
        """Lista de arestas (u, v, peso); no caso não-direcionado, cada aresta aparece uma vez."""
        vistas = []
        for u, viz in self.adj.items():
            for v, p in viz.items():
                if self.dirigido or u < v:
                    vistas.append((u, v, p))
        return vistas

    @classmethod
    def de_csv(cls, caminho_arestas, caminho_vertices=None, dirigido=False, coluna_peso="distancia"):
        g = cls(dirigido)
        if caminho_vertices:
            with open(caminho_vertices, newline="") as f:
                for linha in csv.DictReader(f):
                    g.adicionar_vertice(linha["id"])
        with open(caminho_arestas, newline="") as f:
            for linha in csv.DictReader(f):
                g.adicionar_aresta(linha["origem"], linha["destino"], float(linha[coluna_peso]))
        return g
