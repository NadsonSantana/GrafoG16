"""Funções auxiliares dos testes. O NetworkX é usado apenas como oráculo de referência."""

import random

import networkx as nx

from src.grafo import Grafo


def para_nx(g):
    h = nx.DiGraph() if g.dirigido else nx.Graph()
    h.add_nodes_from(g.adj)
    h.add_weighted_edges_from(g.arestas())
    return h


def grafo_aleatorio(n, p, semente, dirigido=False, desconexo=True):
    """Grafo aleatório ponderado; com desconexo=True, não há arestas entre as duas metades."""
    rng = random.Random(semente)
    g = Grafo(dirigido)
    for v in range(n):
        g.adicionar_vertice(v)
    for u in range(n):
        for v in range(n):
            if u != v and (dirigido or u < v) and rng.random() < p:
                if desconexo and (u < n // 2) != (v < n // 2):
                    continue
                g.adicionar_aresta(u, v, round(rng.uniform(0.1, 5.0), 3))
    return g


def caminho(*vertices):
    g = Grafo()
    for u, v in zip(vertices, vertices[1:]):
        g.adicionar_aresta(u, v)
    return g
