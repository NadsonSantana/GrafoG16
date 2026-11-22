"""Dijkstra, Bellman-Ford e Floyd-Warshall (puro e numpy)."""

import math

import networkx as nx
import numpy as np

from src import algoritmos as alg
from src.grafo import Grafo
from tests.auxiliares import para_nx


def test_bellman_ford_peso_negativo():
    g = Grafo(dirigido=True)
    g.adicionar_aresta("a", "b", 4)
    g.adicionar_aresta("a", "c", 2)
    g.adicionar_aresta("c", "b", -3)
    dist, ciclo = alg.bellman_ford(g, "a")
    assert dist["b"] == -1 and not ciclo
    g.adicionar_aresta("b", "a", 0.5)  # ciclo a-c-b-a com peso -0.5
    assert alg.bellman_ford(g, "a")[1]


def test_dijkstra_bellman_floyd(grafo):
    h = para_nx(grafo)
    vs, matriz = alg.floyd_warshall(grafo)
    idx = {v: i for i, v in enumerate(vs)}
    for s in list(grafo.adj)[:10]:
        esperado = nx.single_source_dijkstra_path_length(h, s)
        dj = alg.dijkstra(grafo, s)
        bf, _ = alg.bellman_ford(grafo, s)
        assert dj.keys() == esperado.keys() == bf.keys()
        for v, d in esperado.items():
            assert math.isclose(dj[v], d) and math.isclose(bf[v], d)
            assert math.isclose(matriz[idx[s]][idx[v]], d)


def test_floyd_warshall_numpy_igual_ao_puro(grafo):
    vs, puro = alg.floyd_warshall(grafo)
    vs_np, vetorizado = alg.floyd_warshall_numpy(grafo)
    assert vs == vs_np
    assert np.allclose(np.array(puro), vetorizado)


def test_todas_as_fontes(grafo):
    dj, bf = alg.dijkstra_todos(grafo), alg.bellman_ford_todos(grafo)
    for s in dj:
        assert dj[s].keys() == bf[s].keys()
        assert all(math.isclose(dj[s][t], bf[s][t]) for t in dj[s])
