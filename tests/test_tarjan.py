"""Tarjan: componentes fortemente conexas e pontes/articulações."""

import networkx as nx

from src import algoritmos as alg
from src.grafo import Grafo
from tests.auxiliares import grafo_aleatorio, para_nx


def test_scc_caso_conhecido():
    g = Grafo(dirigido=True)
    for u, v in [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 4), (6, 6)]:
        g.adicionar_aresta(u, v)
    assert sorted(sorted(c) for c in alg.tarjan(g)) == [[1, 2, 3], [4, 5], [6]]


def test_scc_aleatorio():
    for semente in range(5):
        g = grafo_aleatorio(50, 0.05, semente, dirigido=True)
        esperado = sorted(sorted(c) for c in nx.strongly_connected_components(para_nx(g)))
        assert sorted(sorted(c) for c in alg.tarjan(g)) == esperado


def test_pontes_articulacoes():
    for semente in range(5):
        g = grafo_aleatorio(60, 0.06, semente)
        h = para_nx(g)
        pontes, art = alg.pontes_articulacoes(g)
        assert {frozenset(p) for p in pontes} == {frozenset(p) for p in nx.bridges(h)}
        assert set(art) == set(nx.articulation_points(h))
