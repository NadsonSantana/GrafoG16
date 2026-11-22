"""Prim e Kruskal (floresta geradora mínima)."""

import math

import networkx as nx

from src import algoritmos as alg
from tests.auxiliares import para_nx


def test_prim_kruskal(grafo):
    peso_nx = sum(d["weight"] for *_, d in nx.minimum_spanning_edges(para_nx(grafo), data=True))
    pr, kr = alg.prim(grafo), alg.kruskal(grafo)
    n_comp = len(alg.componentes_conexas(grafo))
    assert len(pr) == len(kr) == grafo.num_vertices() - n_comp
    assert math.isclose(sum(w for *_, w in pr), peso_nx)
    assert math.isclose(sum(w for *_, w in kr), peso_nx)


def test_union_find():
    uf = alg.UnionFind(range(5))
    assert uf.unir(0, 1) and uf.unir(2, 3) and not uf.unir(1, 0)
    assert uf.achar(0) == uf.achar(1) != uf.achar(2)


def test_floresta_do_grafo_real(real):
    assert len(alg.kruskal(real)) == len(alg.prim(real)) == 451 - 62
