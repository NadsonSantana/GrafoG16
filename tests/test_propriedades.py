"""Parte III: centralidades, eficiência, geradores e lei de potência."""

import random

import networkx as nx
import numpy as np
import pytest

from src import propriedades as prop
from tests.auxiliares import grafo_aleatorio, para_nx


def test_betweenness(grafo):
    esperado = nx.betweenness_centrality(para_nx(grafo), normalized=False)
    obtido = prop.betweenness(grafo)
    for v in esperado:
        assert obtido[v] == pytest.approx(esperado[v])


def test_pagerank():
    for dirigido in (False, True):
        g = grafo_aleatorio(40, 0.1, 7, dirigido=dirigido)
        esperado = nx.pagerank(para_nx(g), weight=None, tol=1e-10, max_iter=1000)
        obtido = prop.pagerank(g)
        for v in esperado:
            assert obtido[v] == pytest.approx(esperado[v], abs=1e-6)


def test_eficiencia(grafo):
    assert prop.eficiencia_global(grafo) == pytest.approx(nx.global_efficiency(para_nx(grafo)))


def test_geradores():
    rng = random.Random(0)
    er = prop.gerar_erdos_renyi(100, 250, rng)
    assert er.num_vertices() == 100 and er.num_arestas() == 250
    ba = prop.gerar_barabasi_albert(100, 2, rng)
    assert ba.num_vertices() == 100 and ba.num_arestas() == 3 + 97 * 2
    ws = prop.gerar_watts_strogatz(100, 4, 0.1, rng)
    assert ws.num_vertices() == 100 and ws.num_arestas() == 200


def test_mv_recupera_gamma_conhecido():
    # amostra de uma lei de potência discreta (zeta) com gamma = 2,5 e kmin = 1
    graus = np.random.default_rng(42).zipf(2.5, size=20000).astype(float)
    r = prop.ajuste_mv(graus, 1)
    assert r["gamma"] == pytest.approx(2.5, abs=3 * r["erro_padrao"] + 0.01)
    assert r["ks"] < 0.02


def test_melhor_kmin_exige_cauda_nao_degenerada():
    graus = np.array([1.0] * 100 + [2.0] * 60 + [3.0] * 50 + [4.0] * 40)
    r = prop.melhor_kmin(graus, cauda_minima=50, valores_minimos=3)
    assert r["kmin"] <= 2


def test_aderencia_reconhece_poisson():
    # grafo G(n, m) grande: graus ~ Poisson, sem subdispersão
    g = prop.gerar_erdos_renyi(3000, 6000, random.Random(1))
    ad = prop.aderencia_graus(g, [list(g.adj)])
    assert ad["indice_dispersao"] == pytest.approx(1, abs=0.1)
    assert ad["poisson"]["p"] > 0.001
    assert sum(ad["observados"]) == 3000
    assert sum(ad["esperados_binomial"]) == pytest.approx(3000)
