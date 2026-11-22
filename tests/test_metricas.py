"""Métricas estruturais da Parte I."""

import math

import networkx as nx
import pytest

from src import metricas as met
from tests.auxiliares import caminho, para_nx


def test_excentricidades_caminho():
    assert met.excentricidades(caminho(1, 2, 3, 4)) == {1: 3, 2: 2, 3: 2, 4: 3}


def test_triangulo():
    g = caminho(1, 2, 3, 1)
    assert met.numero_triangulos(g) == 1
    assert met.clusterizacao_media(g) == 1.0
    assert met.transitividade(g) == 1.0


def test_metricas(grafo):
    h = para_nx(grafo)
    assert math.isclose(met.densidade(grafo), nx.density(h))
    assert math.isclose(met.clusterizacao_media(grafo), nx.average_clustering(h))
    assert math.isclose(met.transitividade(grafo), nx.transitivity(h))
    assert met.numero_triangulos(grafo) == sum(nx.triangles(h).values()) // 3
    maior = max(nx.connected_components(h), key=len)
    lcc = h.subgraph(maior)
    assert met.caminhos_medios(grafo.subgrafo(maior)) == pytest.approx(
        nx.average_shortest_path_length(lcc))
    exc = met.excentricidades(grafo.subgrafo(maior))
    assert max(exc.values()) == nx.diameter(lcc)
    assert min(exc.values()) == nx.radius(lcc)


def test_metricas_do_grafo_real(real):
    m = met.metricas_estruturais(real)
    assert (m["vertices"], m["arestas"], m["num_componentes"]) == (451, 951, 62)
    assert m["grau_min"] == 1 and m["grau_max"] == 7
