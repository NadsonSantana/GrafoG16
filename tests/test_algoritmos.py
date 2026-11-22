"""Testes dos algoritmos e métricas. O NetworkX é usado apenas como oráculo de referência."""

import math
import random

import networkx as nx
import pytest

from src import algoritmos as alg
from src import metricas as met
from src import propriedades as prop
from src.benchmark import intervalo_confianca
from src.grafo import Grafo


def para_nx(g):
    h = nx.DiGraph() if g.dirigido else nx.Graph()
    h.add_nodes_from(g.adj)
    h.add_weighted_edges_from(g.arestas())
    return h


def grafo_aleatorio(n, p, semente, dirigido=False, desconexo=True):
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


@pytest.fixture(params=[1, 2, 3])
def grafo(request):
    return grafo_aleatorio(40, 0.15, request.param)


@pytest.fixture
def real():
    return Grafo.de_csv("data/processed/arestas.csv", "data/processed/vertices.csv")


# ---------------------------------------------------------------- casos conhecidos

def test_caminho_simples():
    g = Grafo()
    for u, v in [(1, 2), (2, 3), (3, 4)]:
        g.adicionar_aresta(u, v)
    assert alg.bfs(g, 1) == {1: 0, 2: 1, 3: 2, 4: 3}
    assert alg.dfs(g, 1) == [1, 2, 3, 4]
    assert alg.verificar_euleriano(g)["tipo"] == "caminho"
    assert met.excentricidades(g) == {1: 3, 2: 2, 3: 2, 4: 3}


def test_ciclo_euleriano_e_triangulo():
    g = Grafo()
    for u, v in [(1, 2), (2, 3), (3, 1)]:
        g.adicionar_aresta(u, v)
    assert alg.verificar_euleriano(g)["tipo"] == "circuito"
    assert met.numero_triangulos(g) == 1
    assert met.clusterizacao_media(g) == 1.0


def test_estrela_nao_euleriana():
    g = Grafo()
    for v in range(1, 5):
        g.adicionar_aresta(0, v)
    assert alg.verificar_euleriano(g)["tipo"] == "nao_euleriano"


def test_bellman_ford_peso_negativo():
    g = Grafo(dirigido=True)
    g.adicionar_aresta("a", "b", 4)
    g.adicionar_aresta("a", "c", 2)
    g.adicionar_aresta("c", "b", -3)
    dist, ciclo = alg.bellman_ford(g, "a")
    assert dist["b"] == -1 and not ciclo
    g.adicionar_aresta("b", "a", 0.5)  # ciclo a-c-b-a com peso -0.5
    assert alg.bellman_ford(g, "a")[1]


def test_tarjan_caso_conhecido():
    g = Grafo(dirigido=True)
    for u, v in [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 4), (6, 6)]:
        g.adicionar_aresta(u, v)
    sccs = sorted(sorted(c) for c in alg.tarjan(g))
    assert sccs == [[1, 2, 3], [4, 5], [6]]


# ---------------------------------------------------------------- comparação com oráculo

def test_componentes(grafo):
    esperado = sorted(sorted(c) for c in nx.connected_components(para_nx(grafo)))
    assert sorted(sorted(c) for c in alg.componentes_conexas(grafo)) == esperado


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


def test_tarjan_aleatorio():
    for semente in range(5):
        g = grafo_aleatorio(50, 0.05, semente, dirigido=True)
        esperado = sorted(sorted(c) for c in nx.strongly_connected_components(para_nx(g)))
        assert sorted(sorted(c) for c in alg.tarjan(g)) == esperado


def test_prim_kruskal(grafo):
    peso_nx = sum(d["weight"] for *_, d in nx.minimum_spanning_edges(para_nx(grafo), data=True))
    pr, kr = alg.prim(grafo), alg.kruskal(grafo)
    n_comp = len(alg.componentes_conexas(grafo))
    assert len(pr) == len(kr) == grafo.num_vertices() - n_comp
    assert math.isclose(sum(w for *_, w in pr), peso_nx)
    assert math.isclose(sum(w for *_, w in kr), peso_nx)


def test_metricas(grafo):
    h = para_nx(grafo)
    assert math.isclose(met.densidade(grafo), nx.density(h))
    assert math.isclose(met.clusterizacao_media(grafo), nx.average_clustering(h))
    assert met.numero_triangulos(grafo) == sum(nx.triangles(h).values()) // 3
    maior = max(nx.connected_components(h), key=len)
    lcc = h.subgraph(maior)
    assert met.caminhos_medios(grafo.subgrafo(maior)) == pytest.approx(
        nx.average_shortest_path_length(lcc))
    exc = met.excentricidades(grafo.subgrafo(maior))
    assert max(exc.values()) == nx.diameter(lcc)
    assert min(exc.values()) == nx.radius(lcc)


def test_betweenness(grafo):
    esperado = nx.betweenness_centrality(para_nx(grafo), normalized=False)
    obtido = prop.betweenness(grafo)
    for v in esperado:
        assert obtido[v] == pytest.approx(esperado[v])


def test_pontes_articulacoes():
    for semente in range(5):
        g = grafo_aleatorio(60, 0.06, semente)
        h = para_nx(g)
        pontes, art = alg.pontes_articulacoes(g)
        assert {frozenset(p) for p in pontes} == {frozenset(p) for p in nx.bridges(h)}
        assert set(art) == set(nx.articulation_points(h))


def test_pagerank():
    for dirigido in (False, True):
        g = grafo_aleatorio(40, 0.1, 7, dirigido=dirigido)
        h = para_nx(g)
        esperado = nx.pagerank(h, weight=None, tol=1e-10, max_iter=1000)
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


# ---------------------------------------------------------------- grafo real

def test_grafo_real(real):
    assert real.num_vertices() == 451
    assert len(alg.componentes_conexas(real)) == 62
    assert len(alg.bfs_completa(real)) == len(alg.dfs_completa(real)) == 451
    assert len(alg.kruskal(real)) == 451 - 62


def test_intervalo_confianca():
    ic = intervalo_confianca([1.0, 2.0, 3.0, 4.0, 5.0])
    assert ic["distribuicao"] == "t" and ic["media"] == 3.0
    assert ic["ic_inf"] < 3.0 < ic["ic_sup"]
    assert intervalo_confianca([1.0, 2.0] * 15)["distribuicao"] == "Z"
