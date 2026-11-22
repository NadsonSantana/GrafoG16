"""BFS, DFS e componentes conexas."""

import networkx as nx

from src import algoritmos as alg
from tests.auxiliares import caminho, para_nx


def test_bfs_dfs_caminho():
    g = caminho(1, 2, 3, 4)
    assert alg.bfs(g, 1) == {1: 0, 2: 1, 3: 2, 4: 3}
    assert alg.dfs(g, 1) == [1, 2, 3, 4]


def test_bfs_distancias(grafo):
    h = para_nx(grafo)
    for s in list(grafo.adj)[:10]:
        assert alg.bfs(grafo, s) == nx.single_source_shortest_path_length(h, s)


def test_dfs_visita_a_componente(grafo):
    h = para_nx(grafo)
    for s in list(grafo.adj)[:10]:
        ordem = alg.dfs(grafo, s)
        assert len(ordem) == len(set(ordem))
        assert set(ordem) == nx.node_connected_component(h, s)


def test_componentes(grafo):
    esperado = sorted(sorted(c) for c in nx.connected_components(para_nx(grafo)))
    assert sorted(sorted(c) for c in alg.componentes_conexas(grafo)) == esperado


def test_varreduras_completas(grafo):
    assert len(alg.bfs_completa(grafo)) == len(alg.dfs_completa(grafo)) == grafo.num_vertices()


def test_grafo_real(real):
    assert real.num_vertices() == 451
    assert len(alg.componentes_conexas(real)) == 62
    assert len(alg.bfs_completa(real)) == len(alg.dfs_completa(real)) == 451
