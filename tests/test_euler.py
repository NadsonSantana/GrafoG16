"""Verificação de Eulerianidade."""

from src import algoritmos as alg
from src.grafo import Grafo
from tests.auxiliares import caminho


def test_caminho_euleriano():
    assert alg.verificar_euleriano(caminho(1, 2, 3, 4))["tipo"] == "caminho"


def test_circuito_euleriano():
    assert alg.verificar_euleriano(caminho(1, 2, 3, 1))["tipo"] == "circuito"


def test_estrela_nao_euleriana():
    g = Grafo()
    for v in range(1, 5):
        g.adicionar_aresta(0, v)
    r = alg.verificar_euleriano(g)
    assert r["tipo"] == "nao_euleriano" and r["vertices_grau_impar"] == 4


def test_desconexo_nao_euleriano():
    g = caminho(1, 2, 3, 1)
    for u, v in [(4, 5), (5, 6), (6, 4)]:
        g.adicionar_aresta(u, v)
    r = alg.verificar_euleriano(g)
    assert r["tipo"] == "nao_euleriano" and not r["conexo"]
