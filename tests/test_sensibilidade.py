"""Sensibilidade ao limiar, réplicas das partidas (escala) e robustez."""

import math
import random

from src.comum import carregar
from src.parte2 import replicar_partidas
from src.propriedades import robustez
from src.sensibilidade import grafos_por_limiar, matrizes_atencao, medidas


def test_limiar_adotado_reproduz_o_grafo_tratado(real):
    g, arcos = grafos_por_limiar(matrizes_atencao(), lambda n: 1.0 / n)
    assert g.num_vertices() == real.num_vertices() == 451
    assert g.num_arestas() == real.num_arestas() and arcos == 1241
    assert set(map(frozenset, ((u, v) for u, v, _ in g.arestas()))) == \
        set(map(frozenset, ((u, v) for u, v, _ in real.arestas())))


def test_reciprocidade_e_limiares_monotonos():
    m = matrizes_atencao()
    g1, a1 = grafos_por_limiar(m, lambda n: 1.0 / n)
    assert round(medidas(g1, a1, repeticoes=1)["reciprocidade"], 3) == 0.467
    g2, _ = grafos_por_limiar(m, lambda n: 2.0 / n)
    assert set(map(frozenset, ((u, v) for u, v, _ in g2.arestas()))) <= \
        set(map(frozenset, ((u, v) for u, v, _ in g1.arestas())))


def test_replicas_preservam_estrutura_das_partidas():
    g, d = carregar()
    h, hd = replicar_partidas(g, d, 124, random.Random(0))
    assert h.num_vertices() > 2 * 451 - 60 and hd.num_vertices() == h.num_vertices()
    graus_copias = sorted(len(viz) for v, viz in h.adj.items() if v.startswith("c"))
    assert max(graus_copias) <= 7
    # até 62 partidas: subgrafo das primeiras partidas reais
    sub, _ = replicar_partidas(g, d, 8, random.Random(0))
    assert all(not v.startswith("c") for v in sub.adj)


def test_robustez_usa_teto_de_cinco_por_cento(real):
    r = robustez(real, repeticoes=2, fracoes_curva=(0.10,))
    assert r["r"] == math.ceil(0.05 * 451) == 23
    assert r["curva"]["fracoes"] == [0.0, 0.10]
    assert r["curva"]["aleatoria"][0]["A_media"] == 1.0
