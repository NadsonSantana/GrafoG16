"""Testes do tratamento: atenção, limiar 1/N, simetrização, pesos e dados tratados."""

import csv
import json
import math
import os

import numpy as np
import pytest

from src.tratamento import DESTINO, arestas_partida, matriz_atencao, nome_vertice


def partida_sintetica(alvos):
    """Monta o array T x n x (n+1) a partir de alvos[t][i] = índice do alvo de i no
    instante t (0 = laptop, j = participante j)."""
    t, n = len(alvos), len(alvos[0])
    dados = np.zeros((t, n, n + 1), dtype=int)
    for k, linha in enumerate(alvos):
        for i, alvo in enumerate(linha):
            dados[k, i, alvo] = 1
    return dados


# 3 participantes, 4 instantes. Linha k: alvo de P1, P2, P3 no instante k.
ALVOS = [
    [2, 1, 0],
    [2, 1, 1],
    [2, 3, 1],
    [0, 1, 2],
]


def test_nome_vertice():
    assert nome_vertice(5, 3) == "g05_p3"


def test_matriz_atencao_frações_e_diagonal():
    f = matriz_atencao(partida_sintetica(ALVOS))
    # P1 olha P2 em 3 de 4 instantes; P2 olha P1 em 3 de 4; P3 olha P1 em 2 e P2 em 1.
    esperado = np.array([[0, 0.75, 0], [0.75, 0, 0.25], [0.5, 0.25, 0]])
    assert np.allclose(f, esperado)
    assert np.all(np.diag(f) == 0)


def test_atencao_soma_no_maximo_um():
    f = matriz_atencao(partida_sintetica(ALVOS))
    assert np.all(f.sum(axis=1) <= 1 + 1e-12)  # o restante é atenção ao laptop


def test_limiar_um_sobre_n():
    # n = 3 participantes -> N = 3 alvos (2 colegas + laptop) -> limiar 1/3
    arcos, _ = arestas_partida(0, matriz_atencao(partida_sintetica(ALVOS)))
    pares = {(a["origem"], a["destino"]) for a in arcos}
    assert pares == {("g00_p1", "g00_p2"), ("g00_p2", "g00_p1"), ("g00_p3", "g00_p1")}
    # P2 -> P3 (0,25) e P3 -> P2 (0,25) ficam abaixo de 1/3


def test_limiar_inclusivo():
    f = np.array([[0, 0.5], [0.5, 0]])  # n = 2 -> limiar 1/2, exatamente atingido
    arcos, arestas = arestas_partida(1, f)
    assert len(arcos) == 2 and len(arestas) == 1


def test_simetrizacao_e_pesos():
    _, arestas = arestas_partida(0, matriz_atencao(partida_sintetica(ALVOS)))
    por_par = {(a["origem"], a["destino"]): a for a in arestas}
    assert set(por_par) == {("g00_p1", "g00_p2"), ("g00_p1", "g00_p3")}
    # {P1, P3} existe porque P3 -> P1 passa o limiar, mesmo sem P1 -> P3
    a13 = por_par[("g00_p1", "g00_p3")]
    assert a13["atencao"] == pytest.approx((0 + 0.5) / 2)
    assert a13["distancia"] == pytest.approx(-math.log(0.25), abs=1e-6)
    assert por_par[("g00_p1", "g00_p2")]["atencao"] == pytest.approx(0.75)


def test_distancias_positivas():
    _, arestas = arestas_partida(0, matriz_atencao(partida_sintetica(ALVOS)))
    assert all(a["distancia"] > 0 for a in arestas)


# ---------------------------------------------------------------- dados tratados

def ler_csv(nome):
    with open(os.path.join(DESTINO, nome), newline="") as f:
        return list(csv.DictReader(f))


def test_etapas_conferem_com_snap():
    with open(os.path.join(DESTINO, "etapas_tratamento.json")) as f:
        e = json.load(f)
    # valores publicados na página do dataset no SNAP
    assert e["vertices"] == 451 and e["partidas"] == 62
    assert e["arestas_temporais"] + e["arestas_temporais_laptop"] == 3_126_993
    assert e["pares_apos_limiar"] == len(ler_csv("arestas_dirigidas.csv"))
    assert e["arestas_nao_dirigidas"] == len(ler_csv("arestas.csv"))


def test_arestas_tratadas_consistentes():
    vertices = {v["id"] for v in ler_csv("vertices.csv")}
    arestas = ler_csv("arestas.csv")
    pares = set()
    for a in arestas:
        u, v = a["origem"], a["destino"]
        assert u in vertices and v in vertices and u < v  # sem laços, cada par uma vez
        assert u.split("_")[0] == v.split("_")[0]  # nunca entre partidas
        assert 0 < float(a["atencao"]) < 1
        assert float(a["distancia"]) == pytest.approx(-math.log(float(a["atencao"])), abs=1e-4)
        pares.add((u, v))
    assert len(pares) == len(arestas)
    # toda aresta não-direcionada vem de ao menos um arco dirigido
    arcos = {(a["origem"], a["destino"]) for a in ler_csv("arestas_dirigidas.csv")}
    assert all((u, v) in arcos or (v, u) in arcos for u, v in pares)


def test_atencao_de_todos_os_pares_gera_os_arcos():
    pares = ler_csv("atencao_pares.csv")
    assert len(pares) == 2870  # todos os pares ordenados intra-partida
    acima = {(nome_vertice(int(p["partida"]), int(p["origem"])),
              nome_vertice(int(p["partida"]), int(p["destino"])))
             for p in pares if float(p["atencao"]) * int(p["participantes"]) >= 1 - 1e-9}
    arcos = {(a["origem"], a["destino"]) for a in ler_csv("arestas_dirigidas.csv")}
    assert acima == arcos
