"""Sensibilidade das conclusões ao limiar de atenção do tratamento.

Reconstrói o grafo a partir de data/processed/atencao_pares.csv (atenção de todos os
pares ordenados, antes do limiar) para vários limiares e recalcula as medidas que
sustentam as respostas da Parte III. Não precisa dos dados brutos."""

import csv
import os
import random

import numpy as np

from .algoritmos import componentes_conexas
from .comum import DADOS
from .grafo import Grafo
from .metricas import caminhos_medios, clusterizacao_media, excentricidades
from .propriedades import gerar_er_por_blocos
from .tratamento import arestas_partida

# (rótulo, limiar em função do número n de alvos da partida)
LIMIARES = [
    ("0,5/N", lambda n: 0.5 / n),
    ("1/N (adotado)", lambda n: 1.0 / n),
    ("1,5/N", lambda n: 1.5 / n),
    ("2/N", lambda n: 2.0 / n),
    ("f ≥ 0,10", lambda n: 0.10),
]


def matrizes_atencao(caminho=None):
    """{partida: matriz n x n de atenção f(u, v)} a partir de atencao_pares.csv."""
    caminho = caminho or os.path.join(DADOS, "atencao_pares.csv")
    matrizes = {}
    with open(caminho, newline="") as f:
        for linha in csv.DictReader(f):
            p, n = int(linha["partida"]), int(linha["participantes"])
            m = matrizes.setdefault(p, np.zeros((n, n)))
            m[int(linha["origem"]) - 1, int(linha["destino"]) - 1] = float(linha["atencao"])
    return matrizes


def grafos_por_limiar(matrizes, limiar_de_n):
    """Grafo não-direcionado e número de arcos para um limiar (mesmo tratamento do pipeline)."""
    g = Grafo()
    arcos = 0
    for p, atencao in sorted(matrizes.items()):
        n = atencao.shape[0]
        for i in range(n):
            g.adicionar_vertice(f"g{p:02d}_p{i + 1}")
        a, e = arestas_partida(p, atencao, limiar_de_n(n))
        arcos += len(a)
        for x in e:
            g.adicionar_aresta(x["origem"], x["destino"], x["distancia"])
    return g, arcos


def medidas(g, arcos, repeticoes=10, semente=42):
    """Medidas do grafo e sigma contra o aleatório que preserva as partidas."""
    rng = random.Random(semente)
    n, m = g.num_vertices(), g.num_arestas()
    comps = componentes_conexas(g)
    nao_triviais = [c for c in comps if len(c) > 1]
    # blocos = partidas (com limiares altos, uma partida pode se partir em várias componentes)
    partidas = {}
    for v in g.adj:
        partidas.setdefault(v.split("_")[0], []).append(v)
    c_obs, l_obs = clusterizacao_media(g), caminhos_medios(g)
    cs, ls = [], []
    for _ in range(repeticoes):
        r = gerar_er_por_blocos(g, list(partidas.values()), rng)
        cs.append(clusterizacao_media(r))
        ls.append(caminhos_medios(r))
    c_r, l_r = float(np.mean(cs)), float(np.mean(ls))
    recip = 2 * (arcos - m)  # arcos cujo inverso também existe (arcos - m = pares recíprocos)
    return {
        "arcos": arcos, "arestas": m, "grau_medio": 2 * m / n,
        "isolados": sum(1 for v in g.adj if not g.adj[v]),
        "componentes": len(comps),
        "diametro_max": max(max(excentricidades(g, c).values()) for c in nao_triviais),
        "C": c_obs, "L": l_obs,
        "reciprocidade": recip / arcos if arcos else 0.0,
        "sigma_blocos": (c_obs / c_r) / (l_obs / l_r) if c_r and l_r else float("nan"),
    }


def sensibilidade_limiar(repeticoes=10):
    matrizes = matrizes_atencao()
    linhas = []
    for rotulo, limiar in LIMIARES:
        g, arcos = grafos_por_limiar(matrizes, limiar)
        linhas.append({"limiar": rotulo, **medidas(g, arcos, repeticoes)})
    return linhas
