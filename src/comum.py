"""Utilitarios compartilhados pelos scripts da Parte III."""
from pathlib import Path

import networkx as nx
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR_PROC = RAIZ / "data" / "processed"
DIR_FIG = RAIZ / "figures"

# Paleta de referencia (skill dataviz)
AZUL = "#2a78d6"
AZUL_CLARO = "#9ec5f4"
LARANJA = "#eb6834"
AQUA = "#1baf7a"
AMARELO = "#eda100"
TEXTO = "#0b0b0b"
TEXTO_2 = "#52514e"
GRADE = "#e5e4e0"


def carregar_grafo():
    """G16 tratado, nao dirigido (com os 451 vertices, inclusive isolados)."""
    nos = pd.read_csv(DIR_PROC / "g16_nos.csv")
    arestas = pd.read_csv(DIR_PROC / "g16_nao_dirigido.csv")
    G = nx.Graph()
    for r in nos.itertuples(index=False):
        G.add_node(r.no, rede=r.rede, participante=r.participante)
    for r in arestas.itertuples(index=False):
        G.add_edge(r.u, r.v, atencao=r.atencao, prob=r.prob, distancia=r.distancia)
    return G


def carregar_dirigido():
    nos = pd.read_csv(DIR_PROC / "g16_nos.csv")
    arcos = pd.read_csv(DIR_PROC / "g16_dirigido.csv")
    D = nx.DiGraph()
    for r in nos.itertuples(index=False):
        D.add_node(r.no, rede=r.rede, participante=r.participante)
    for r in arcos.itertuples(index=False):
        D.add_edge(r.u, r.v, atencao=r.atencao, prob=r.prob)
    return D


def adjacencia(G):
    """Lista de adjacencia como dict de conjuntos (independente do NetworkX)."""
    return {v: set(G.neighbors(v)) for v in G}


def estilo(ax, grade="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(TEXTO_2)
    ax.tick_params(colors=TEXTO_2, labelsize=8)
    if grade:
        ax.grid(axis=grade, color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)
