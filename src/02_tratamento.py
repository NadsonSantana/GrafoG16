"""Etapa 2 - Tratamento do G16 (comm-f2f-Resistance).

Pipeline:
  1. Temporal -> estatico: para cada partida e cada par dirigido (u, v), a
     atencao a_uv e a fracao dos quadros (1/3 s) em que u olha para v (versao
     binaria do SNAP). A versao ponderada (probabilidade media) e guardada
     como atributo auxiliar.
  2. Remocao do "LAPTOP" (nao e participante) e de auto-loops (inexistentes).
  3. Limiar relativo: o arco u->v e mantido se a_uv >= 1/(n-1), isto e, se u
     olha para v mais do que olharia sob atencao uniforme entre os n-1 colegas.
     O limiar se ajusta ao tamanho da partida (5 a 8 jogadores).
  4. Dirigido -> nao dirigido: atencao mutua m_uv = (a_uv + a_vu) / 2 com o
     mesmo limiar 1/(n-1). Distancia da aresta = 1 - m_uv (atencao alta =
     vertices "proximos"), usada nos algoritmos de caminho minimo e AGM.
  5. Vertices recebem rotulo global "R<rede>_P<participante>"; todos os 451
     sao mantidos (inclusive isolados). Partidas nao compartilham jogadores,
     logo o grafo e desconexo por construcao.

Saidas em data/processed/ e figures/.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR_RAW = RAIZ / "data" / "raw" / "comm-f2f-Resistance"
DIR_PROC = RAIZ / "data" / "processed"
DIR_FIG = RAIZ / "figures"

COL = re.compile(r"^P(\d+)_TO_P(\d+)$")
FATORES_SENSIBILIDADE = [0.5, 0.75, 1.0, 1.25, 1.5]

# Paleta de referencia (skill dataviz): serie 1 e tintas de texto
AZUL = "#2a78d6"
TEXTO = "#0b0b0b"
TEXTO_2 = "#52514e"


def rotulo(rede, p):
    return f"R{rede}_P{p}"


def atencao_partida(rede, n):
    """Matrizes n x n de atencao media (binaria e ponderada) entre participantes."""
    b = pd.read_csv(DIR_RAW / "network" / f"network{rede}.csv")
    w = pd.read_csv(DIR_RAW / "network" / f"network{rede}_weighted.csv")
    A_bin = np.zeros((n, n))
    A_w = np.zeros((n, n))
    for c in b.columns:
        m = COL.match(c)
        if not m:
            continue  # TIME e colunas *_TO_LAPTOP
        i, j = int(m.group(1)), int(m.group(2))
        if i > n or j > n or i == j:
            continue
        A_bin[i - 1, j - 1] = b[c].mean()
        A_w[i - 1, j - 1] = w[c].mean()
    return A_bin, A_w


def construir(lista, matrizes, fator=1.0):
    """Grafo dirigido e nao dirigido para um limiar fator * 1/(n-1)."""
    D, G = nx.DiGraph(), nx.Graph()
    for rede, n in lista.itertuples(index=False):
        A_bin, A_w = matrizes[rede]
        tau = fator / (n - 1)
        for p in range(1, n + 1):
            for H in (D, G):
                H.add_node(rotulo(rede, p), rede=rede, participante=p)
        for i in range(n):
            for j in range(n):
                if i != j and A_bin[i, j] >= tau:
                    D.add_edge(rotulo(rede, i + 1), rotulo(rede, j + 1),
                               atencao=A_bin[i, j], prob=A_w[i, j])
            for j in range(i + 1, n):
                m = (A_bin[i, j] + A_bin[j, i]) / 2
                if m >= tau:
                    G.add_edge(rotulo(rede, i + 1), rotulo(rede, j + 1),
                               atencao=m,
                               prob=(A_w[i, j] + A_w[j, i]) / 2,
                               distancia=1 - m)
    return D, G


def resumo(G):
    comps = list(nx.connected_components(G))
    return {
        "vertices": G.number_of_nodes(),
        "arestas": G.number_of_edges(),
        "grau_medio": 2 * G.number_of_edges() / G.number_of_nodes(),
        "isolados": nx.number_of_isolates(G),
        "componentes": len(comps),
        "maior_componente": max(len(c) for c in comps),
    }


def figura_atencao(pares):
    razao = pares["atencao_mutua"] / pares["uniforme"]
    fig, ax = plt.subplots(figsize=(7, 4), dpi=150)
    ax.hist(razao, bins=40, color=AZUL, edgecolor="white", linewidth=1)
    ax.axvline(1.0, color=TEXTO_2, linestyle="--", linewidth=1.5)
    ax.text(1.03, ax.get_ylim()[1] * 0.92,
            f"limiar = 1/(n-1)\n{(razao >= 1).mean():.0%} dos pares mantidos",
            color=TEXTO, fontsize=9, va="top")
    ax.set_xlabel("Atenção mútua / atenção uniforme  [m_uv · (n−1)]", color=TEXTO)
    ax.set_ylabel("Pares de participantes", color=TEXTO)
    ax.set_title("Distribuição da atenção mútua agregada (1.435 pares, 62 partidas)",
                 color=TEXTO, fontsize=11, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(TEXTO_2)
    ax.spines["bottom"].set_color(TEXTO_2)
    ax.tick_params(colors=TEXTO_2)
    ax.grid(axis="y", color="#e5e4e0", linewidth=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(DIR_FIG / "distribuicao_atencao.png")
    plt.close(fig)


def main():
    DIR_PROC.mkdir(parents=True, exist_ok=True)
    DIR_FIG.mkdir(parents=True, exist_ok=True)
    lista = pd.read_csv(DIR_RAW / "network_list.csv")
    matrizes = {rede: atencao_partida(rede, n) for rede, n in lista.itertuples(index=False)}

    # Tabela de todos os pares (antes do limiar), para a figura e o relatorio
    linhas = []
    for rede, n in lista.itertuples(index=False):
        A_bin, A_w = matrizes[rede]
        for i in range(n):
            for j in range(i + 1, n):
                linhas.append({
                    "rede": rede, "n": n, "u": rotulo(rede, i + 1), "v": rotulo(rede, j + 1),
                    "atencao_uv": A_bin[i, j], "atencao_vu": A_bin[j, i],
                    "atencao_mutua": (A_bin[i, j] + A_bin[j, i]) / 2,
                    "prob_mutua": (A_w[i, j] + A_w[j, i]) / 2,
                    "uniforme": 1 / (n - 1),
                })
    pares = pd.DataFrame(linhas)
    pares.to_csv(DIR_PROC / "pares_atencao.csv", index=False)
    figura_atencao(pares)

    # Sensibilidade ao limiar
    sens = []
    for f in FATORES_SENSIBILIDADE:
        D, G = construir(lista, matrizes, f)
        sens.append({"fator_limiar": f, "arcos_dirigidos": D.number_of_edges(), **resumo(G)})
    sens = pd.DataFrame(sens)
    sens.to_csv(DIR_PROC / "sensibilidade_limiar.csv", index=False)

    # Grafo final (fator 1.0)
    D, G = construir(lista, matrizes, 1.0)
    nos = pd.DataFrame([{"no": v, **d} for v, d in G.nodes(data=True)])
    nos.to_csv(DIR_PROC / "g16_nos.csv", index=False)
    pd.DataFrame([{"u": u, "v": v, **d} for u, v, d in G.edges(data=True)]) \
        .to_csv(DIR_PROC / "g16_nao_dirigido.csv", index=False)
    pd.DataFrame([{"u": u, "v": v, **d} for u, v, d in D.edges(data=True)]) \
        .to_csv(DIR_PROC / "g16_dirigido.csv", index=False)

    pd.set_option("display.width", 200)
    print("Atencao mutua agregada (todos os pares):")
    print(pares[["atencao_mutua", "prob_mutua"]].describe().round(4).T)
    print(f"correlacao binaria x ponderada: {pares.atencao_mutua.corr(pares.prob_mutua):.3f}")
    print("\nSensibilidade ao limiar (fator x 1/(n-1)):")
    print(sens.round(2).to_string(index=False))
    r = resumo(G)
    print("\nGrafo final nao dirigido:", r)
    print("Grafo final dirigido:", D.number_of_nodes(), "vertices,", D.number_of_edges(), "arcos")
    print("Arquivos: data/processed/{g16_nos,g16_nao_dirigido,g16_dirigido,pares_atencao,"
          "sensibilidade_limiar}.csv, figures/distribuicao_atencao.png")


if __name__ == "__main__":
    main()
