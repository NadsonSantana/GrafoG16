"""Etapa 5c - Parte III, questao 3: robustez do G16 (Apendice C).

Politica (grafo desconexo, maior componente com 8 vertices): a analise e
feita no grafo inteiro (N = 451); remover 5% = r = 23 vertices.

Metricas apos remover vertices e arestas incidentes:
  A  tamanho da maior componente (absoluto e relativo ao original, S/S0)
  B  numero de componentes
  C  distancias: caminho medio entre pares conectados e eficiencia global
     E = media de 1/d(u,v) sobre todos os pares (1/inf = 0)
  D  fracao de vertices isolados (entre os remanescentes)
  E* (extra) fracao de pares conectados preservada em relacao ao original,
     mais sensivel que A num grafo de componentes pequenas.

(a) Remocao aleatoria: 50 repeticoes, r vertices uniformes.
(b) Ataque direcionado com ranking fixo calculado no grafo original:
    grau (hubs) e betweenness (vertices ponte). Ambos sao comparados com a
    distribuicao aleatoria (z-score e percentil empirico).
Extra: curvas de degradacao para fracoes removidas de 0 a 30%.
"""
import math
import random
from collections import deque

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from comum import (AQUA, AZUL, AZUL_CLARO, DIR_FIG, DIR_PROC, LARANJA, TEXTO, TEXTO_2,
                   carregar_grafo, estilo)

SEMENTE = 224215945
FRACAO = 0.05
REPETICOES = 50
FRACOES_CURVA = np.round(np.arange(0, 0.301, 0.025), 3)
REPETICOES_CURVA = 30

METRICAS = [
    ("A_maior_componente", "A · maior componente (vértices)"),
    ("A_relativo", "A · S / S₀"),
    ("B_componentes", "B · nº de componentes"),
    ("C_caminho_medio", "C · caminho médio (pares conectados)"),
    ("C_eficiencia", "C · eficiência global"),
    ("D_isolados", "D · fração de isolados"),
    ("E_pares_conectados", "E* · pares conectados / original"),
]


def bfs_dist(adj, s):
    dist = {s: 0}
    fila = deque([s])
    while fila:
        u = fila.popleft()
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                fila.append(v)
    return dist


def metricas(G, S0, pares0):
    adj = {v: set(G.neighbors(v)) for v in G}
    n = len(adj)
    comps = list(nx.connected_components(G))
    S = max(len(c) for c in comps)
    soma_d, soma_inv, pares = 0, 0.0, 0
    for s in adj:
        for t, d in bfs_dist(adj, s).items():
            if d > 0:
                soma_d += d
                soma_inv += 1 / d
                pares += 1
    return {
        "A_maior_componente": S,
        "A_relativo": S / S0,
        "B_componentes": len(comps),
        "C_caminho_medio": soma_d / pares if pares else float("nan"),
        "C_eficiencia": soma_inv / (n * (n - 1)),
        "D_isolados": sum(1 for v in adj if not adj[v]) / n,
        "E_pares_conectados": (pares // 2) / pares0,
    }


def pares_conectados(G):
    return sum(len(c) * (len(c) - 1) // 2 for c in nx.connected_components(G))


def remover(G, alvos):
    H = G.copy()
    H.remove_nodes_from(alvos)
    return H


def rankings(G):
    bc = nx.betweenness_centrality(G, normalized=True)
    grau = dict(G.degree())
    por_grau = sorted(G, key=lambda v: (-grau[v], -bc[v], v))
    por_bc = sorted(G, key=lambda v: (-bc[v], -grau[v], v))
    return por_grau, por_bc, grau, bc


# ------------------------------------------------------- figuras
def figura_boxplots(aleat, original, ataques):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), dpi=150)
    estilos = {"grau": (LARANJA, "D"), "betweenness": (AQUA, "^")}
    for ax, (chave, titulo) in zip(axes.flat, METRICAS):
        dados = aleat[chave].to_numpy()
        ax.boxplot(dados, positions=[0], widths=0.5, patch_artist=True, showfliers=True,
                   boxprops=dict(facecolor=AZUL_CLARO, edgecolor=AZUL),
                   medianprops=dict(color=AZUL, linewidth=2),
                   whiskerprops=dict(color=AZUL), capprops=dict(color=AZUL),
                   flierprops=dict(marker="o", markersize=3, markerfacecolor=AZUL,
                                   markeredgecolor="white"))
        mu, sd = dados.mean(), dados.std(ddof=1)
        ax.errorbar([0.42], [mu], yerr=[sd], fmt="o", color=TEXTO, ms=5, capsize=4,
                    elinewidth=1.2)
        ax.axhline(original[chave], color=TEXTO_2, linestyle="--", linewidth=1)
        for i, (nome, vals) in enumerate(ataques.items()):
            cor, marc = estilos[nome]
            ax.scatter([1 + 0.5 * i], [vals[chave]], color=cor, marker=marc, s=60, zorder=3,
                       edgecolor="white", linewidth=0.8)
        ax.set_xticks([0, 1, 1.5])
        ax.set_xticklabels(["aleatória", "grau", "betw."], fontsize=8)
        ax.set_xlim(-0.5, 2)
        ax.set_title(titulo, color=TEXTO, fontsize=9, loc="left")
        estilo(ax)
    leg = axes.flat[-1]
    leg.axis("off")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    leg.legend(handles=[
        Patch(facecolor=AZUL_CLARO, edgecolor=AZUL, label=f"remoção aleatória ({REPETICOES} rep.)"),
        Line2D([], [], color=TEXTO, marker="o", linestyle="none", label="média ± desvio padrão"),
        Line2D([], [], color=LARANJA, marker="D", linestyle="none", label="ataque: top-5% grau"),
        Line2D([], [], color=AQUA, marker="^", linestyle="none",
               label="ataque: top-5% betweenness"),
        Line2D([], [], color=TEXTO_2, linestyle="--", label="grafo original"),
    ], loc="center left", frameon=False, fontsize=9)
    fig.suptitle("Robustez do G16: remoção de 5% dos vértices (r = 23 de 451)", color=TEXTO,
                 fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte3_robustez_boxplots.png")
    plt.close(fig)


def figura_curvas(curva):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), dpi=150)
    paineis = [("E_pares_conectados", "E* · pares conectados / original"),
               ("C_eficiencia", "C · eficiência global"),
               ("D_isolados", "D · fração de isolados")]
    cores = {"aleatória": AZUL, "grau": LARANJA, "betweenness": AQUA}
    for ax, (chave, titulo) in zip(axes, paineis):
        for estrat, cor in cores.items():
            s = curva[curva.estrategia == estrat].sort_values("fracao")
            x = s.fracao * 100
            ax.plot(x, s[f"{chave}_media"], color=cor, linewidth=2, marker="o", ms=3,
                    label=estrat)
            if estrat == "aleatória":
                ax.fill_between(x, s[f"{chave}_media"] - s[f"{chave}_dp"],
                                s[f"{chave}_media"] + s[f"{chave}_dp"], color=cor, alpha=0.15,
                                linewidth=0)
        ax.axvline(5, color=TEXTO_2, linestyle=":", linewidth=1)
        ax.set_xlabel("% de vértices removidos", color=TEXTO)
        ax.set_title(titulo, color=TEXTO, fontsize=10, loc="left")
        estilo(ax, grade="both")
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("Degradação progressiva (ranking fixo; faixa = ±1 dp da remoção aleatória, "
                 f"{REPETICOES_CURVA} rep.)", color=TEXTO, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte3_robustez_curvas.png")
    plt.close(fig)


def main():
    rng = random.Random(SEMENTE)
    G = carregar_grafo()
    N = G.number_of_nodes()
    r = math.ceil(FRACAO * N)
    S0 = max(len(c) for c in nx.connected_components(G))
    pares0 = pares_conectados(G)
    original = metricas(G, S0, pares0)

    vertices = list(G)
    aleat = pd.DataFrame([metricas(remover(G, rng.sample(vertices, r)), S0, pares0)
                          for _ in range(REPETICOES)])
    aleat.to_csv(DIR_PROC / "parte3_robustez_aleatoria.csv", index=False)

    por_grau, por_bc, grau, bc = rankings(G)
    ataques = {"grau": metricas(remover(G, por_grau[:r]), S0, pares0),
               "betweenness": metricas(remover(G, por_bc[:r]), S0, pares0)}
    pd.DataFrame([{"vertice": v, "grau": grau[v], "betweenness": bc[v],
                   "rank_grau": por_grau.index(v) + 1, "rank_betweenness": por_bc.index(v) + 1}
                  for v in sorted(set(por_grau[:r]) | set(por_bc[:r]))]) \
        .to_csv(DIR_PROC / "parte3_robustez_alvos.csv", index=False)

    linhas = []
    for chave, titulo in METRICAS:
        a = aleat[chave]
        mu, sd = a.mean(), a.std(ddof=1)
        linha = {"metrica": titulo, "original": original[chave], "aleatoria_media": mu,
                 "aleatoria_dp": sd, "aleatoria_min": a.min(), "aleatoria_max": a.max()}
        for nome, vals in ataques.items():
            v = vals[chave]
            linha[f"ataque_{nome}"] = v
            linha[f"z_{nome}"] = (v - mu) / sd if sd > 0 else float("nan")
            linha[f"percentil_{nome}"] = (a < v).mean() * 100
        linhas.append(linha)
    comp = pd.DataFrame(linhas)
    comp.to_csv(DIR_PROC / "parte3_robustez_comparacao.csv", index=False)
    figura_boxplots(aleat, original, ataques)

    # Curvas de degradacao
    linhas = []
    for f in FRACOES_CURVA:
        k = int(round(f * N))
        amostras = pd.DataFrame([metricas(remover(G, rng.sample(vertices, k)), S0, pares0)
                                 for _ in range(REPETICOES_CURVA if k else 1)])
        linha = {"estrategia": "aleatória", "fracao": f, "removidos": k}
        for chave, _ in METRICAS:
            linha[f"{chave}_media"] = amostras[chave].mean()
            linha[f"{chave}_dp"] = amostras[chave].std(ddof=1) if k else 0.0
        linhas.append(linha)
        for nome, ranking in (("grau", por_grau), ("betweenness", por_bc)):
            m = metricas(remover(G, ranking[:k]), S0, pares0)
            linhas.append({"estrategia": nome, "fracao": f, "removidos": k,
                           **{f"{c}_media": m[c] for c, _ in METRICAS},
                           **{f"{c}_dp": 0.0 for c, _ in METRICAS}})
    curva = pd.DataFrame(linhas)
    curva.to_csv(DIR_PROC / "parte3_robustez_curvas.csv", index=False)
    figura_curvas(curva)

    pd.set_option("display.width", 250)
    print(f"N = {N}, r = {r}, S0 = {S0}, pares conectados originais = {pares0}")
    print(f"Sobreposição top-{r} grau × betweenness: {len(set(por_grau[:r]) & set(por_bc[:r]))}")
    print(f"Graus dos alvos por grau: {sorted((grau[v] for v in por_grau[:r]), reverse=True)}")
    print(f"Graus dos alvos por betweenness: {sorted((grau[v] for v in por_bc[:r]), reverse=True)}")
    print()
    print(comp.round(4).to_string(index=False))
    print("\nCurva (E* pares conectados / original):")
    print(curva.pivot(index="fracao", columns="estrategia", values="E_pares_conectados_media")
          .round(3).to_string())


if __name__ == "__main__":
    main()
