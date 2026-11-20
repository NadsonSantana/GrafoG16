"""Etapa 5a - Parte III, questao 1: o G16 tem propriedade small-world?

Implementacao propria (sem funcoes prontas de small-world, clustering ou
caminho medio): coeficiente de clusterizacao local, caminho medio por BFS e
gerador Erdos-Renyi G(n, m).

Criterio (Apendice A): L ~ L_rand e C >> C_rand, resumido por
sigma = (C / C_rand) / (L / L_rand); sigma > 1 indica small-world.

Duas comparacoes:
  (1) Por componente (politica adotada): cada componente nao trivial e
      comparada com grafos ER conexos de mesmo n e m (amostragem por rejeicao).
  (2) Global "ingenua": grafo inteiro x ER G(451, 548), para mostrar por que
      ela e enganosa num grafo formado por 62 partidas desconexas.
"""
import random
from collections import deque
from itertools import combinations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from comum import (AZUL, AZUL_CLARO, DIR_FIG, DIR_PROC, LARANJA, TEXTO, TEXTO_2,
                   adjacencia, carregar_grafo, estilo)

SEMENTE = 224215945
AMOSTRAS_COMP = 200
AMOSTRAS_GLOBAL = 50
TENTATIVAS_CONEXO = 5000


# ------------------------------------------------------- medidas proprias
def clusterizacao_local(adj, v):
    viz = adj[v]
    k = len(viz)
    if k < 2:
        return 0.0
    ligacoes = sum(1 for a, b in combinations(viz, 2) if b in adj[a])
    return 2 * ligacoes / (k * (k - 1))


def clusterizacao_media(adj):
    """Media sobre todos os vertices (grau < 2 conta como 0, convencao usual)."""
    return sum(clusterizacao_local(adj, v) for v in adj) / len(adj)


def distancias_bfs(adj, origem):
    dist = {origem: 0}
    fila = deque([origem])
    while fila:
        u = fila.popleft()
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                fila.append(v)
    return dist


def caminho_medio(adj):
    """Media das distancias entre pares (ordenados) conectados e distintos."""
    soma, pares = 0, 0
    for s in adj:
        d = distancias_bfs(adj, s)
        soma += sum(d.values())
        pares += len(d) - 1
    return soma / pares if pares else float("nan"), pares


def er_gnm(n, m, rng):
    """Erdos-Renyi G(n, m): m arestas sorteadas uniformemente sem reposicao."""
    todos = list(combinations(range(n), 2))
    adj = {i: set() for i in range(n)}
    for u, v in rng.sample(todos, m):
        adj[u].add(v)
        adj[v].add(u)
    return adj


def conexo(adj):
    return len(distancias_bfs(adj, next(iter(adj)))) == len(adj)


def er_conexo(n, m, rng):
    for _ in range(TENTATIVAS_CONEXO):
        adj = er_gnm(n, m, rng)
        if conexo(adj):
            return adj
    return None


# ------------------------------------------------------- analises
def por_componente(G, rng):
    linhas = []
    for c in sorted(nx.connected_components(G), key=lambda c: (-len(c), min(c))):
        n = len(c)
        if n < 2:
            continue
        H = G.subgraph(c)
        m = H.number_of_edges()
        adj = adjacencia(H)
        rede = G.nodes[next(iter(c))]["rede"]
        C = clusterizacao_media(adj)
        L, _ = caminho_medio(adj)
        linha = {"rede": rede, "n": n, "m": m, "densidade": 2 * m / (n * (n - 1)),
                 "C": C, "L": L}
        completo = m == n * (n - 1) // 2
        if n < 4 or completo:
            linha["situacao"] = "excluida (n < 4 ou completa: ER idêntico)"
            linhas.append(linha)
            continue
        Cs, Ls = [], []
        for _ in range(AMOSTRAS_COMP):
            R = er_conexo(n, m, rng)
            if R is None:
                break
            Cs.append(clusterizacao_media(R))
            Ls.append(caminho_medio(R)[0])
        if len(Cs) < AMOSTRAS_COMP:
            linha["situacao"] = "excluida (ER conexo raro)"
            linhas.append(linha)
            continue
        Cs, Ls = np.array(Cs), np.array(Ls)
        C_r, L_r = Cs.mean(), Ls.mean()
        linha.update({
            "situacao": "analisada", "C_rand": C_r, "C_rand_dp": Cs.std(ddof=1),
            "L_rand": L_r, "L_rand_dp": Ls.std(ddof=1),
            "C_sobre_Crand": C / C_r if C_r > 0 else np.nan,
            "L_sobre_Lrand": L / L_r,
            "sigma": (C / C_r) / (L / L_r) if C_r > 0 else np.nan,
            # p-valor empirico unilateral: fracao dos ER com C >= C observado
            "p_C_maior": (np.sum(Cs >= C) + 1) / (len(Cs) + 1),
            "p_L_menor": (np.sum(Ls <= L) + 1) / (len(Ls) + 1),
        })
        linhas.append(linha)
    return pd.DataFrame(linhas)


def global_ingenua(G, rng):
    adj = adjacencia(G)
    n, m = G.number_of_nodes(), G.number_of_edges()
    C = clusterizacao_media(adj)
    L, _ = caminho_medio(adj)
    Cs, Ls = [], []
    for _ in range(AMOSTRAS_GLOBAL):
        R = er_gnm(n, m, rng)
        Cs.append(clusterizacao_media(R))
        Ls.append(caminho_medio(R)[0])
    C_r, L_r = np.mean(Cs), np.mean(Ls)
    k = 2 * m / n
    return {"n": n, "m": m, "grau_medio": k, "C": C, "L": L,
            "C_rand": C_r, "L_rand": L_r,
            "C_rand_teorico": k / n, "L_rand_teorico": np.log(n) / np.log(k),
            "sigma": (C / C_r) / (L / L_r)}


# ------------------------------------------------------- figura
def figura(df):
    a = df[df.situacao == "analisada"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), dpi=150)

    ax = axes[0]
    ax.scatter(a.C_rand, a.C, s=28, color=AZUL, edgecolor="white", linewidth=0.8, zorder=3)
    lim = [0, max(a.C.max(), a.C_rand.max()) * 1.08]
    ax.plot(lim, lim, color=TEXTO_2, linestyle="--", linewidth=1)
    ax.set_xlim(lim), ax.set_ylim(lim)
    ax.set_xlabel("C_rand (média de 200 ER conexos)", color=TEXTO)
    ax.set_ylabel("C observado", color=TEXTO)
    ax.set_title("Clusterização por componente", color=TEXTO, fontsize=10, loc="left")
    ax.text(0.97, 0.03, f"acima da diagonal: {(a.C > a.C_rand).sum()}/{len(a)}",
            transform=ax.transAxes, ha="right", fontsize=8, color=TEXTO_2)
    estilo(ax, grade="both")

    ax = axes[1]
    ax.scatter(a.L_rand, a.L, s=28, color=AZUL, edgecolor="white", linewidth=0.8, zorder=3)
    lim = [min(a.L.min(), a.L_rand.min()) * 0.95, max(a.L.max(), a.L_rand.max()) * 1.05]
    ax.plot(lim, lim, color=TEXTO_2, linestyle="--", linewidth=1)
    ax.set_xlim(lim), ax.set_ylim(lim)
    ax.set_xlabel("L_rand (média de 200 ER conexos)", color=TEXTO)
    ax.set_ylabel("L observado", color=TEXTO)
    ax.set_title("Caminho médio por componente", color=TEXTO, fontsize=10, loc="left")
    ax.text(0.97, 0.03, f"acima da diagonal: {(a.L > a.L_rand).sum()}/{len(a)}",
            transform=ax.transAxes, ha="right", fontsize=8, color=TEXTO_2)
    estilo(ax, grade="both")

    ax = axes[2]
    s = a.sigma.dropna()
    ax.hist(s, bins=20, color=AZUL, edgecolor="white", linewidth=1)
    ax.axvline(1, color=TEXTO_2, linestyle="--", linewidth=1.2)
    ax.axvline(s.median(), color=LARANJA, linewidth=2)
    ax.text(s.max(), ax.get_ylim()[1] * 0.95, f"mediana = {s.median():.2f}",
            color=LARANJA, fontsize=8, va="top", ha="right")
    ax.set_xlabel("σ = (C/C_rand) / (L/L_rand)", color=TEXTO)
    ax.set_ylabel("Componentes", color=TEXTO)
    ax.set_title("Índice small-world σ (σ > 1 ⇒ small-world)", color=TEXTO, fontsize=10,
                 loc="left")
    estilo(ax)

    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte3_small_world.png")
    plt.close(fig)


def main():
    rng = random.Random(SEMENTE)
    G = carregar_grafo()

    # Conferencia das implementacoes proprias com o NetworkX
    adj = adjacencia(G)
    assert np.isclose(clusterizacao_media(adj), nx.average_clustering(G))
    H = G.subgraph(max(nx.connected_components(G), key=len))
    assert np.isclose(caminho_medio(adjacencia(H))[0], nx.average_shortest_path_length(H))

    df = por_componente(G, rng)
    df.to_csv(DIR_PROC / "parte3_small_world_componentes.csv", index=False)
    a = df[df.situacao == "analisada"]

    agregado = {
        "componentes_analisadas": len(a),
        "componentes_excluidas": int((df.situacao != "analisada").sum()),
        "C_medio": a.C.mean(), "C_rand_medio": a.C_rand.mean(),
        "L_medio": a.L.mean(), "L_rand_medio": a.L_rand.mean(),
        "densidade_media": a.densidade.mean(),
    }
    agregado["C_sobre_Crand"] = agregado["C_medio"] / agregado["C_rand_medio"]
    agregado["L_sobre_Lrand"] = agregado["L_medio"] / agregado["L_rand_medio"]
    agregado["sigma_agregado"] = agregado["C_sobre_Crand"] / agregado["L_sobre_Lrand"]
    agregado["sigma_mediana"] = a.sigma.median()
    agregado["sigma_maior_1"] = int((a.sigma > 1).sum())
    agregado["C_signif_maior_5pct"] = int((a.p_C_maior < 0.05).sum())
    agregado["L_signif_menor_5pct"] = int((a.p_L_menor < 0.05).sum())

    glob = global_ingenua(G, rng)
    pd.DataFrame([{"analise": "por componente (agregado)", **agregado},
                  {"analise": "global ingênua", **glob}]) \
        .to_csv(DIR_PROC / "parte3_small_world_resumo.csv", index=False)
    figura(df)

    print("Implementações próprias conferidas com NetworkX (C médio e L).\n")
    print("Por componente:")
    print(df.situacao.value_counts().to_string())
    for k, v in agregado.items():
        print(f"  {k:24s} {v:.4f}" if isinstance(v, float) else f"  {k:24s} {v}")
    print("\nσ por componente (analisadas):")
    print(a.sigma.describe().round(3).to_string())
    print("\nGlobal ingênua (grafo inteiro × ER G(451, 548)):")
    for k, v in glob.items():
        print(f"  {k:16s} {v:.4f}" if isinstance(v, float) else f"  {k:16s} {v}")


if __name__ == "__main__":
    main()
