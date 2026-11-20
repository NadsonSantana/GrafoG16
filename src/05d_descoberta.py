"""Etapa 5d - Parte III, questao 4: descoberta mais interessante.

Hipotese: a atencao em cada partida se concentra em um "jogador focal"
(quem todos olham acima do uniforme), e esse papel e persistente no tempo.

Analises (versao dirigida e dados temporais brutos):
  1. Partidas com um jogador cujo grau de entrada dirigido e n-1
     (todos os demais o olham acima de 1/(n-1)).
  2. Concentracao da atencao recebida x emitida (indice de Gini por partida).
  3. Reciprocidade dos arcos x esperado ao acaso (densidade dirigida).
  4. Persistencia: cada partida dividida em 3 tercos; o jogador mais olhado
     e o mesmo nos 3 tercos? Comparado com o nulo de permutacao (rotulos
     embaralhados por terco), e correlacao de Spearman entre 1o e 3o terco.
"""
import re

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

from comum import (AZUL, AZUL_CLARO, DIR_FIG, DIR_PROC, LARANJA, RAIZ, TEXTO, TEXTO_2,
                   carregar_dirigido, estilo)

SEMENTE = 224215945
DIR_RAW = RAIZ / "data" / "raw" / "comm-f2f-Resistance"
COL = re.compile(r"^P(\d+)_TO_P(\d+)$")
TERCOS = 3
PERMUTACOES = 10_000


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if x.sum() == 0:
        return 0.0
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


def matrizes_por_terco(rede, n):
    b = pd.read_csv(DIR_RAW / "network" / f"network{rede}.csv")
    partes = np.array_split(np.arange(len(b)), TERCOS)
    Ms = []
    for idx in partes:
        M = np.zeros((n, n))
        bloco = b.iloc[idx]
        for c in b.columns:
            m = COL.match(c)
            if not m:
                continue
            i, j = int(m.group(1)), int(m.group(2))
            if i <= n and j <= n and i != j:
                M[i - 1, j - 1] = bloco[c].mean()
        Ms.append(M)
    return Ms


def main():
    rng = np.random.default_rng(SEMENTE)
    D = carregar_dirigido()
    pares = pd.read_csv(DIR_PROC / "pares_atencao.csv")
    lista = pd.read_csv(DIR_RAW / "network_list.csv")

    # Matriz de atencao agregada por partida (a partir de pares_atencao.csv)
    linhas = []
    for rede, n in lista.itertuples(index=False):
        A = np.zeros((n, n))
        for r in pares[pares.rede == rede].itertuples(index=False):
            i = int(r.u.split("_P")[1]) - 1
            j = int(r.v.split("_P")[1]) - 1
            A[i, j], A[j, i] = r.atencao_uv, r.atencao_vu
        recebida = A.sum(axis=0)  # atencao que cada jogador recebe
        emitida = A.sum(axis=1)   # fracao do tempo olhando para colegas
        nos = [f"R{rede}_P{p}" for p in range(1, n + 1)]
        grau_ent = [D.in_degree(v) for v in nos]
        focal = int(np.argmax(recebida))

        Ms = matrizes_por_terco(rede, n)
        tops = [int(np.argmax(M.sum(axis=0))) for M in Ms]
        rec_t = [M.sum(axis=0) for M in Ms]
        rho = stats.spearmanr(rec_t[0], rec_t[-1]).statistic

        linhas.append({
            "rede": rede, "n": n,
            "max_grau_entrada": max(grau_ent),
            "tem_focal_total": max(grau_ent) == n - 1,
            "focal": nos[focal], "focal_grau_entrada": grau_ent[focal],
            "parcela_focal": recebida[focal] / recebida.sum(),
            "parcela_uniforme": 1 / n,
            "gini_recebida": gini(recebida), "gini_emitida": gini(emitida),
            "top_tercos": "-".join(str(t + 1) for t in tops),
            "mesmo_top_3_tercos": len(set(tops)) == 1,
            "top_terco_1_igual_3": tops[0] == tops[-1],
            "spearman_terco1_terco3": rho,
        })
    df = pd.DataFrame(linhas)
    df.to_csv(DIR_PROC / "parte3_descoberta_partidas.csv", index=False)

    # Nulo de permutacao para a persistencia: em cada terco o "top" e um
    # jogador sorteado uniformemente -> P(mesmo top nos 3) = 1/n^2 por partida.
    esperado_3 = (1 / df.n ** 2).sum()
    esperado_13 = (1 / df.n).sum()
    sims = np.array([sum(rng.integers(n) == rng.integers(n) == rng.integers(n) for n in df.n)
                     for _ in range(PERMUTACOES)])
    p_3 = (np.sum(sims >= df.mesmo_top_3_tercos.sum()) + 1) / (PERMUTACOES + 1)

    # Reciprocidade x acaso
    arcos = D.number_of_edges()
    possiveis = sum(n * (n - 1) for n in lista.NUMBER_OF_PARTICIPANTS)
    rec = nx.reciprocity(D)

    # Gini recebida x emitida (teste pareado de Wilcoxon)
    w = stats.wilcoxon(df.gini_recebida, df.gini_emitida, alternative="greater")

    resumo = {
        "partidas": len(df),
        "partidas_com_focal_total (grau entrada = n-1)": int(df.tem_focal_total.sum()),
        "parcela_media_focal": df.parcela_focal.mean(),
        "parcela_media_uniforme": df.parcela_uniforme.mean(),
        "razao_parcela_focal_uniforme": (df.parcela_focal / df.parcela_uniforme).mean(),
        "gini_recebida_media": df.gini_recebida.mean(),
        "gini_emitida_media": df.gini_emitida.mean(),
        "wilcoxon_gini_rec>emit_p": w.pvalue,
        "mesmo_top_3_tercos_obs": int(df.mesmo_top_3_tercos.sum()),
        "mesmo_top_3_tercos_esperado_acaso": esperado_3,
        "mesmo_top_3_tercos_p_permutacao": p_3,
        "top_terco1_igual_terco3_obs": int(df.top_terco_1_igual_3.sum()),
        "top_terco1_igual_terco3_esperado": esperado_13,
        "spearman_terco1_terco3_mediana": df.spearman_terco1_terco3.median(),
        "reciprocidade_observada": rec,
        "reciprocidade_esperada_acaso (densidade)": arcos / possiveis,
    }
    pd.DataFrame(list(resumo.items()), columns=["medida", "valor"]) \
        .to_csv(DIR_PROC / "parte3_descoberta_resumo.csv", index=False)

    # Figura
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), dpi=150)
    ax = axes[0]
    ordem = df.sort_values("parcela_focal", ascending=False).reset_index(drop=True)
    x = np.arange(len(ordem))
    ax.bar(x, ordem.parcela_focal, color=[AZUL if f else AZUL_CLARO for f in ordem.tem_focal_total],
           width=0.85, edgecolor="white", linewidth=0.5)
    ax.plot(x, ordem.parcela_uniforme, color=TEXTO, linestyle="--", linewidth=1.2,
            drawstyle="steps-mid", zorder=3)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color=AZUL, label=f"todos o olham acima do uniforme ({int(df.tem_focal_total.sum())})"),
        Patch(color=AZUL_CLARO, label=f"demais partidas ({int((~df.tem_focal_total).sum())})"),
        Line2D([], [], color=TEXTO, linestyle="--", label="parcela uniforme 1/n"),
    ], fontsize=8, frameon=False, loc="upper right")
    ax.set_ylim(0, ordem.parcela_focal.max() * 1.25)
    ax.set_xticks([])
    ax.set_xlabel("62 partidas (ordenadas)", color=TEXTO)
    ax.set_ylabel("parcela da atenção recebida\npelo jogador mais olhado", color=TEXTO)
    ax.set_title("Jogador focal por partida", color=TEXTO, fontsize=10, loc="left")
    estilo(ax)

    ax = axes[1]
    obs = int(df.mesmo_top_3_tercos.sum())
    ax.hist(sims, bins=np.arange(sims.max() + 2) - 0.5, color=AZUL_CLARO, edgecolor="white",
            density=True, label="nulo: top sorteado ao acaso")
    ax.axvline(obs, color=LARANJA, linewidth=2.5, label=f"observado = {obs}")
    ax.set_xlim(-0.5, max(obs, sims.max()) + 2)
    ax.set_xlabel("partidas com o mesmo jogador mais olhado nos 3 terços", color=TEXTO)
    ax.set_ylabel("probabilidade", color=TEXTO)
    ax.set_title(f"Persistência do foco no tempo (p = {p_3:.1e}, {PERMUTACOES:,} simulações)"
                 .replace(",", "."),
                 color=TEXTO, fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False)
    estilo(ax)
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte3_descoberta.png")
    plt.close(fig)

    for k, v in resumo.items():
        print(f"  {k:48s} {v:.4f}" if isinstance(v, float) else f"  {k:48s} {v}")


if __name__ == "__main__":
    main()
