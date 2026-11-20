"""Etapa 5e - Parte III, questao 5 (EXTRA): o G16 se aproxima de algum modelo
classico (Erdos-Renyi, Barabasi-Albert, Watts-Strogatz)?

Como o grafo e a uniao de componentes pequenas (uma por partida), cada
componente nao trivial (n >= 4) e comparada com 200 grafos de cada modelo de
mesmo tamanho:
  ER  G(n, m) com o mesmo numero de arestas;
  BA  m_BA arestas por novo vertice, m_BA escolhido para m_BA(n - m_BA) ~ m;
  WS  anel com k vizinhos (k par mais proximo de 2m/n, minimo 2) e p = 0,1.
Estatisticas por componente: clusterizacao media C, caminho medio L (pares
conectados), desvio padrao do grau e centralizacao de grau de Freeman
(estrela = 1). Para cada modelo: media das estatisticas sobre as componentes,
z-score do observado e distancia de variacao total (TV) entre a distribuicao
de graus agregada observada e a do modelo.
"""
import random
from collections import Counter

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from comum import (AMARELO, AQUA, AZUL, DIR_FIG, DIR_PROC, LARANJA, TEXTO, TEXTO_2,
                   carregar_grafo, estilo)

SEMENTE = 224215945
AMOSTRAS = 200
P_WS = 0.1
MODELOS = ["ER", "BA", "WS"]
CORES = {"observado": AZUL, "ER": LARANJA, "BA": AQUA, "WS": AMARELO}
ESTATS = [("C", "clusterização média C"), ("L", "caminho médio L"),
          ("dp_grau", "desvio padrão do grau"), ("centralizacao", "centralização de grau")]


def estatisticas(H):
    n = H.number_of_nodes()
    graus = np.array([d for _, d in H.degree()])
    soma, pares = 0, 0
    for _, dist in nx.all_pairs_shortest_path_length(H):
        soma += sum(dist.values())
        pares += len(dist) - 1
    return {
        "C": nx.average_clustering(H),
        "L": soma / pares if pares else np.nan,
        "dp_grau": graus.std(),
        "centralizacao": (graus.max() - graus).sum() / ((n - 1) * (n - 2)),
    }, graus


def gerar(modelo, n, m, rng):
    semente = rng.randrange(2 ** 31)
    if modelo == "ER":
        return nx.gnm_random_graph(n, m, seed=semente)
    if modelo == "BA":
        m_ba = min(range(1, n), key=lambda q: (abs(q * (n - q) - m), q))
        return nx.barabasi_albert_graph(n, m_ba, seed=semente)
    k = max(2, 2 * round(m / n))
    k = min(k, n - 1 if (n - 1) % 2 == 0 else n - 2)
    return nx.watts_strogatz_graph(n, k, P_WS, seed=semente)


def distancia_tv(c1, c2):
    t1, t2 = sum(c1.values()), sum(c2.values())
    return 0.5 * sum(abs(c1.get(k, 0) / t1 - c2.get(k, 0) / t2) for k in set(c1) | set(c2))


def main():
    rng = random.Random(SEMENTE)
    G = carregar_grafo()
    comps = [c for c in nx.connected_components(G) if len(c) >= 4]

    obs_linhas, graus_obs = [], Counter()
    mod_linhas = {mo: [] for mo in MODELOS}
    graus_mod = {mo: Counter() for mo in MODELOS}
    arestas_mod = {mo: [] for mo in MODELOS}
    for c in comps:
        H = G.subgraph(c)
        n, m = H.number_of_nodes(), H.number_of_edges()
        e, g = estatisticas(H)
        obs_linhas.append({"n": n, "m": m, **e})
        graus_obs.update(g)
        for mo in MODELOS:
            amostras = []
            for _ in range(AMOSTRAS):
                R = gerar(mo, n, m, rng)
                er, gr = estatisticas(R)
                amostras.append(er)
                graus_mod[mo].update(gr)
                arestas_mod[mo].append(R.number_of_edges() - m)
            am = pd.DataFrame(amostras)
            mod_linhas[mo].append(am.mean().to_dict())

    obs = pd.DataFrame(obs_linhas)
    linhas = [{"modelo": "observado", **obs[[k for k, _ in ESTATS]].mean().to_dict(),
               "dif_media_arestas": 0.0, "TV_graus": 0.0}]
    for mo in MODELOS:
        md = pd.DataFrame(mod_linhas[mo])
        linha = {"modelo": mo, **md.mean().to_dict(),
                 "dif_media_arestas": float(np.mean(arestas_mod[mo])),
                 "TV_graus": distancia_tv(graus_obs, graus_mod[mo])}
        # z pareado: diferenca observada - modelo por componente
        for k, _ in ESTATS:
            d = (obs[k] - md[k]).dropna()
            linha[f"z_{k}"] = d.mean() / (d.std(ddof=1) / np.sqrt(len(d)))
        linhas.append(linha)
    res = pd.DataFrame(linhas)
    res.to_csv(DIR_PROC / "parte3_modelos.csv", index=False)

    # Figura: 4 estatisticas + distribuicao de graus
    fig, axes = plt.subplots(1, 5, figsize=(17, 3.8), dpi=150)
    nomes = ["observado"] + MODELOS
    for ax, (k, titulo) in zip(axes[:4], ESTATS):
        vals = res.set_index("modelo").loc[nomes, k]
        ax.bar(range(4), vals, color=[CORES[n] for n in nomes], width=0.7, edgecolor="white")
        for i, v in enumerate(vals):
            ax.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8, color=TEXTO_2)
        ax.set_xticks(range(4))
        ax.set_xticklabels(nomes, fontsize=8)
        ax.set_title(titulo, color=TEXTO, fontsize=10, loc="left")
        estilo(ax)
    ax = axes[4]
    ks = range(0, max(max(graus_obs), *(max(graus_mod[m]) for m in MODELOS)) + 1)
    tot_o = sum(graus_obs.values())
    ax.plot(ks, [graus_obs.get(k, 0) / tot_o for k in ks], color=CORES["observado"],
            linewidth=2.5, marker="o", ms=4, label="observado", zorder=3)
    for mo in MODELOS:
        t = sum(graus_mod[mo].values())
        tv = res.set_index("modelo").loc[mo, "TV_graus"]
        ax.plot(ks, [graus_mod[mo].get(k, 0) / t for k in ks], color=CORES[mo], linewidth=2,
                marker="o", ms=3, label=f"{mo} (TV = {tv:.2f})")
    ax.set_xlabel("grau k", color=TEXTO)
    ax.set_title("P(k) nas componentes (n ≥ 4)", color=TEXTO, fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False)
    estilo(ax)
    fig.suptitle(f"Observado × modelos clássicos por componente ({len(comps)} componentes, "
                 f"{AMOSTRAS} grafos por modelo; WS com p = {P_WS})", color=TEXTO, fontsize=11,
                 x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte3_modelos.png")
    plt.close(fig)

    pd.set_option("display.width", 220)
    print(f"{len(comps)} componentes com n >= 4 comparadas")
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
