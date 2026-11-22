"""Parte I: métricas estruturais, distribuição de graus e visualização do grafo.

Uso: python -m src.parte1 [--so-figuras]"""

import os

import networkx as nx  # usado apenas para o layout da visualização

from . import algoritmos as alg
from .comum import (AZUL, DADOS, TINTA_2, carregar, eixo_graus_log, ler_json, plt, preparar,
                    salvar_figura, salvar_json, so_figuras)
from .grafo import Grafo
from .metricas import metricas_estruturais


def parte1(g):
    m = metricas_estruturais(g)
    tamanhos = m["tamanhos_componentes"]
    m["tamanhos_componentes_freq"] = {t: tamanhos.count(t) for t in sorted(set(tamanhos))}
    salvar_json("parte1_metricas.json", m)
    figura_graus(m["distribuicao_graus"], g.num_vertices())
    figura_visualizacao(g)
    return m


def figura_graus(dist, n):
    ks = list(dist)
    pk = [dist[k] / n for k in ks]
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.4, 2.4))
    a.bar(ks, [dist[k] for k in ks], color=AZUL, width=0.7)
    for k in ks:
        a.text(k, dist[k] + 2, str(dist[k]), ha="center", fontsize=7, color=TINTA_2)
    a.set_xlabel("grau $k$")
    a.set_ylabel("número de vértices")
    a.set_title("Escala linear", fontsize=9)
    b.loglog(ks, pk, "o", color=AZUL, markersize=5)
    eixo_graus_log(b, ks)
    b.set_xlabel("grau $k$ (log)")
    b.set_ylabel("$P(k)$ (log)")
    b.set_title("Escala log-log", fontsize=9)
    fig.tight_layout()
    salvar_figura(fig, "distribuicao_graus")


def figura_visualizacao(g):
    """Grafo completo (62 partidas) e zoom em uma partida com espessura = atenção."""
    h = nx.Graph()
    h.add_nodes_from(g.adj)
    h.add_edges_from((u, v) for u, v, _ in g.arestas())
    comps = alg.componentes_conexas(g)
    comps.sort(key=lambda c: min(c))
    pos = {}
    colunas = 16
    for i, comp in enumerate(comps):
        bloco = nx.Graph()
        bloco.add_nodes_from(sorted(comp))  # ordem fixa: o layout não depende de hash
        bloco.add_edges_from(h.subgraph(comp).edges)
        sub = nx.spring_layout(bloco, seed=1)
        cx, cy = (i % colunas) * 3.0, -(i // colunas) * 3.0
        for v, (x, y) in sub.items():
            pos[v] = (cx + x, cy + y)
    fig, ax = plt.subplots(figsize=(7.0, 2.0))
    nx.draw_networkx_edges(h, pos, ax=ax, width=0.4, edge_color="#9a9890", alpha=0.7)
    nx.draw_networkx_nodes(h, pos, ax=ax, node_size=6, node_color=AZUL, linewidths=0)
    ax.set_axis_off()
    salvar_figura(fig, "grafo_completo")

    # zoom: partida 0 na versão dirigida ponderada
    dirigidas = Grafo.de_csv(os.path.join(DADOS, "arestas_dirigidas.csv"), dirigido=True,
                             coluna_peso="atencao")
    sub = dirigidas.subgrafo([v for v in dirigidas.adj if v.startswith("g00_")])
    d = nx.DiGraph()
    d.add_weighted_edges_from(sub.arestas())
    posz = nx.circular_layout(sorted(d.nodes))
    fig, ax = plt.subplots(figsize=(3.2, 3.0))
    larguras = [0.5 + 6 * d[u][v]["weight"] for u, v in d.edges]
    nx.draw_networkx_edges(d, posz, ax=ax, width=larguras, edge_color="#9a9890",
                           arrowsize=8, connectionstyle="arc3,rad=0.12", node_size=420)
    nx.draw_networkx_nodes(d, posz, ax=ax, node_size=420, node_color=AZUL, linewidths=0)
    nx.draw_networkx_labels(d, posz, ax=ax, labels={v: v.split("_")[1].upper() for v in d},
                            font_size=7, font_color="white")
    ax.set_axis_off()
    salvar_figura(fig, "partida0")


def figuras(g):
    """Redesenha as figuras da Parte I a partir de results/parte1_metricas.json."""
    m = ler_json("parte1_metricas.json")
    figura_graus({int(k): v for k, v in m["distribuicao_graus"].items()}, m["vertices"])
    figura_visualizacao(g)


def main():
    preparar()
    g, _ = carregar()
    if so_figuras():
        print("Parte I - redesenhando figuras")
        figuras(g)
        return
    print("Parte I - métricas estruturais")
    parte1(g)


if __name__ == "__main__":
    main()
