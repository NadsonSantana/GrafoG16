"""Executa todas as análises (Partes I, II e III) e salva resultados e figuras.

Uso: python -m src.executar [--rapido]
"""

import json
import math
import os
import random
import sys

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402  (usado apenas para o layout da visualização)
import numpy as np  # noqa: E402

from . import algoritmos as alg  # noqa: E402
from . import propriedades as prop  # noqa: E402
from .benchmark import medir  # noqa: E402
from .grafo import Grafo  # noqa: E402
from .metricas import metricas_estruturais  # noqa: E402

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DADOS = os.path.join(RAIZ, "data", "processed")
RESULTADOS = os.path.join(RAIZ, "results")
FIGURAS = os.path.join(RAIZ, "figures")

AZUL, LARANJA, AQUA, VERMELHO = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
TINTA, TINTA_2, GRADE = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.edgecolor": TINTA_2,
    "axes.labelcolor": TINTA, "xtick.color": TINTA_2, "ytick.color": TINTA_2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.6,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.dpi": 200,
})


def eixo_graus_log(ax, ks):
    """Marcas inteiras legíveis no eixo x logarítmico de graus."""
    ax.set_xticks(list(ks), [str(int(k)) for k in ks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def carregar():
    g = Grafo.de_csv(os.path.join(DADOS, "arestas.csv"), os.path.join(DADOS, "vertices.csv"))
    d = Grafo.de_csv(os.path.join(DADOS, "arestas_dirigidas.csv"),
                     os.path.join(DADOS, "vertices.csv"), dirigido=True)
    return g, d


def salvar_json(nome, dados):
    with open(os.path.join(RESULTADOS, nome), "w") as f:
        json.dump(dados, f, indent=2, ensure_ascii=False, default=str)


# ---------------------------------------------------------------- Parte I

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
    fig.savefig(os.path.join(FIGURAS, "distribuicao_graus.pdf"))
    plt.close(fig)


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
        sub = nx.spring_layout(h.subgraph(comp), seed=1)
        cx, cy = (i % colunas) * 3.0, -(i // colunas) * 3.0
        for v, (x, y) in sub.items():
            pos[v] = (cx + x, cy + y)
    fig, ax = plt.subplots(figsize=(7.0, 2.0))
    nx.draw_networkx_edges(h, pos, ax=ax, width=0.4, edge_color="#9a9890", alpha=0.7)
    nx.draw_networkx_nodes(h, pos, ax=ax, node_size=6, node_color=AZUL, linewidths=0)
    ax.set_axis_off()
    fig.savefig(os.path.join(FIGURAS, "grafo_completo.pdf"))
    fig.savefig(os.path.join(FIGURAS, "grafo_completo.png"))
    plt.close(fig)

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
    fig.savefig(os.path.join(FIGURAS, "partida0.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- Parte II

def parte2(g, d, rapido=False):
    n_rapido = 5 if rapido else 50
    n_bf = 5 if rapido else 30
    n_fw = 3 if rapido else 10
    casos = [
        ("BFS", "O(V+E)", lambda: alg.bfs_completa(g), n_rapido),
        ("DFS", "O(V+E)", lambda: alg.dfs_completa(g), n_rapido),
        ("Eulerianidade", "O(V+E)", lambda: alg.verificar_euleriano(g), n_rapido),
        ("Dijkstra (todas as fontes)", "O(V(V+E)\\log V)", lambda: alg.dijkstra_todos(g), n_rapido),
        ("Bellman-Ford (todas as fontes)", "O(V^2 E)", lambda: alg.bellman_ford_todos(g), n_bf),
        ("Floyd-Warshall", "O(V^3)", lambda: alg.floyd_warshall(g), n_fw),
        ("Tarjan", "O(V+E)", lambda: alg.tarjan(d), n_rapido),
        ("Prim", "O(E\\log V)", lambda: alg.prim(g), n_rapido),
        ("Kruskal", "O(E\\log V)", lambda: alg.kruskal(g), n_rapido),
    ]
    tempos = {}
    saidas = {}
    for nome, complexidade, funcao, n in casos:
        resultado, est = medir(funcao, n)
        est["complexidade"] = complexidade
        tempos[nome] = est
        saidas[nome] = resultado
        print(f"  {nome:32s} n={n:3d} media={est['media']:.6f}s")

    # resumo dos resultados de cada algoritmo
    euler_por_partida = [alg.verificar_euleriano(g.subgrafo(c))["tipo"]
                         for c in alg.componentes_conexas(g)]
    vs, fw = saidas["Floyd-Warshall"]
    idx = {v: i for i, v in enumerate(vs)}
    dj = saidas["Dijkstra (todas as fontes)"]
    bf = saidas["Bellman-Ford (todas as fontes)"]
    concordam = all(math.isclose(dj[s][t], fw[idx[s]][idx[t]], abs_tol=1e-9)
                    and math.isclose(bf[s][t], dj[s][t], abs_tol=1e-9)
                    for s in dj for t in dj[s])
    finitos = [fw[i][j] for i in range(len(vs)) for j in range(len(vs))
               if i != j and fw[i][j] < math.inf]
    sccs = saidas["Tarjan"]
    pr, kr = saidas["Prim"], saidas["Kruskal"]
    # exemplo de caminho mínimo: maior distância finita encontrada
    i_max, j_max = max(((i, j) for i in range(len(vs)) for j in range(len(vs))
                        if i != j and fw[i][j] < math.inf), key=lambda p: fw[p[0]][p[1]])
    resumo = {
        "bfs_visitados": len(saidas["BFS"]),
        "dfs_visitados": len(saidas["DFS"]),
        "euler_global": saidas["Eulerianidade"],
        "euler_por_partida": {t: euler_por_partida.count(t) for t in set(euler_por_partida)},
        "pares_alcancaveis": len(finitos),
        "distancia_media": float(np.mean(finitos)),
        "distancia_max": float(max(finitos)),
        "par_distancia_max": [vs[i_max], vs[j_max]],
        "dijkstra_bf_fw_concordam": concordam,
        "bellman_ford_ciclo_negativo": False,
        "num_scc": len(sccs),
        "scc_tamanhos": {t: [len(c) for c in sccs].count(t)
                         for t in sorted({len(c) for c in sccs})},
        "scc_nao_triviais": sum(1 for c in sccs if len(c) > 1),
        "partidas_fortemente_conexas": sum(
            1 for c in sccs if len(c) > 1 and len(c) == sum(
                1 for v in g.adj if v.split("_")[0] == c[0].split("_")[0])),
        "prim_arestas": len(pr), "prim_peso": sum(w for *_, w in pr),
        "kruskal_arestas": len(kr), "kruskal_peso": sum(w for *_, w in kr),
    }
    salvar_json("parte2_tempos.json", tempos)
    salvar_json("parte2_resultados.json", resumo)
    figura_tempos(tempos)
    return tempos, resumo


def figura_tempos(tempos):
    nomes = list(tempos)
    medias = [tempos[k]["media"] * 1000 for k in nomes]
    erros = [(tempos[k]["ic_sup"] - tempos[k]["media"]) * 1000 for k in nomes]
    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    y = np.arange(len(nomes))
    ax.barh(y, medias, xerr=erros, color=AZUL, height=0.6,
            error_kw={"ecolor": TINTA, "elinewidth": 0.8, "capsize": 2})
    ax.set_yticks(y, [n.replace(" (todas as fontes)", "*") for n in nomes])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("tempo médio (ms, escala log) com IC 95%")
    ax.grid(axis="y", visible=False)
    for yi, m, e in zip(y, medias, erros):
        ax.text((m + e) * 1.2, yi, f"{m:.3g} ms", va="center", fontsize=7, color=TINTA_2)
    ax.set_xlim(right=max(medias) * 12)
    fig.savefig(os.path.join(FIGURAS, "tempos.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- Parte III

def parte3(g, rapido=False):
    blocos = alg.componentes_conexas(g)
    reps = 5 if rapido else 30
    sw = prop.small_world(g, blocos, repeticoes=reps)
    lp = prop.lei_de_potencia(g)
    rob = prop.robustez(g, repeticoes=reps)
    modelos = prop.comparar_modelos(g, repeticoes=3 if rapido else 10)
    por_partida = small_world_por_partida(g, blocos)
    salvar_json("parte3_small_world.json", {"global": sw, "por_partida": por_partida})
    salvar_json("parte3_lei_potencia.json", lp)
    salvar_json("parte3_robustez.json", rob)
    salvar_json("parte3_modelos.json", modelos)
    figura_lei_potencia(lp)
    figura_robustez(rob)
    return sw, lp, rob, modelos


def small_world_por_partida(g, blocos):
    """C e L de cada partida (componente) isoladamente."""
    from .metricas import caminhos_medios, clusterizacao_media, densidade
    linhas = []
    for c in blocos:
        s = g.subgrafo(c)
        linhas.append({"partida": c[0].split("_")[0], "n": len(c), "m": s.num_arestas(),
                       "densidade": densidade(s), "C": clusterizacao_media(s),
                       "L": caminhos_medios(s)})
    return {
        "densidade_media": float(np.mean([x["densidade"] for x in linhas])),
        "C_medio": float(np.mean([x["C"] for x in linhas])),
        "L_medio": float(np.mean([x["L"] for x in linhas])),
        "partidas_completas": sum(1 for x in linhas if x["densidade"] == 1.0),
        "linhas": linhas,
    }


def figura_lei_potencia(lp):
    ks = np.array(lp["ks"])
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.4, 2.5))
    a.loglog(ks, lp["pk"], "o", color=AZUL, markersize=5, label="observado")
    xs = np.linspace(ks.min(), ks.max(), 50)
    a.loglog(xs, lp["constante"] * xs ** (-lp["gamma_mq"]), "-", color=LARANJA, linewidth=2,
             label=f"ajuste $\\propto k^{{{-lp['gamma_mq']:.2f}}}$ ($R^2$={lp['r2']:.2f})")
    a.loglog(ks, [lp["poisson"][str(int(k))] if str(int(k)) in lp["poisson"]
                  else lp["poisson"][int(k)] for k in ks],
             "s--", color=AQUA, markersize=4, linewidth=1.2, label="Poisson($\\langle k\\rangle$)")
    a.set_xlabel("grau $k$ (log)")
    a.set_ylabel("$P(k)$ (log)")
    a.legend(fontsize=7, loc="lower right")
    eixo_graus_log(a, ks)
    cc = lp["ccdf"]
    kc = sorted(int(k) for k in cc)
    b.loglog(kc, [cc[k] if k in cc else cc[str(k)] for k in kc], "o-", color=AZUL,
             markersize=5, linewidth=2)
    b.set_xlabel("grau $k$ (log)")
    b.set_ylabel("$P(K \\geq k)$ (log)")
    b.set_title("CCDF", fontsize=9)
    eixo_graus_log(b, kc)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURAS, "lei_potencia.pdf"))
    plt.close(fig)


def figura_robustez(rob):
    metricas = [("A_partidas", "A: fração no núcleo\nda partida"), ("B_componentes", "B: nº de\ncomponentes"),
                ("C_eficiencia", "C: eficiência\nglobal"), ("D_isolados", "D: fração de\nisolados")]
    fig, eixos = plt.subplots(1, 4, figsize=(6.6, 2.4))
    for ax, (chave, titulo) in zip(eixos, metricas):
        valores = [r[chave] for r in rob["aleatoria"]]
        ax.boxplot(valores, widths=0.5, patch_artist=True,
                   boxprops={"facecolor": "#cfe0f5", "edgecolor": AZUL},
                   medianprops={"color": AZUL, "linewidth": 2},
                   whiskerprops={"color": AZUL}, capprops={"color": AZUL},
                   flierprops={"markeredgecolor": AZUL, "markersize": 3})
        media, dp = np.mean(valores), np.std(valores, ddof=1)
        ax.errorbar(1.35, media, yerr=dp, fmt="o", color=TINTA, markersize=4, capsize=2,
                    elinewidth=0.8)
        ax.plot(1, rob["original"][chave], "_", color=TINTA_2, markersize=18, markeredgewidth=1.5)
        ax.plot(1.75, rob["grau"][chave], "D", color=VERMELHO, markersize=6)
        ax.plot(2.05, rob["betweenness"][chave], "^", color=LARANJA, markersize=7)
        ax.set_xlim(0.5, 2.4)
        ax.set_xticks([1, 1.9], ["aleat.", "direc."])
        ax.set_title(titulo, fontsize=8)
        ax.grid(axis="x", visible=False)
    handles = [
        plt.Line2D([], [], color=TINTA_2, marker="_", linestyle="", markersize=12, label="original"),
        plt.Line2D([], [], color=TINTA, marker="o", linestyle="", label="média ± dp (aleatória)"),
        plt.Line2D([], [], color=VERMELHO, marker="D", linestyle="", label="ataque por grau"),
        plt.Line2D([], [], color=LARANJA, marker="^", linestyle="", label="ataque por betweenness"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7, bbox_to_anchor=(0.5, -0.08))
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURAS, "robustez.pdf"))
    plt.close(fig)


def main():
    rapido = "--rapido" in sys.argv
    random.seed(42)
    os.makedirs(RESULTADOS, exist_ok=True)
    os.makedirs(FIGURAS, exist_ok=True)
    g, d = carregar()
    print("Parte I - métricas estruturais")
    parte1(g)
    print("Parte II - algoritmos")
    parte2(g, d, rapido)
    print("Parte III - propriedades da rede")
    parte3(g, rapido)
    print(f"Resultados em {RESULTADOS} e figuras em {FIGURAS}")


if __name__ == "__main__":
    main()
