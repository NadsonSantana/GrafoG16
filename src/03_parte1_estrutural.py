"""Etapa 3 - Parte I: analise estrutural obrigatoria do G16 tratado.

Politica para grafo desconexo (85 componentes, maior com 8 vertices):
medidas globais sobre o grafo inteiro; medidas de distancia (diametro, raio,
caminho medio) calculadas dentro de cada componente e agregadas. Distancias
em numero de saltos (grafo nao ponderado).
"""
import math
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR_PROC = RAIZ / "data" / "processed"
DIR_FIG = RAIZ / "figures"

AZUL = "#2a78d6"
TEXTO = "#0b0b0b"
TEXTO_2 = "#52514e"
GRADE = "#e5e4e0"


def carregar_grafo():
    nos = pd.read_csv(DIR_PROC / "g16_nos.csv")
    arestas = pd.read_csv(DIR_PROC / "g16_nao_dirigido.csv")
    G = nx.Graph()
    for r in nos.itertuples(index=False):
        G.add_node(r.no, rede=r.rede, participante=r.participante)
    for r in arestas.itertuples(index=False):
        G.add_edge(r.u, r.v, atencao=r.atencao, prob=r.prob, distancia=r.distancia)
    return G


def estilo(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(TEXTO_2)
    ax.tick_params(colors=TEXTO_2)
    ax.grid(axis="y", color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)


def figura_graus(dist, n):
    ks = list(range(0, max(dist) + 1))
    pk = [dist.get(k, 0) / n for k in ks]
    fig, ax = plt.subplots(figsize=(7, 4), dpi=150)
    barras = ax.bar(ks, pk, color=AZUL, width=0.8, edgecolor="white", linewidth=1)
    for b, k in zip(barras, ks):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.005,
                str(dist.get(k, 0)), ha="center", va="bottom", fontsize=8, color=TEXTO_2)
    ax.set_xticks(ks)
    ax.set_xlabel("Grau k", color=TEXTO)
    ax.set_ylabel("P(k)", color=TEXTO)
    ax.set_title("Distribuição de graus (rótulo = nº de vértices)", color=TEXTO,
                 fontsize=11, loc="left")
    estilo(ax)
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte1_distribuicao_graus.png")
    plt.close(fig)


def figura_componentes(tamanhos):
    cont = Counter(tamanhos)
    ts = list(range(1, max(cont) + 1))
    fig, ax = plt.subplots(figsize=(7, 4), dpi=150)
    barras = ax.bar(ts, [cont.get(t, 0) for t in ts], color=AZUL, width=0.8,
                    edgecolor="white", linewidth=1)
    for b, t in zip(barras, ts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.3,
                str(cont.get(t, 0)), ha="center", va="bottom", fontsize=8, color=TEXTO_2)
    ax.set_xticks(ts)
    ax.set_xlabel("Tamanho da componente (vértices)", color=TEXTO)
    ax.set_ylabel("Número de componentes", color=TEXTO)
    ax.set_title("Tamanho das 85 componentes conexas", color=TEXTO, fontsize=11, loc="left")
    estilo(ax)
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte1_tamanho_componentes.png")
    plt.close(fig)


def figura_grafo(G, comps):
    """Grafo inteiro: cada componente com layout proprio, em grade por tamanho."""
    comps = sorted(comps, key=lambda c: (-len(c), min(c)))
    colunas = 12
    linhas = math.ceil(len(comps) / colunas)
    pos = {}
    for idx, c in enumerate(comps):
        H = G.subgraph(c)
        local = nx.spring_layout(H, seed=42, weight="atencao") if len(c) > 1 else {next(iter(c)): (0, 0)}
        cx, cy = idx % colunas, -(idx // colunas)
        for v, (x, y) in local.items():
            pos[v] = (cx * 2.6 + x, cy * 2.6 + y)
    larg = [0.5 + 4 * G.edges[e]["atencao"] for e in G.edges]
    fig, ax = plt.subplots(figsize=(14, 1.3 * linhas + 1), dpi=150)
    nx.draw_networkx_edges(G, pos, ax=ax, width=larg, edge_color="#9ec5f4")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=18, node_color=AZUL,
                           edgecolors="white", linewidths=0.5)
    ax.set_title("G16 tratado: 451 vértices, 548 arestas, 85 componentes "
                 "(espessura = atenção mútua; ordenado por tamanho)",
                 color=TEXTO, fontsize=11, loc="left")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte1_grafo.png")
    plt.close(fig)


def main():
    G = carregar_grafo()
    n, m = G.number_of_nodes(), G.number_of_edges()
    graus = [d for _, d in G.degree()]
    dist_graus = Counter(graus)
    comps = list(nx.connected_components(G))
    tamanhos = sorted((len(c) for c in comps), reverse=True)

    # Distancias por componente (saltos)
    linhas_comp = []
    soma_dist, pares_conectados = 0, 0
    for c in sorted(comps, key=lambda c: (-len(c), min(c))):
        H = G.subgraph(c)
        k = len(c)
        if k > 1:
            exc = nx.eccentricity(H)
            diam, raio = max(exc.values()), min(exc.values())
            L = nx.average_shortest_path_length(H)
            pares = k * (k - 1) // 2
            soma_dist += L * pares
            pares_conectados += pares
        else:
            diam, raio, L = 0, 0, float("nan")
        linhas_comp.append({
            "rede": G.nodes[next(iter(c))]["rede"], "vertices": k,
            "arestas": H.number_of_edges(), "diametro": diam, "raio": raio,
            "caminho_medio": L, "densidade": nx.density(H) if k > 1 else float("nan"),
        })
    comp_df = pd.DataFrame(linhas_comp)
    comp_df.to_csv(DIR_PROC / "parte1_componentes.csv", index=False)
    nao_triv = comp_df[comp_df.vertices > 1]

    tri_por_no = nx.triangles(G)
    medidas = [
        ("Número de vértices", n),
        ("Número de arestas", m),
        ("Grau mínimo", min(graus)),
        ("Grau máximo", max(graus)),
        ("Grau médio", round(2 * m / n, 4)),
        ("Densidade (grafo inteiro)", round(nx.density(G), 6)),
        ("Densidade média por componente (>1 vértice)", round(nao_triv.densidade.mean(), 4)),
        ("Número de componentes conexas", len(comps)),
        ("Componentes não triviais (>1 vértice)", len(nao_triv)),
        ("Vértices isolados", nx.number_of_isolates(G)),
        ("Maior componente (vértices)", tamanhos[0]),
        ("Diâmetro (máx. entre componentes)", int(comp_df.diametro.max())),
        ("Diâmetro médio (componentes não triviais)", round(nao_triv.diametro.mean(), 4)),
        ("Raio (grafo inteiro, desconexo)", "∞"),
        ("Raio mínimo / máximo (componentes não triviais)",
         f"{int(nao_triv.raio.min())} / {int(nao_triv.raio.max())}"),
        ("Raio médio (componentes não triviais)", round(nao_triv.raio.mean(), 4)),
        ("Caminho médio (pares conectados)", round(soma_dist / pares_conectados, 4)),
        ("Caminho médio (média das componentes)", round(nao_triv.caminho_medio.mean(), 4)),
        ("Coef. de clusterização médio (todos os vértices)", round(nx.average_clustering(G), 4)),
        ("Coef. de clusterização médio (grau ≥ 2)",
         round(sum(nx.clustering(G, v) for v in G if G.degree(v) >= 2)
               / sum(1 for v in G if G.degree(v) >= 2), 4)),
        ("Transitividade global", round(nx.transitivity(G), 4)),
        ("Número de triângulos", sum(tri_por_no.values()) // 3),
    ]
    med_df = pd.DataFrame(medidas, columns=["medida", "valor"])
    med_df.to_csv(DIR_PROC / "parte1_medidas.csv", index=False)
    pd.DataFrame(sorted(dist_graus.items()), columns=["grau", "vertices"]) \
        .assign(P=lambda d: d.vertices / n) \
        .to_csv(DIR_PROC / "parte1_distribuicao_graus.csv", index=False)

    checklist = pd.DataFrame([
        ("Número de vértices", "Sim", ""),
        ("Número de arestas", "Sim", ""),
        ("Grau mínimo, máximo e médio", "Sim", ""),
        ("Distribuição de graus", "Sim", "tabela + figura"),
        ("Densidade", "Sim", "global e por componente"),
        ("Número de componentes conexas", "Sim", ""),
        ("Tamanho de cada componente", "Sim", "parte1_componentes.csv + figura"),
        ("Diâmetro", "Sim", "por componente; global = máximo"),
        ("Raio", "Sim", "infinito no grafo desconexo; reportado por componente"),
        ("Comprimento médio dos caminhos", "Sim", "sobre pares conectados"),
        ("Coeficiente de clusterização médio", "Sim", ""),
        ("Número de triângulos", "Sim", "viável: grafo pequeno"),
        ("Visualização do grafo", "Sim", "grafo completo (451 vértices)"),
    ], columns=["medida", "computada", "observacao"])
    checklist.to_csv(DIR_PROC / "parte1_checklist.csv", index=False)

    figura_graus(dist_graus, n)
    figura_componentes(tamanhos)
    figura_grafo(G, comps)

    print(med_df.to_string(index=False))
    print("\nDistribuição de graus:", dict(sorted(dist_graus.items())))
    print("Tamanhos das componentes (tamanho: quantidade):",
          dict(sorted(Counter(tamanhos).items())))
    print("\nComponentes não triviais por diâmetro:",
          dict(sorted(Counter(nao_triv.diametro).items())))
    print("Componentes não triviais por raio:", dict(sorted(Counter(nao_triv.raio).items())))
    print("\nArquivos: data/processed/parte1_{medidas,componentes,distribuicao_graus,checklist}.csv;"
          " figures/parte1_{distribuicao_graus,tamanho_componentes,grafo}.png")


if __name__ == "__main__":
    main()
