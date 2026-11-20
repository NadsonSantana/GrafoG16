"""Etapa 4 - Parte II: algoritmos da disciplina no G16 tratado.

1. Executa as implementacoes proprias (src/algoritmos.py) e valida cada
   resultado contra o NetworkX.
2. Tempo real no grafo completo: 30 amostras por algoritmo (cada amostra =
   tempo medio de um lote de execucoes que dura >= 5 ms)
   (n >= 30 -> IC 95% pela Normal, z).
3. Escalabilidade: subgrafos com as k primeiras partidas, 10 repeticoes por
   ponto (n < 30 e sigma desconhecido -> IC 95% pela t-Student), comparando o
   tempo observado com a funcao de custo teorica.
"""
import gc
import math
import time
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

import algoritmos as alg

RAIZ = Path(__file__).resolve().parents[1]
DIR_PROC = RAIZ / "data" / "processed"
DIR_FIG = RAIZ / "figures"

REPS_COMPLETO = 30
REPS_ESCALA = 10
PARTIDAS_ESCALA = [8, 16, 24, 32, 40, 48, 56, 62]
CONFIANCA = 0.95
LOTE_MIN_S = 0.005

AZUL = "#2a78d6"
AZUL_CLARO = "#9ec5f4"
TEXTO = "#0b0b0b"
TEXTO_2 = "#52514e"
GRADE = "#e5e4e0"


# ------------------------------------------------------------- dados
def carregar():
    nos = pd.read_csv(DIR_PROC / "g16_nos.csv")
    und = pd.read_csv(DIR_PROC / "g16_nao_dirigido.csv")
    dir_ = pd.read_csv(DIR_PROC / "g16_dirigido.csv")
    return nos, und, dir_


def montar(nos, und, dir_, partidas=None):
    """Estruturas de entrada dos algoritmos (opcionalmente so as k primeiras partidas)."""
    if partidas is not None:
        manter = set(nos[nos.rede < partidas].no)
        nos = nos[nos.no.isin(manter)]
        und = und[und.u.isin(manter)]
        dir_ = dir_[dir_.u.isin(manter)]
    V = list(nos.no)
    adj = {v: {} for v in V}
    arestas = []
    for u, v, w in zip(und.u, und.v, und.distancia):
        adj[u][v] = w
        adj[v][u] = w
        arestas.append((u, v, w))
    arcos = [(u, v, w) for u, v, w in arestas] + [(v, u, w) for u, v, w in arestas]
    adjD = {v: {} for v in V}
    for u, v, a in zip(dir_.u, dir_.v, dir_.atencao):
        adjD[u][v] = 1 - a
    # Origem: vertice de maior grau (desempate pelo rotulo)
    fonte = max(V, key=lambda x: (len(adj[x]), x))
    return {"V": V, "adj": adj, "arestas": arestas, "arcos": arcos, "adjD": adjD,
            "fonte": fonte, "n": len(V), "m": len(arestas), "mD": len(dir_)}


def para_nx(d):
    G = nx.Graph()
    G.add_nodes_from(d["V"])
    G.add_weighted_edges_from(d["arestas"], weight="distancia")
    D = nx.DiGraph()
    D.add_nodes_from(d["V"])
    D.add_weighted_edges_from([(u, v, w) for u in d["adjD"] for v, w in d["adjD"][u].items()],
                              weight="distancia")
    return G, D


# ------------------------------------------------------------- algoritmos
ALGORITMOS = [
    # nome, chamada, complexidade teorica, funcao de custo (n, m, mD)
    ("BFS", lambda d: alg.bfs_completa(d["adj"]), "O(V + E)",
     lambda n, m, mD: n + m),
    ("DFS", lambda d: alg.dfs_completa(d["adj"]), "O(V + E)",
     lambda n, m, mD: n + m),
    ("Eulerianidade", lambda d: alg.eulerianidade(d["adj"]), "O(V + E)",
     lambda n, m, mD: n + m),
    ("Dijkstra", lambda d: alg.dijkstra(d["adj"], d["fonte"]), "O((V + E) log V)",
     lambda n, m, mD: (n + m) * math.log2(n)),
    ("Bellman-Ford", lambda d: alg.bellman_ford(d["V"], d["arcos"], d["fonte"]), "O(V · E)",
     lambda n, m, mD: n * 2 * m),
    ("Floyd-Warshall", lambda d: alg.floyd_warshall(d["V"], d["adj"]), "O(V³)",
     lambda n, m, mD: n ** 3),
    ("Tarjan (CFC, dirigido)", lambda d: alg.tarjan_scc(d["adjD"]), "O(V + E)",
     lambda n, m, mD: n + mD),
    ("Tarjan (articulações/pontes)", lambda d: alg.tarjan_articulacoes_pontes(d["adj"]),
     "O(V + E)", lambda n, m, mD: n + m),
    ("Kruskal", lambda d: alg.kruskal(d["V"], d["arestas"]), "O(E log E)",
     lambda n, m, mD: m * math.log2(m)),
    ("Prim", lambda d: alg.prim(d["adj"]), "O(E log V)",
     lambda n, m, mD: m * math.log2(n)),
]


# ------------------------------------------------------------- estatistica
def medir(fn, reps, lote_min=LOTE_MIN_S):
    """Tempo de parede por execucao; coletor de lixo desligado na medicao.
    Como no timeit, cada amostra executa fn em lote (k vezes, k calibrado para
    o lote durar >= lote_min) e registra o tempo medio por execucao, reduzindo
    o efeito da resolucao do relogio e de interrupcoes do sistema em
    algoritmos muito rapidos. Retorna (amostras em s, k)."""
    t0 = time.perf_counter()
    fn()  # aquecimento (nao contabilizado) e calibracao
    um = time.perf_counter() - t0
    k = max(1, math.ceil(lote_min / max(um, 1e-9)))
    t = []
    gc.collect()
    gc.disable()
    try:
        for _ in range(reps):
            t0 = time.perf_counter()
            for _ in range(k):
                fn()
            t.append((time.perf_counter() - t0) / k)
    finally:
        gc.enable()
    return np.array(t), k


def intervalo(amostra, confianca=CONFIANCA):
    """IC da media: z se n >= 30; t-Student se n < 30 (sigma desconhecido)."""
    n = len(amostra)
    media = amostra.mean()
    s = amostra.std(ddof=1)
    if n >= 30:
        crit, dist = stats.norm.ppf(1 - (1 - confianca) / 2), "Normal (z)"
    else:
        crit, dist = stats.t.ppf(1 - (1 - confianca) / 2, n - 1), f"t-Student (gl={n - 1})"
    meia = crit * s / math.sqrt(n)
    return media, s, media - meia, media + meia, crit, dist


# ------------------------------------------------------------- validacao
def validar(d):
    G, D = para_nx(d)
    adj, fonte = d["adj"], d["fonte"]
    res, saidas = [], {}

    comp_ref = sorted(sorted(c) for c in nx.connected_components(G))
    bfs_c = alg.bfs_completa(adj)
    dfs_c = alg.dfs_completa(adj)
    res.append(("BFS", sorted(sorted(c) for c in bfs_c) == comp_ref))
    res.append(("DFS", sorted(sorted(c) for c in dfs_c) == comp_ref))
    _, dist_bfs, _ = alg.bfs(adj, fonte)
    res.append(("BFS (distâncias em saltos)",
                dist_bfs == dict(nx.single_source_shortest_path_length(G, fonte))))
    saidas["BFS"] = f"{len(bfs_c)} árvores de busca; a partir de {fonte}: {len(dist_bfs)} vértices alcançados, profundidade máx. {max(dist_bfs.values())}"
    saidas["DFS"] = f"{len(dfs_c)} árvores de busca (igual ao nº de componentes)"

    eu = alg.eulerianidade(adj)
    res.append(("Eulerianidade (grafo inteiro)",
                (eu["classe"] == "euleriano") == nx.is_eulerian(G)))
    classes = Counter()
    ok_comp = True
    for c in nx.connected_components(G):
        if len(c) == 1:
            continue
        sub = {u: {v: w for v, w in adj[u].items()} for u in c}
        cl = alg.eulerianidade(sub)["classe"]
        classes[cl] += 1
        H = G.subgraph(c)
        ok_comp &= (cl == "euleriano") == nx.is_eulerian(H)
        ok_comp &= (cl in ("euleriano", "semi-euleriano")) == nx.has_eulerian_path(H)
    res.append(("Eulerianidade (por componente)", ok_comp))
    saidas["Eulerianidade"] = (f"grafo inteiro: {eu['classe']} ({eu['impares']} vértices de grau ímpar; "
                               f"arestas em uma única componente: {eu['conexo_nas_arestas']}); "
                               f"componentes não triviais: {dict(classes)}")

    dist_dj, pai_dj = alg.dijkstra(adj, fonte)
    ref = nx.single_source_dijkstra_path_length(G, fonte, weight="distancia")
    res.append(("Dijkstra", dist_dj.keys() == ref.keys()
                and all(math.isclose(dist_dj[k], ref[k]) for k in ref)))
    alvo = max(dist_dj, key=dist_dj.get)
    caminho, x = [], alvo
    while x is not None:
        caminho.append(x)
        x = pai_dj[x]
    saidas["Dijkstra"] = (f"origem {fonte}: {len(dist_dj)} alcançados; mais distante {alvo} "
                          f"(d = {dist_dj[alvo]:.4f}) via {' → '.join(reversed(caminho))}")

    dist_bf, _, neg = alg.bellman_ford(d["V"], d["arcos"], fonte)
    fin = {k: v for k, v in dist_bf.items() if v != alg.INF}
    ref_bf = nx.single_source_bellman_ford_path_length(G, fonte, weight="distancia")
    res.append(("Bellman-Ford", fin.keys() == ref_bf.keys()
                and all(math.isclose(fin[k], ref_bf[k]) for k in ref_bf) and not neg))
    saidas["Bellman-Ford"] = (f"mesmas distâncias do Dijkstra: "
                              f"{all(math.isclose(fin[k], dist_dj[k]) for k in dist_dj)}; "
                              f"ciclo negativo: {neg}")

    idx, M = alg.floyd_warshall(d["V"], adj)
    ref_fw = dict(nx.all_pairs_dijkstra_path_length(G, weight="distancia"))
    ok = True
    for u in d["V"]:
        linha = M[idx[u]]
        finitos = {d["V"][j] for j in np.flatnonzero(np.isfinite(linha))}
        ok &= finitos == ref_fw[u].keys()
        ok &= all(math.isclose(linha[idx[v]], ref_fw[u][v], abs_tol=1e-12) for v in ref_fw[u])
    res.append(("Floyd-Warshall", ok))
    finitos = M[np.isfinite(M) & (M > 0)]
    saidas["Floyd-Warshall"] = (f"matriz {M.shape[0]}×{M.shape[1]}; {finitos.size // 2} pares conectados; "
                                f"distância ponderada média {finitos.mean():.4f}, máxima {finitos.max():.4f}")

    sccs = alg.tarjan_scc(d["adjD"])
    res.append(("Tarjan (CFC)", sorted(sorted(c) for c in sccs)
                == sorted(sorted(c) for c in nx.strongly_connected_components(D))))
    tam = Counter(len(c) for c in sccs)
    saidas["Tarjan (CFC, dirigido)"] = (f"{len(sccs)} componentes fortemente conexas; "
                                        f"tamanhos {dict(sorted(tam.items()))}")

    art, pontes = alg.tarjan_articulacoes_pontes(adj)
    res.append(("Tarjan (articulações)", art == set(nx.articulation_points(G))))
    res.append(("Tarjan (pontes)", {frozenset(p) for p in pontes}
                == {frozenset(p) for p in nx.bridges(G)}))
    saidas["Tarjan (articulações/pontes)"] = f"{len(art)} pontos de articulação; {len(pontes)} pontes"

    kr = alg.kruskal(d["V"], d["arestas"])
    pr = alg.prim(adj)
    ref_mst = nx.minimum_spanning_tree(G, weight="distancia")
    peso_ref = ref_mst.size(weight="distancia")
    res.append(("Kruskal", math.isclose(sum(w for *_, w in kr), peso_ref)
                and len(kr) == ref_mst.number_of_edges()))
    res.append(("Prim", math.isclose(sum(w for *_, w in pr), peso_ref)
                and len(pr) == ref_mst.number_of_edges()))
    saidas["Kruskal"] = (f"floresta com {len(kr)} arestas (= V − nº componentes = "
                         f"{d['n']} − {nx.number_connected_components(G)}); peso total {peso_ref:.4f}")
    saidas["Prim"] = f"mesma floresta em peso: {math.isclose(sum(w for *_, w in pr), sum(w for *_, w in kr))}"
    return res, saidas


# ------------------------------------------------------------- figuras
def estilo(ax, grade="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(TEXTO_2)
    ax.tick_params(colors=TEXTO_2, labelsize=8)
    ax.grid(axis=grade, color=GRADE, linewidth=0.8)
    ax.set_axisbelow(True)


def figura_tempos(df):
    df = df.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
    y = np.arange(len(df))
    media = df["media_ms"].to_numpy()
    erro = np.vstack([media - df["ic_inf_ms"].to_numpy(), df["ic_sup_ms"].to_numpy() - media])
    ax.barh(y, media, color=AZUL, height=0.6, edgecolor="white")
    ax.errorbar(media, y, xerr=erro, fmt="none", ecolor=TEXTO, elinewidth=1, capsize=3)
    for yi, mi in zip(y, media):
        ax.text(mi * 1.6, yi, f"{mi:.3g} ms", va="center", fontsize=8, color=TEXTO_2)
    ax.set_yticks(y)
    ax.set_yticklabels(df["algoritmo"], color=TEXTO)
    ax.set_xscale("log")
    ax.set_xlim(right=ax.get_xlim()[1] * 4)
    ax.set_xlabel("Tempo médio por execução (ms, escala log) — barra = IC 95% (z, n = 30)",
                  color=TEXTO)
    ax.set_title("Tempo real no G16 tratado (V = 451, E = 548)", color=TEXTO,
                 fontsize=11, loc="left")
    estilo(ax, grade="x")
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte2_tempos.png")
    plt.close(fig)


def figura_escala(esc):
    nomes = list(dict.fromkeys(esc.algoritmo))
    cols = 5
    lin = math.ceil(len(nomes) / cols)
    fig, axes = plt.subplots(lin, cols, figsize=(15, 3.2 * lin), dpi=150)
    for ax, nome in zip(axes.flat, nomes):
        s = esc[esc.algoritmo == nome]
        x = s.custo.to_numpy()
        y = s.media_ms.to_numpy()
        ax.errorbar(x, y, yerr=[y - s.ic_inf_ms, s.ic_sup_ms - y], fmt="o", color=AZUL,
                    ms=4, elinewidth=1, capsize=2)
        a = (x @ y) / (x @ x)  # reta pela origem: t = a * custo
        r2 = 1 - ((y - a * x) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        xs = np.linspace(0, x.max() * 1.05, 50)
        ax.plot(xs, a * xs, color=AZUL_CLARO, linewidth=2, zorder=0)
        ajuste = f"R² = {r2:.3f}" if r2 >= 0 else "tempo ~constante"
        ax.set_title(f"{nome}\n{s.complexidade.iloc[0]}  ·  {ajuste}", fontsize=9,
                     color=TEXTO, loc="left")
        ax.set_xlabel("custo teórico f(V, E)", fontsize=8, color=TEXTO_2)
        ax.set_ylabel("ms", fontsize=8, color=TEXTO_2)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        estilo(ax)
    for ax in list(axes.flat)[len(nomes):]:
        ax.axis("off")
    fig.suptitle("Escalabilidade: tempo observado × custo teórico (subgrafos com 8…62 partidas; "
                 "IC 95% t-Student, n = 10; reta pela origem)", fontsize=11, color=TEXTO, x=0.01,
                 ha="left")
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte2_escalabilidade.png")
    plt.close(fig)


# ------------------------------------------------------------- main
def main():
    nos, und, dir_ = carregar()
    d = montar(nos, und, dir_)
    print(f"G16: V = {d['n']}, E = {d['m']} (dirigido: {d['mD']} arcos); origem = {d['fonte']}\n")

    validacoes, saidas = validar(d)
    print("Validação contra NetworkX:")
    for nome, ok in validacoes:
        print(f"  {'OK ' if ok else 'FALHOU'} {nome}")
    pd.DataFrame(validacoes, columns=["teste", "ok"]).to_csv(DIR_PROC / "parte2_validacao.csv",
                                                             index=False)
    print("\nResultados dos algoritmos:")
    for nome, txt in saidas.items():
        print(f"  {nome}: {txt}")
    pd.DataFrame(list(saidas.items()), columns=["algoritmo", "resultado"]) \
        .to_csv(DIR_PROC / "parte2_resultados.csv", index=False)

    linhas = []
    for nome, fn, comp, _ in ALGORITMOS:
        t, lote = medir(lambda: fn(d), REPS_COMPLETO)
        t = t * 1000
        media, s, inf, sup, crit, dist = intervalo(t)
        linhas.append({"algoritmo": nome, "complexidade": comp, "n": len(t), "lote": lote,
                       "media_ms": media,
                       "desvio_ms": s, "ic_inf_ms": inf, "ic_sup_ms": sup,
                       "valor_critico": crit, "distribuicao": dist})
    tempos = pd.DataFrame(linhas)
    tempos.to_csv(DIR_PROC / "parte2_tempos.csv", index=False)
    print(f"\nTempo real no grafo completo ({REPS_COMPLETO} repetições, IC {CONFIANCA:.0%}):")
    print(tempos.drop(columns=["valor_critico"]).to_string(
        index=False, float_format=lambda v: f"{v:.4f}"))
    figura_tempos(tempos)

    linhas = []
    for k in PARTIDAS_ESCALA:
        dk = montar(nos, und, dir_, partidas=k)
        for nome, fn, comp, custo in ALGORITMOS:
            t, lote = medir(lambda: fn(dk), REPS_ESCALA)
            t = t * 1000
            media, s, inf, sup, crit, dist = intervalo(t)
            linhas.append({"algoritmo": nome, "complexidade": comp, "partidas": k, "lote": lote,
                           "V": dk["n"],
                           "E": dk["m"], "custo": custo(dk["n"], dk["m"], dk["mD"]),
                           "n": len(t), "media_ms": media, "desvio_ms": s, "ic_inf_ms": inf,
                           "ic_sup_ms": sup, "distribuicao": dist})
    esc = pd.DataFrame(linhas)
    esc.to_csv(DIR_PROC / "parte2_escalabilidade.csv", index=False)
    figura_escala(esc)
    print(f"\nEscalabilidade ({REPS_ESCALA} repetições por ponto, {esc.distribuicao.iloc[0]}):")
    print(esc.pivot(index="partidas", columns="algoritmo", values="media_ms")
          .round(4).to_string())
    print("\nArquivos: data/processed/parte2_{validacao,resultados,tempos,escalabilidade}.csv; "
          "figures/parte2_{tempos,escalabilidade}.png")


if __name__ == "__main__":
    main()
