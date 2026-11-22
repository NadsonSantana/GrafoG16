"""Parte II: algoritmos da disciplina no grafo completo (tempo real com IC 95%) e
experimento de escala da complexidade observada (ajuste log-log sobre a mediana
de cada ponto, robusta a execuções atípicas; média, desvio e IC 95% também são salvos).

Uso: python -m src.parte2 [--rapido | --so-figuras]"""

import math
import random

import numpy as np

from . import algoritmos as alg
from .benchmark import medir
from .comum import (AQUA, AZUL, LARANJA, TINTA, TINTA_2, VERMELHO, carregar, matplotlib,
                    ler_json, plt, preparar, salvar_figura, salvar_json, so_figuras)


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
        ("Floyd-Warshall (numpy)", "O(V^3)", lambda: alg.floyd_warshall_numpy(g), n_fw),
        ("Tarjan", "O(V+E)", lambda: alg.tarjan(d), n_rapido),
        ("Tarjan (pontes e articulações)", "O(V+E)", lambda: alg.pontes_articulacoes(g), n_rapido),
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
    vs_np, fw_np = saidas["Floyd-Warshall (numpy)"]
    fw_numpy_concorda = vs_np == vs and bool(np.allclose(np.array(fw), fw_np))
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
    pontes, articulacoes = saidas["Tarjan (pontes e articulações)"]
    # exemplo de caminho mínimo: maior distância finita encontrada
    i_max, j_max = max(((i, j) for i in range(len(vs)) for j in range(len(vs))
                        if i != j and fw[i][j] < math.inf), key=lambda p: fw[p[0]][p[1]])
    resumo = {
        "bfs_visitados": len(saidas["BFS"]),
        "dfs_visitados": len(saidas["DFS"]),
        "euler_global": saidas["Eulerianidade"],
        "euler_por_partida": {t: euler_por_partida.count(t) for t in sorted(set(euler_por_partida))},
        "pares_alcancaveis": len(finitos),
        "distancia_media": float(np.mean(finitos)),
        "distancia_max": float(max(finitos)),
        "par_distancia_max": [vs[i_max], vs[j_max]],
        "dijkstra_bf_fw_concordam": concordam,
        "fw_numpy_concorda": fw_numpy_concorda,
        "bellman_ford_ciclo_negativo": False,
        "num_scc": len(sccs),
        "scc_tamanhos": {t: [len(c) for c in sccs].count(t)
                         for t in sorted({len(c) for c in sccs})},
        "scc_nao_triviais": sum(1 for c in sccs if len(c) > 1),
        "partidas_fortemente_conexas": sum(
            1 for c in sccs if len(c) > 1 and len(c) == sum(
                1 for v in g.adj if v.split("_")[0] == c[0].split("_")[0])),
        "pontes": len(pontes),
        "articulacoes": len(articulacoes),
        "partidas_com_ponte": len({u.split("_")[0] for u, _ in pontes}),
        "partidas_com_articulacao": len({v.split("_")[0] for v in articulacoes}),
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
    salvar_figura(fig, "tempos")


TEORICA = {  # nome -> (rótulo LaTeX, f(V, E), usa versão dirigida)
    "BFS": ("O(V+E)", lambda v, e: v + e, False),
    "DFS": ("O(V+E)", lambda v, e: v + e, False),
    "Eulerianidade": ("O(V+E)", lambda v, e: v + e, False),
    "Dijkstra*": ("O(V(V+E)\\log V)", lambda v, e: v * (v + e) * math.log2(v), False),
    "Bellman-Ford*": ("O(V^2E)", lambda v, e: v * v * e, False),
    "Floyd-Warshall": ("O(V^3)", lambda v, e: v ** 3, False),
    "Floyd-Warshall (numpy)": ("O(V^3)", lambda v, e: v ** 3, False),
    "Tarjan": ("O(V+E)", lambda v, e: v + e, True),
    "Pontes/articulações": ("O(V+E)", lambda v, e: v + e, False),
    "Prim": ("O(E\\log V)", lambda v, e: e * math.log2(v), False),
    "Kruskal": ("O(E\\log V)", lambda v, e: e * math.log2(v), False),
}
FUNCOES = {
    "BFS": alg.bfs_completa, "DFS": alg.dfs_completa, "Eulerianidade": alg.verificar_euleriano,
    "Dijkstra*": alg.dijkstra_todos, "Bellman-Ford*": alg.bellman_ford_todos,
    "Floyd-Warshall": alg.floyd_warshall, "Floyd-Warshall (numpy)": alg.floyd_warshall_numpy,
    "Tarjan": alg.tarjan,
    "Pontes/articulações": alg.pontes_articulacoes,
    "Prim": alg.prim, "Kruskal": alg.kruskal,
}


# número de partidas de cada ponto do experimento de escala: até 62 são as primeiras
# partidas reais; acima disso, as 62 reais mais cópias de partidas sorteadas com reposição
PARTIDAS_ESCALA = [8, 16, 31, 62, 124, 248, 496, 992, 1984, 3968]
LIMITE_PARTIDAS = {"Bellman-Ford*": 124, "Floyd-Warshall": 248, "Floyd-Warshall (numpy)": 124}
REPETICOES_REDUZIDAS = ("Bellman-Ford*", "Floyd-Warshall", "Floyd-Warshall (numpy)")


def replicar_partidas(g, d, quantidade, rng):
    """União de 'quantidade' partidas, nas versões não-direcionada e dirigida.

    Até 62, usa as primeiras partidas reais. Acima disso, acrescenta às 62 reais cópias de
    partidas sorteadas com reposição, com vértices renomeados (cXXXXX_pY). Cada cópia
    preserva a estrutura e os pesos de uma partida real, de modo que V e E crescem
    mantendo o tipo de grafo analisado."""
    partidas = {}
    for v in g.adj:
        partidas.setdefault(v.split("_")[0], []).append(v)
    nomes = sorted(partidas)
    if quantidade <= len(nomes):
        vs = [v for p in nomes[:quantidade] for v in partidas[p]]
        return g.subgrafo(vs), d.subgrafo(vs)
    h, hd = g.subgrafo(list(g.adj)), d.subgrafo(list(d.adj))
    for j in range(quantidade - len(nomes)):
        p = rng.choice(nomes)
        novo = {v: f"c{j:05d}_{v.split('_', 1)[1]}" for v in partidas[p]}
        for u in partidas[p]:
            h.adj[novo[u]] = {novo[v]: w for v, w in g.adj[u].items()}
            hd.adj[novo[u]] = {novo[v]: w for v, w in d.adj[u].items()}
    return h, hd


def parte2_escala(g, d, rapido=False, semente=42):
    """Complexidade observada: cada algoritmo é medido em grafos com 8 a 3.968 partidas
    (V de ~60 a ~29 mil), e ajusta-se log(tempo) ~ b log f(V, E) e log(tempo) ~ a log V."""
    rng = random.Random(semente)
    tamanhos = [16, 62, 124] if rapido else PARTIDAS_ESCALA
    linhas = []
    for k in tamanhos:
        sub_g, sub_d = replicar_partidas(g, d, k, rng)
        for nome, (_, f, dirigido) in TEORICA.items():
            if k > LIMITE_PARTIDAS.get(nome, math.inf):
                continue
            h = sub_d if dirigido else sub_g
            n = 3 if rapido else (10 if nome in REPETICOES_REDUZIDAS else 30)
            _, est = medir(lambda: FUNCOES[nome](h), n)
            v_, e_ = h.num_vertices(), h.num_arestas()
            linhas.append({"algoritmo": nome, "partidas": k, "real": k <= 62, "V": v_, "E": e_,
                           "f_teorica": f(v_, e_), **est})
        print(f"  escala: {k} partidas ({sub_g.num_vertices()} vértices)")
    ajuste = {}
    for nome, (rotulo, _, _) in TEORICA.items():
        pts = [x for x in linhas if x["algoritmo"] == nome]
        x = np.log([p["f_teorica"] for p in pts])
        y = np.log([p["mediana"] for p in pts])  # mediana: robusta a execuções atípicas
        b, _ = np.polyfit(x, y, 1)
        xv = np.log([p["V"] for p in pts])
        a, _ = np.polyfit(xv, y, 1)
        # curva teórica c * f(V, E) (inclinação fixa em 1 na escala log)
        log_c = float(np.mean(y - x))
        r2_c = 1 - float(np.sum((y - x - log_c) ** 2) / np.sum((y - y.mean()) ** 2))
        ajuste[nome] = {"complexidade": rotulo, "expoente": float(b),
                        "r2": float(np.corrcoef(x, y)[0, 1] ** 2),
                        "expoente_v": float(a), "r2_v": float(np.corrcoef(xv, y)[0, 1] ** 2),
                        "log_c": log_c, "r2_teorica": r2_c,
                        "v_min": min(p["V"] for p in pts), "v_max": max(p["V"] for p in pts)}
    salvar_json("parte2_escala.json", {"linhas": linhas, "ajuste": ajuste})
    figura_escala(linhas, ajuste)
    return ajuste


def figura_escala(linhas, ajuste):
    """Um painel por algoritmo: tempo médio (IC 95%) contra V em log-log e a curva
    teórica c * f(V, E) ajustada. Marcadores cheios: partidas reais; vazios: réplicas."""
    fig, eixos = plt.subplots(3, 4, figsize=(7.0, 5.4))
    eixos = eixos.ravel()
    for ax, nome in zip(eixos, TEORICA):
        pts = sorted((x for x in linhas if x["algoritmo"] == nome), key=lambda p: p["V"])
        vs = np.array([p["V"] for p in pts])
        medias = np.array([p["media"] * 1000 for p in pts])
        inf = np.array([max(p["ic_inf"], p["media"] * 1e-3) * 1000 for p in pts])
        sup = np.array([p["ic_sup"] * 1000 for p in pts])
        ax.errorbar(vs, medias, yerr=[medias - inf, sup - medias], fmt="none", ecolor=AZUL,
                    elinewidth=0.8, capsize=1.5)
        reais = np.array([p["real"] for p in pts])
        ax.plot(vs[reais], medias[reais], "o", color=AZUL, markersize=3.5)
        ax.plot(vs[~reais], medias[~reais], "o", markerfacecolor="white", color=AZUL,
                markersize=3.5)
        a = ajuste[nome]
        teorica = np.exp(a["log_c"]) * np.array([p["f_teorica"] for p in pts]) * 1000
        ax.plot(vs, teorica, "--", color=LARANJA, linewidth=1.2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{nome}\n$f={a['complexidade']}$, $b={a['expoente']:.2f}$".replace(".", ","),
                     fontsize=6.5)
        ax.tick_params(labelsize=6)
        ax.grid(which="minor", visible=False)
    for ax in eixos[8:11]:
        ax.set_xlabel("$V$ (log)", fontsize=7)
    for ax in eixos[::4]:
        ax.set_ylabel("tempo (ms, log)", fontsize=7)
    legenda = eixos[-1]
    legenda.axis("off")
    legenda.legend(handles=[
        plt.Line2D([], [], color=AZUL, marker="o", linestyle="", markersize=4,
                   label="partidas reais"),
        plt.Line2D([], [], color=AZUL, marker="o", markerfacecolor="white", linestyle="",
                   markersize=4, label="réplicas das partidas"),
        plt.Line2D([], [], color=AZUL, marker="|", linestyle="", markersize=8, label="IC 95%"),
        plt.Line2D([], [], color=LARANJA, linestyle="--", label="$c \\cdot f(V, E)$ ajustada"),
    ], loc="center", fontsize=7)
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    salvar_figura(fig, "escala")


def figuras():
    """Redesenha as figuras da Parte II a partir dos JSONs de results/."""
    figura_tempos(ler_json("parte2_tempos.json"))
    escala = ler_json("parte2_escala.json")
    figura_escala(escala["linhas"], escala["ajuste"])


def main():
    rapido = preparar()
    if so_figuras():
        print("Parte II - redesenhando figuras")
        figuras()
        return
    g, d = carregar()
    print("Parte II - algoritmos")
    parte2(g, d, rapido)
    parte2_escala(g, d, rapido)


if __name__ == "__main__":
    main()
