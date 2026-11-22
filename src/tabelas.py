"""Gera as tabelas do relatório (results/tab_*.tex) e macros com valores citados no
texto (results/valores.tex) a partir dos JSONs de results/.

Uso: python -m src.tabelas   (também é chamado ao final de python -m src.executar)
"""

import json
import os
import statistics

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTADOS = os.path.join(RAIZ, "results")
DADOS = os.path.join(RAIZ, "data", "processed")


def ler(nome, pasta=RESULTADOS):
    with open(os.path.join(pasta, nome)) as f:
        return json.load(f)


def num(x, casas=4):
    """Formata no padrão brasileiro: 1.234,5678."""
    if isinstance(x, int) or (isinstance(x, float) and x.is_integer() and casas == 0):
        return f"{int(x):,}".replace(",", ".")
    texto = f"{x:,.{casas}f}"
    return texto.replace(",", "X").replace(".", ",").replace("X", ".")


def cientifico(x, casas=2):
    mantissa, expoente = f"{x:.{casas}e}".split("e")
    return f"${mantissa.replace('.', '{,}')} \\times 10^{{{int(expoente)}}}$"


def tabular(colunas, cabecalho, linhas, rodape=None):
    partes = [f"\\begin{{tabular}}{{{colunas}}}", "  \\toprule",
              "  " + " & ".join(cabecalho) + " \\\\", "  \\midrule"]
    partes += ["  " + " & ".join(str(c) for c in linha) + " \\\\" for linha in linhas]
    partes.append("  \\bottomrule")
    if rodape:
        partes.append(f"  \\multicolumn{{{len(cabecalho)}}}{{l}}{{\\footnotesize {rodape}}} \\\\")
    partes.append("\\end{tabular}")
    return "\n".join(partes) + "\n"


def salvar(nome, conteudo):
    with open(os.path.join(RESULTADOS, nome), "w") as f:
        f.write(conteudo)


# ---------------------------------------------------------------- tabelas

def tab_etapas():
    e = ler("etapas_tratamento.json", DADOS)
    v = num(e["vertices"])
    linhas = [
        ("Grafo temporal bruto (com laptop)", v,
         num(e["arestas_temporais"] + e["arestas_temporais_laptop"])),
        ("Sem laptop", v, num(e["arestas_temporais"])),
        ("Estático (pares que ocorrem)", v, num(e["pares_estaticos"]) + " arcos"),
        ("Após limiar $f \\geq 1/N$", v, num(e["pares_apos_limiar"]) + " arcos"),
        ("Não-direcionado (final)", v, num(e["arestas_nao_dirigidas"])),
    ]
    salvar("tab_etapas.tex", tabular("lrr", ["Etapa", "Vértices", "Arestas"], linhas))


def tab_metricas():
    m = ler("parte1_metricas.json")
    sw = ler("parte3_small_world.json")["por_partida"]
    tam = m["tamanhos_componentes_freq"]
    exato = "exato"
    linhas = [
        ("Número de vértices", num(m["vertices"]), exato),
        ("Número de arestas", num(m["arestas"]), exato),
        ("Grau mín. / máx.", f"{m['grau_min']} / {m['grau_max']}", exato),
        ("Grau médio / mediano", f"{num(m['grau_medio'])} / {num(m['grau_mediano'], 0)}", exato),
        ("Distribuição de graus", "Fig.~\\ref{fig:graus}", "exato (todos os vértices)"),
        ("Densidade global", cientifico(m["densidade"]), exato),
        ("Densidade média intra-partida", num(sw["densidade_media"]), "média das 62 partidas"),
        ("Número de componentes", num(m["num_componentes"]), "BFS"),
        ("Tamanho das componentes", f"{min(map(int, tam))} a {max(map(int, tam))}",
         "BFS (distribuição no texto)"),
        ("Tamanho da LCC", num(m["tamanho_lcc"]), "BFS"),
        ("Diâmetro (LCC / máx.)", f"{m['diametro_lcc']} / {m['diametro_max_componentes']}",
         "BFS de todos os vértices"),
        ("Raio (LCC / mín.)", f"{m['raio_lcc']} / {m['raio_min_componentes']}",
         "BFS de todos os vértices"),
        ("Comp. médio (LCC / global)",
         f"{num(m['caminho_medio_lcc'])} / {num(m['caminho_medio_global'])}",
         "BFS de todos os vértices"),
        ("Clusterização média", num(m["clusterizacao_media"]), "exato (vizinhanças)"),
        ("Número de triângulos", num(m["triangulos"]), "exato (vizinhanças)"),
        ("Visualização", "Fig.~\\ref{fig:teaser}, \\ref{fig:partida0}", "grafo completo"),
    ]
    salvar("tab_metricas.tex", tabular("llp{0.30\\columnwidth}",
                                       ["Métrica", "Valor", "Método"], linhas))


def tab_tempos():
    t = ler("parte2_tempos.json")
    linhas = []
    for nome, e in t.items():
        ic = f"[{num(e['ic_inf'] * 1000, 3)}; {num(e['ic_sup'] * 1000, 3)}]"
        if e["distribuicao"] == "t":
            ic += "*"
        linhas.append((nome, f"${e['complexidade']}$", e["n"], num(e["media"] * 1000, 3),
                       num(e["desvio"] * 1000, 3), ic))
    n_t = [e["n"] for e in t.values() if e["distribuicao"] == "t"]
    rodape = (f"* IC com distribuição t-Student ($n = {n_t[0]}$, {n_t[0] - 1} graus de liberdade); "
              "os demais, com distribuição Normal ($n \\geq 30$).") if n_t else None
    salvar("tab_tempos.tex", tabular(
        "llrrrc", ["Algoritmo", "Complexidade", "$n$", "Média (ms)", "Desvio (ms)", "IC 95\\% (ms)"],
        linhas, rodape))


def tab_escala():
    a = ler("parte2_escala.json")["ajuste"]
    linhas = [(nome, f"${x['complexidade']}$", num(x["expoente"], 2), num(x["r2"], 3),
               num(x["expoente_v"], 2), num(x["r2_v"], 3)) for nome, x in a.items()]
    salvar("tab_escala.tex", tabular(
        "llrrrr", ["Algoritmo", "$f(V,E)$", "$b$", "$R^2$", "$a$", "$R^2$"], linhas,
        "* a partir de todas as fontes. $b$: tempo $\\propto f(V,E)^b$; $a$: tempo $\\propto V^a$."))


def tab_small_world():
    sw = ler("parte3_small_world.json")["global"]
    g, b, t = sw["er_global"], sw["er_blocos"], sw["teorico"]
    linhas = [
        ("Clusterização $C$", num(sw["C"]), num(g["C_rand"]), num(b["C_rand"]), num(t["C_rand"])),
        ("Comp. médio $L$", num(sw["L"]), num(g["L_rand"]), num(b["L_rand"]), num(t["L_rand"])),
        ("Componentes", "62", num(g["componentes_rand"], 1), num(b["componentes_rand"], 1), "--"),
        ("$C/C_{rand}$", "--", num(g["gamma"], 1), num(b["gamma"], 2),
         num(sw["C"] / t["C_rand"], 1)),
        ("$L/L_{rand}$", "--", num(g["lambda"], 2), num(b["lambda"], 2),
         num(sw["L"] / t["L_rand"], 2)),
        ("$\\sigma$", "--", f"\\textbf{{{num(g['sigma'], 1)}}}", f"\\textbf{{{num(b['sigma'], 2)}}}",
         num((sw["C"] / t["C_rand"]) / (sw["L"] / t["L_rand"]), 1)),
    ]
    salvar("tab_smallworld.tex", tabular(
        "lrrrr", ["Métrica", "Observado", "$G(n,m)$", "Por blocos", "Teórico"], linhas))


def tab_potencia():
    lp = ler("parte3_lei_potencia.json")
    mv = lp["mv"]
    linhas = [
        ("Regressão log-log de $P(k)$: inclinação", num(-lp["gamma_mq"], 2)),
        ("Regressão log-log: $R^2$", num(lp["r2"], 2)),
        ("MV: $k_{min}$ (mínima distância KS)", mv["kmin"]),
        ("MV: $\\gamma \\pm$ erro padrão", f"{num(mv['gamma'], 2)} $\\pm$ {num(mv['erro_padrao'], 2)}"),
        ("MV: vértices na cauda",
         f"{mv['n_cauda']} ({sum(1 for k in lp['ks'] if k >= mv['kmin'])} valores de $k$)"),
        ("MV: distância KS", num(mv["ks"], 3)),
        ("$\\langle k \\rangle$ / Var($k$)", f"{num(lp['k_medio'], 2)} / {num(lp['variancia_k'], 2)}"),
        ("$k_{max} / \\langle k \\rangle$", num(lp["razao_kmax_kmedio"], 2)),
    ]
    salvar("tab_potencia.tex", tabular("lr", ["Medida", "Valor"], linhas))


def tab_robustez():
    r = ler("parte3_robustez.json")
    al = r["aleatoria"]

    def mdp(chave, f=1.0, casas=3):
        xs = [x[chave] * f for x in al]
        return f"{num(statistics.mean(xs), casas)} $\\pm$ {num(statistics.stdev(xs), casas)}"

    o, gr, bc = r["original"], r["grau"], r["betweenness"]
    linhas = [
        ("A: LCC", o["A_lcc"], mdp("A_lcc", casas=0), gr["A_lcc"], bc["A_lcc"]),
        ("A: fração no núcleo", num(o["A_partidas"], 3), mdp("A_partidas"),
         num(gr["A_partidas"], 3), num(bc["A_partidas"], 3)),
        ("B: componentes", o["B_componentes"], mdp("B_componentes", casas=1),
         gr["B_componentes"], bc["B_componentes"]),
        ("C: comp. médio", num(o["C_caminho_medio"], 3), mdp("C_caminho_medio"),
         num(gr["C_caminho_medio"], 3), num(bc["C_caminho_medio"], 3)),
        ("C: eficiência ($\\times 10^{-2}$)", num(o["C_eficiencia"] * 100, 3),
         mdp("C_eficiencia", 100), num(gr["C_eficiencia"] * 100, 3), num(bc["C_eficiencia"] * 100, 3)),
        ("D: isolados (\\%)", num(o["D_isolados"] * 100, 2), mdp("D_isolados", 100, 2),
         num(gr["D_isolados"] * 100, 2), num(bc["D_isolados"] * 100, 2)),
    ]
    salvar("tab_robustez.tex", tabular(
        "lrrrr", ["Métrica", "Original", "Aleatória", "Grau", "Betw."], linhas))


def tab_modelos():
    m = ler("parte3_modelos.json")
    linhas = []
    for nome in ("Observado", "Erdos-Renyi", "Barabasi-Albert", "Watts-Strogatz"):
        x = m[nome]
        rotulo = {"Erdos-Renyi": "Erdős–Rényi", "Barabasi-Albert": "Barabási–Albert"}.get(nome, nome)
        linhas.append((rotulo, num(x["grau_max"], 1 if nome != "Observado" else 0),
                       num(x["var_grau"], 2), num(x["C"], 3), num(x["L"], 2),
                       num(x["componentes"], 1 if nome != "Observado" else 0)))
    salvar("tab_modelos.tex", tabular(
        "lrrrrr", ["Modelo", "$k_{max}$", "Var($k$)", "$C$", "$L$", "Comp."], linhas))


def valores():
    """Macros com valores que variam a cada execução e são citados no texto."""
    t = ler("parte2_tempos.json")
    a = ler("parte2_escala.json")["ajuste"]
    ms = {k: v["media"] * 1000 for k, v in t.items()}
    lineares = ["BFS", "DFS", "Eulerianidade", "Tarjan", "Prim", "Kruskal"]
    macros = {
        "tempoMaxLinear": num(max(ms[k] for k in lineares), 1),
        "tempoBF": num(ms["Bellman-Ford (todas as fontes)"], 1),
        "tempoFW": num(ms["Floyd-Warshall"], 1),
        "razaoBFFW": num(ms["Bellman-Ford (todas as fontes)"] / ms["Floyd-Warshall"], 1),
        "expMinLinear": num(min(a[k]["expoente"] for k in ("BFS", "DFS", "Tarjan")), 2),
        "expMaxLinear": num(max(a[k]["expoente"] for k in ("BFS", "DFS", "Tarjan")), 2),
        "expVDijkstra": num(a["Dijkstra*"]["expoente_v"], 2),
        "expVBF": num(a["Bellman-Ford*"]["expoente_v"], 2),
        "expVFW": num(a["Floyd-Warshall"]["expoente_v"], 2),
        "expFW": num(a["Floyd-Warshall"]["expoente"], 2),
        "expPrim": num(a["Prim"]["expoente"], 2),
        "expKruskal": num(a["Kruskal"]["expoente"], 2),
        "expEuler": num(a["Eulerianidade"]["expoente"], 2),
    }
    sw = ler("parte3_small_world.json")["global"]
    modelos = ler("parte3_modelos.json")
    macros.update({
        "sigmaGlobal": num(sw["er_global"]["sigma"], 0),
        "sigmaBlocos": num(sw["er_blocos"]["sigma"], 2),
        "excessoC": num((sw["er_blocos"]["gamma"] - 1) * 100, 0),
        "fracaoCBlocos": num(sw["er_blocos"]["C_rand"] / sw["C"] * 100, 0),
        "LBlocos": num(sw["er_blocos"]["L_rand"], 2),
        "kmaxBA": num(modelos["Barabasi-Albert"]["grau_max"], 0),
    })
    salvar("valores.tex", "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items()))


def gerar():
    for f in (tab_etapas, tab_metricas, tab_tempos, tab_escala, tab_small_world, tab_potencia,
              tab_robustez, tab_modelos, valores):
        f()


if __name__ == "__main__":
    gerar()
    print(f"Tabelas geradas em {RESULTADOS}")
