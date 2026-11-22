"""Parte III: small-world, lei de potência, robustez e comparação com modelos clássicos.

Uso: python -m src.parte3 [--rapido | --so-figuras]"""


import numpy as np

from . import algoritmos as alg
from . import propriedades as prop
from .sensibilidade import sensibilidade_limiar
from .comum import (AQUA, AZUL, LARANJA, TINTA, TINTA_2, VERMELHO, carregar, eixo_graus_log,
                    ler_json, plt, preparar, salvar_figura, salvar_json, so_figuras)


def parte3(g, d, rapido=False):
    blocos = alg.componentes_conexas(g)
    reps = 5 if rapido else 30
    sw = prop.small_world(g, blocos, repeticoes=reps)
    lp = prop.lei_de_potencia(g)
    lp["aderencia"] = prop.aderencia_graus(g, blocos)
    rob = prop.robustez(g, repeticoes=reps, dirigido=d)
    modelos = prop.comparar_modelos(g, repeticoes=3 if rapido else 10)
    por_partida = small_world_por_partida(g, blocos)
    sens = sensibilidade_limiar(repeticoes=3 if rapido else 10)
    salvar_json("parte3_small_world.json", {"global": sw, "por_partida": por_partida})
    salvar_json("parte3_lei_potencia.json", lp)
    salvar_json("parte3_robustez.json", rob)
    salvar_json("parte3_modelos.json", modelos)
    salvar_json("parte3_sensibilidade.json", sens)
    figura_lei_potencia(lp)
    figura_robustez(rob)
    figura_curva_remocao(rob["curva"])
    figura_modelos(modelos["ccdf_graus"])
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
    salvar_figura(fig, "lei_potencia")


def figura_modelos(ccdfs):
    """CCDF dos graus do grafo observado contra os modelos clássicos (amostras agregadas)."""
    estilos = {"Observado": (AZUL, "o", "observado (G16)"),
               "Erdos-Renyi": (LARANJA, "s", "Erdős–Rényi $G(n,m)$"),
               "Barabasi-Albert": (AQUA, "^", "Barabási–Albert ($m_0 = 2$)"),
               "Watts-Strogatz": (VERMELHO, "D", "Watts–Strogatz ($p = 0{,}1$)")}
    fig, ax = plt.subplots(figsize=(3.3, 2.6))
    for nome, (cor, marcador, rotulo) in estilos.items():
        pontos = sorted((int(k), v) for k, v in ccdfs[nome].items())
        ax.loglog([k for k, _ in pontos], [v for _, v in pontos], marker=marcador, color=cor,
                  markersize=3.5, linewidth=1.5 if nome == "Observado" else 1, label=rotulo)
    ax.set_xlabel("grau $k$ (log)")
    ax.set_ylabel("$P(K \\geq k)$ (log)")
    ax.legend(fontsize=6, loc="lower left")
    salvar_figura(fig, "modelos_ccdf")


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
        ax.plot(1.7, rob["grau"][chave], "D", color=VERMELHO, markersize=6)
        ax.plot(1.95, rob["betweenness"][chave], "^", color=LARANJA, markersize=7)
        ax.plot(2.2, rob["pagerank"][chave], "s", color=AQUA, markersize=6)
        ax.set_xlim(0.5, 2.5)
        ax.set_xticks([1, 1.95], ["aleat.", "direc."])
        ax.set_title(titulo, fontsize=8)
        ax.grid(axis="x", visible=False)
    handles = [
        plt.Line2D([], [], color=TINTA_2, marker="_", linestyle="", markersize=12, label="original"),
        plt.Line2D([], [], color=TINTA, marker="o", linestyle="", label="média ± dp (aleatória)"),
        plt.Line2D([], [], color=VERMELHO, marker="D", linestyle="", label="ataque por grau"),
        plt.Line2D([], [], color=LARANJA, marker="^", linestyle="", label="ataque por betweenness"),
        plt.Line2D([], [], color=AQUA, marker="s", linestyle="", label="ataque por PageRank"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7, bbox_to_anchor=(0.5, -0.14))
    fig.tight_layout()
    salvar_figura(fig, "robustez")


def figura_curva_remocao(curva):
    """S(f): fração no núcleo da partida (A) e fração de isolados (D) em função da fração
    removida, para a remoção aleatória (média ± desvio) e os ataques de ranking fixo."""
    fs = np.array(curva["fracoes"]) * 100
    estilos = {"grau": (VERMELHO, "D", "grau"), "betweenness": (LARANJA, "^", "betweenness"),
               "pagerank": (AQUA, "s", "PageRank")}
    fig, eixos = plt.subplots(1, 2, figsize=(6.4, 2.5))
    for ax, (chave_al, chave_at, rotulo) in zip(
            eixos, [("A", "A", "A: fração no núcleo da partida"), ("D", "D", "D: fração de isolados")]):
        media = np.array([x[f"{chave_al}_media"] for x in curva["aleatoria"]])
        dp = np.array([x[f"{chave_al}_dp"] for x in curva["aleatoria"]])
        ax.fill_between(fs, media - dp, media + dp, color=AZUL, alpha=0.2, linewidth=0)
        ax.plot(fs, media, "o-", color=AZUL, markersize=3.5, linewidth=1.8, label="aleatória (média ± dp)")
        for nome, (cor, mk, rot) in estilos.items():
            ax.plot(fs, [x[chave_at] for x in curva["ataques"][nome]], marker=mk, color=cor,
                    markersize=3.5, linewidth=1.2, label=f"ataque por {rot}")
        if chave_al == "A":
            ax.plot(fs, 1 - fs / 100, ":", color=TINTA_2, linewidth=1, label="limite $1 - f$")
        ax.set_xlabel("vértices removidos $f$ (%)")
        ax.set_title(rotulo, fontsize=8)
        ax.set_xticks(fs[::2])
    eixos[0].legend(fontsize=6, loc="lower left")
    fig.tight_layout()
    salvar_figura(fig, "robustez_curva")


def figuras():
    """Redesenha as figuras da Parte III a partir dos JSONs de results/."""
    figura_lei_potencia(ler_json("parte3_lei_potencia.json"))
    rob = ler_json("parte3_robustez.json")
    figura_robustez(rob)
    figura_curva_remocao(rob["curva"])
    figura_modelos(ler_json("parte3_modelos.json")["ccdf_graus"])


def main():
    rapido = preparar()
    if so_figuras():
        print("Parte III - redesenhando figuras")
        figuras()
        return
    g, d = carregar()
    print("Parte III - propriedades da rede")
    parte3(g, d, rapido)


if __name__ == "__main__":
    main()
