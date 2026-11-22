"""Configuração compartilhada pelas Partes I, II e III: caminhos, estilo dos gráficos,
carga do grafo tratado e gravação de resultados."""

import json
import os
import sys

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .grafo import Grafo  # noqa: E402

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


def ler_json(nome):
    with open(os.path.join(RESULTADOS, nome)) as f:
        return json.load(f)


def so_figuras():
    """--so-figuras: redesenha as figuras a partir dos JSONs de results/, sem recalcular."""
    return "--so-figuras" in sys.argv


def salvar_json(nome, dados):
    with open(os.path.join(RESULTADOS, nome), "w") as f:
        json.dump(dados, f, indent=2, ensure_ascii=False, default=str)


def salvar_figura(fig, nome):
    """Salva a figura em PDF vetorial (relatório) e em PNG (README e notebook)."""
    # sem data de criação/software nos metadados: redesenhar não altera os arquivos
    fig.savefig(os.path.join(FIGURAS, f"{nome}.pdf"), metadata={"CreationDate": None})
    fig.savefig(os.path.join(FIGURAS, f"{nome}.png"), metadata={"Software": None})
    plt.close(fig)


def preparar():
    """Cria as pastas de saída e indica se o modo rápido (--rapido) foi pedido."""
    os.makedirs(RESULTADOS, exist_ok=True)
    os.makedirs(FIGURAS, exist_ok=True)
    return "--rapido" in sys.argv
