"""Etapa 1 - Inspecao dos dados brutos do G16 (comm-f2f-Resistance).

Le os 62 arquivos de cada versao (binaria e ponderada) e gera estatisticas
brutas por partida em data/processed/inspecao_bruta.csv.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
DIR_RAW = RAIZ / "data" / "raw" / "comm-f2f-Resistance"
DIR_NET = DIR_RAW / "network"
SAIDA = RAIZ / "data" / "processed" / "inspecao_bruta.csv"

COL = re.compile(r"^P(\d+)_TO_(P(\d+)|LAPTOP)$")


def colunas_validas(df, n_part):
    """Colunas P_i -> P_j entre participantes reais (i, j <= n_part, i != j)."""
    pares, laptop, auto = [], [], []
    for c in df.columns:
        m = COL.match(c)
        if not m:
            continue
        i = int(m.group(1))
        if i > n_part:
            continue
        if m.group(2) == "LAPTOP":
            laptop.append(c)
        else:
            j = int(m.group(3))
            if j > n_part:
                continue
            (auto if i == j else pares).append(c)
    return pares, laptop, auto


def main():
    lista = pd.read_csv(DIR_RAW / "network_list.csv")
    linhas = []
    for net, n_part in lista.itertuples(index=False):
        b = pd.read_csv(DIR_NET / f"network{net}.csv")
        w = pd.read_csv(DIR_NET / f"network{net}_weighted.csv")
        pares, laptop, auto = colunas_validas(b, n_part)
        tempo = b["TIME"].to_numpy()
        passos = np.diff(tempo)

        bp = b[pares].to_numpy()
        wp = w[pares].to_numpy()
        # Colunas de participantes inexistentes (P > n_part) devem ser vazias
        extras = [c for c in b.columns if c != "TIME" and c not in pares + laptop + auto]
        linhas.append({
            "rede": net,
            "participantes": n_part,
            "instantes": len(b),
            "tempo_min": tempo.min(),
            "tempo_max": tempo.max(),
            "passo_mediano": float(np.median(passos)) if len(passos) else np.nan,
            "linhas_iguais_bin_w": len(b) == len(w),
            "arestas_temporais_bin": int(bp.sum()),
            "pares_dirigidos_possiveis": n_part * (n_part - 1),
            "pares_dirigidos_ativos_bin": int((bp.sum(axis=0) > 0).sum()),
            "olhar_laptop_bin": int(b[laptop].to_numpy().sum()),
            "autoloops_bin": int(b[auto].to_numpy().sum()),
            "soma_extras_bin": float(b[extras].to_numpy().sum()) if extras else 0.0,
            "nan_bin": int(b.isna().sum().sum()),
            "nan_w": int(w.isna().sum().sum()),
            "w_min": float(np.nanmin(wp)),
            "w_max": float(np.nanmax(wp)),
            "w_medio": float(np.nanmean(wp)),
        })
    df = pd.DataFrame(linhas)
    SAIDA.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(SAIDA, index=False)

    pd.set_option("display.width", 200)
    print(df.describe().T[["mean", "min", "max"]].round(3))
    print("\nTotais:")
    print(f"  redes (partidas)              : {len(df)}")
    print(f"  vertices (soma participantes) : {df.participantes.sum()}")
    print(f"  instantes (soma)              : {df.instantes.sum()}")
    print(f"  arestas temporais binarias    : {df.arestas_temporais_bin.sum()} (so P->P)")
    print(f"  olhares para o laptop         : {df.olhar_laptop_bin.sum()}")
    print(f"  P->P + laptop                 : {df.arestas_temporais_bin.sum() + df.olhar_laptop_bin.sum()}")
    print(f"  duracao total (s, max-min+1)  : {(df.tempo_max - df.tempo_min + 1).sum():.0f}")
    print(f"  pares dirigidos ativos/possiveis: {df.pares_dirigidos_ativos_bin.sum()}/{df.pares_dirigidos_possiveis.sum()}")
    print(f"  autoloops={df.autoloops_bin.sum()}  NaN bin={df.nan_bin.sum()}  NaN w={df.nan_w.sum()}  extras={df.soma_extras_bin.sum()}")
    print(f"  linhas bin==w em todas: {df.linhas_iguais_bin_w.all()}")
    print(f"\nSalvo em {SAIDA.relative_to(RAIZ)}")


if __name__ == "__main__":
    main()
