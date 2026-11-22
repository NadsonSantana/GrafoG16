"""Tratamento do dataset comm-f2f-Resistance (SNAP).

Pipeline:
  1. Leitura das 62 redes temporais binárias (quem olha para quem a cada 1/3 s).
  2. Grafo temporal -> estático: para cada par ordenado (u, v) da mesma partida,
     calcula a fração do tempo em que u olha para v (atenção f(u, v)).
  3. Remoção do alvo LAPTOP (não é participante) e de auto-laços.
  4. Limiar de atenção: mantém u -> v somente se f(u, v) >= 1/N, onde N é o
     número de alvos possíveis de u (N-1 participantes + laptop), isto é,
     somente se u olha para v mais do que olharia sob atenção uniforme.
  5. Dirigido -> não-direcionado: {u, v} existe se u -> v ou v -> u existir;
     peso w = (f(u, v) + f(v, u)) / 2 e distância d = -ln(w).
  6. União das 62 partidas em um único grafo (vértice = "gXX_pY").
"""

import csv
import json
import math
import os

import numpy as np
import pandas as pd

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORIGEM = os.path.join(RAIZ, "data", "raw", "comm-f2f-Resistance")
DESTINO = os.path.join(RAIZ, "data", "processed")


def nome_vertice(partida, participante):
    return f"g{partida:02d}_p{participante}"


def carregar_partida(partida, n):
    """Retorna array T x n x (n+1); posição 0 do último eixo é o laptop."""
    df = pd.read_csv(os.path.join(ORIGEM, "network", f"network{partida}.csv"), index_col=0)
    return df.values.reshape(len(df), n, n + 1)


def matriz_atencao(dados):
    """f[i, j]: fração dos instantes em que i olha para j (laptop e auto-laços excluídos)."""
    atencao = dados[:, :, 1:].mean(axis=0)
    np.fill_diagonal(atencao, 0.0)
    return atencao


def foco_por_instante(dados):
    """Em cada instante, fração dos olhares a participantes recebida pelo mais olhado."""
    pessoas = dados[:, :, 1:].astype(float)
    recebidos = pessoas.sum(axis=1)  # T x n
    total = recebidos.sum(axis=1)
    return np.where(total > 0, recebidos.max(axis=1) / np.maximum(total, 1), np.nan)


def foco_uniforme(n, rng, instantes=4000):
    """Valor esperado do foco se cada participante escolhesse um alvo ao acaso entre os
    n alvos possíveis (n-1 colegas + laptop), estimado por simulação."""
    dados = np.zeros((instantes, n, n + 1))
    for i in range(n):
        alvos = rng.integers(0, n, size=instantes)  # 0 = laptop; j >= 1 = colega
        alvos = np.where((alvos >= 1) & (alvos >= i + 1), alvos + 1, alvos)  # pula o próprio i
        dados[np.arange(instantes), i, np.minimum(alvos, n)] = 1
    return float(np.nanmean(foco_por_instante(dados)))


def serie_temporal(partida, dados, foco_ao_acaso, faixas=20):
    """Descrição do grafo original ao longo da partida (tempo normalizado em faixas):
    fração dos olhares dirigidos ao laptop e foco (fração recebida pelo mais olhado)."""
    foco = foco_por_instante(dados)
    linhas = []
    for f, idx in enumerate(np.array_split(np.arange(dados.shape[0]), faixas)):
        linhas.append({"partida": partida, "faixa": f,
                       "laptop": round(float(dados[idx, :, 0].mean()), 6),
                       "foco": round(float(np.nanmean(foco[idx])), 6),
                       "foco_uniforme": round(foco_ao_acaso, 6)})
    return linhas


def arestas_partida(partida, atencao, limiar=None):
    """Aplica o limiar (padrão 1/N, N = n alvos: n-1 colegas + laptop) e a simetrização.

    Retorna (arcos dirigidos, arestas não-direcionadas), cada um como lista de dicionários
    com atenção e distância d = -ln(atenção).
    """
    n = atencao.shape[0]
    if limiar is None:
        limiar = 1.0 / n
    arcos = []
    for i in range(n):
        for j in range(n):
            if i != j and atencao[i, j] >= limiar:
                arcos.append({
                    "origem": nome_vertice(partida, i + 1), "destino": nome_vertice(partida, j + 1),
                    "atencao": round(float(atencao[i, j]), 6),
                    "distancia": round(-math.log(atencao[i, j]), 6),
                })
    arestas = []
    for i in range(n):
        for j in range(i + 1, n):
            if atencao[i, j] >= limiar or atencao[j, i] >= limiar:
                w = (atencao[i, j] + atencao[j, i]) / 2
                arestas.append({
                    "origem": nome_vertice(partida, i + 1), "destino": nome_vertice(partida, j + 1),
                    "atencao": round(float(w), 6), "distancia": round(-math.log(w), 6),
                })
    return arcos, arestas


def tratar():
    os.makedirs(DESTINO, exist_ok=True)
    meta = pd.read_csv(os.path.join(ORIGEM, "network_list.csv"))

    vertices = []
    pares = []
    serie = []
    rng = np.random.default_rng(42)
    dirigidas = []
    nao_dirigidas = []
    partidas = []
    etapas = {
        "snapshots": 0,
        "arestas_temporais": 0,
        "arestas_temporais_laptop": 0,
        "pares_estaticos": 0,
        "pares_apos_limiar": 0,
    }

    for _, linha in meta.iterrows():
        g, n = int(linha.NETWORK), int(linha.NUMBER_OF_PARTICIPANTS)
        dados = carregar_partida(g, n)
        t = dados.shape[0]
        etapas["snapshots"] += t
        etapas["arestas_temporais"] += int(dados[:, :, 1:].sum())
        etapas["arestas_temporais_laptop"] += int(dados[:, :, 0].sum())

        atencao = matriz_atencao(dados)
        etapas["pares_estaticos"] += int((atencao > 0).sum())

        partidas.append({"partida": g, "participantes": n, "duracao_s": round(t / 3, 1),
                         "atencao_laptop": round(float(dados[:, :, 0].mean()), 4)})
        for i in range(n):
            vertices.append({"id": nome_vertice(g, i + 1), "partida": g, "participante": i + 1})
        for i in range(n):
            for j in range(n):
                if i != j:
                    pares.append({"partida": g, "participantes": n, "origem": i + 1, "destino": j + 1,
                                  "atencao": round(float(atencao[i, j]), 6)})
        serie.extend(serie_temporal(g, dados, foco_uniforme(n, rng)))
        arcos, arestas = arestas_partida(g, atencao)
        dirigidas.extend(arcos)
        nao_dirigidas.extend(arestas)

    etapas["pares_apos_limiar"] = len(dirigidas)
    etapas["arestas_nao_dirigidas"] = len(nao_dirigidas)
    etapas["vertices"] = len(vertices)
    etapas["partidas"] = len(partidas)

    _salvar_csv("vertices.csv", vertices)
    _salvar_csv("arestas_dirigidas.csv", dirigidas)
    _salvar_csv("arestas.csv", nao_dirigidas)
    _salvar_csv("partidas.csv", partidas)
    _salvar_csv("atencao_pares.csv", pares)  # todos os pares ordenados, antes do limiar
    _salvar_csv("serie_temporal.csv", serie)  # laptop e foco ao longo de cada partida
    with open(os.path.join(DESTINO, "etapas_tratamento.json"), "w") as f:
        json.dump(etapas, f, indent=2)
    return etapas


def _salvar_csv(nome, linhas):
    with open(os.path.join(DESTINO, nome), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)


if __name__ == "__main__":
    print(json.dumps(tratar(), indent=2))
