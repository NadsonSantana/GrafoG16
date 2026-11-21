"""Análise de propriedades da rede (Parte III): small-world, lei de potência,
robustez e comparação com modelos clássicos. Implementações próprias."""

import math
import random
from collections import deque

import numpy as np

from .algoritmos import bfs, componentes_conexas
from .grafo import Grafo
from .metricas import caminhos_medios, clusterizacao_media, distribuicao_graus, graus


# ---------------------------------------------------------------- geradores

def gerar_erdos_renyi(n, m, rng):
    """G(n, m): m arestas distintas escolhidas uniformemente entre os n(n-1)/2 pares."""
    g = Grafo()
    for v in range(n):
        g.adicionar_vertice(v)
    while g.num_arestas() < m:
        u, v = rng.randrange(n), rng.randrange(n)
        if u != v and v not in g.adj[u]:
            g.adicionar_aresta(u, v)
    return g


def gerar_er_por_blocos(g, blocos, rng):
    """Aleatório equivalente que preserva as partidas: dentro de cada bloco, sorteia
    o mesmo número de arestas que o grafo original possui naquele bloco."""
    h = Grafo()
    for bloco in blocos:
        vs = list(bloco)
        m_bloco = g.subgrafo(vs).num_arestas()
        for v in vs:
            h.adicionar_vertice(v)
        pares = [(vs[i], vs[j]) for i in range(len(vs)) for j in range(i + 1, len(vs))]
        for u, v in rng.sample(pares, m_bloco):
            h.adicionar_aresta(u, v)
    return h


def gerar_barabasi_albert(n, m_novo, rng):
    """Anexação preferencial: cada novo vértice liga-se a m_novo vértices existentes
    com probabilidade proporcional ao grau."""
    g = Grafo()
    for v in range(m_novo + 1):
        for u in range(v):
            g.adicionar_aresta(u, v)
    repeticoes = [v for v in g.adj for _ in range(len(g.adj[v]))]
    for v in range(m_novo + 1, n):
        alvos = set()
        while len(alvos) < m_novo:
            alvos.add(rng.choice(repeticoes))
        for u in alvos:
            g.adicionar_aresta(u, v)
            repeticoes.extend([u, v])
    return g


def gerar_watts_strogatz(n, k, p, rng):
    """Anel regular com k vizinhos por vértice e religação com probabilidade p."""
    g = Grafo()
    for v in range(n):
        g.adicionar_vertice(v)
    for v in range(n):
        for j in range(1, k // 2 + 1):
            g.adicionar_aresta(v, (v + j) % n)
    for v in range(n):
        for j in range(1, k // 2 + 1):
            u = (v + j) % n
            if rng.random() < p and u in g.adj[v]:
                candidatos = [w for w in range(n) if w != v and w not in g.adj[v]]
                if candidatos:
                    w = rng.choice(candidatos)
                    del g.adj[v][u]
                    del g.adj[u][v]
                    g.adicionar_aresta(v, w)
    return g


# ---------------------------------------------------------------- small-world

def small_world(g, blocos, repeticoes=30, semente=42):
    """sigma = (C / C_rand) / (L / L_rand), contra dois aleatórios equivalentes:
    G(n, m) global e G(n, m) preservando as partidas (blocos)."""
    rng = random.Random(semente)
    n, m = g.num_vertices(), g.num_arestas()
    c_obs, l_obs = clusterizacao_media(g), caminhos_medios(g)
    resultado = {"C": c_obs, "L": l_obs}
    for nome, gerador in (("er_global", lambda: gerar_erdos_renyi(n, m, rng)),
                          ("er_blocos", lambda: gerar_er_por_blocos(g, blocos, rng))):
        cs, ls, comps = [], [], []
        for _ in range(repeticoes):
            r = gerador()
            cs.append(clusterizacao_media(r))
            ls.append(caminhos_medios(r))
            comps.append(len(componentes_conexas(r)))
        c_r, l_r = float(np.mean(cs)), float(np.mean(ls))
        resultado[nome] = {
            "C_rand": c_r, "C_rand_dp": float(np.std(cs, ddof=1)),
            "L_rand": l_r, "L_rand_dp": float(np.std(ls, ddof=1)),
            "componentes_rand": float(np.mean(comps)),
            "gamma": c_obs / c_r, "lambda": l_obs / l_r,
            "sigma": (c_obs / c_r) / (l_obs / l_r),
        }
    return resultado


# ---------------------------------------------------------------- lei de potência

def lei_de_potencia(g):
    """Ajuste P(k) ~ C k^-gamma por mínimos quadrados em log-log, CCDF e
    estimador de máxima verossimilhança discreto (aproximação de Clauset et al.)."""
    dist = distribuicao_graus(g)
    n = g.num_vertices()
    ks = np.array([k for k in dist if k > 0], dtype=float)
    pk = np.array([dist[int(k)] / n for k in ks])
    inclinacao, intercepto = np.polyfit(np.log10(ks), np.log10(pk), 1)
    previsto = inclinacao * np.log10(ks) + intercepto
    ss_res = float(np.sum((np.log10(pk) - previsto) ** 2))
    ss_tot = float(np.sum((np.log10(pk) - np.log10(pk).mean()) ** 2))

    lista = np.array([k for k in graus(g).values() if k > 0], dtype=float)
    kmin = lista.min()
    gamma_mv = 1 + len(lista) / np.sum(np.log(lista / (kmin - 0.5)))

    ccdf = {int(k): float((lista >= k).mean()) for k in ks}
    media = float(lista.mean())
    poisson = {int(k): math.exp(-media) * media ** k / math.factorial(int(k)) for k in ks}
    return {
        "ks": ks.tolist(), "pk": pk.tolist(), "ccdf": ccdf, "poisson": poisson,
        "gamma_mq": float(-inclinacao), "constante": float(10 ** intercepto),
        "r2": 1 - ss_res / ss_tot if ss_tot else 0.0, "gamma_mv": float(gamma_mv),
        "k_medio": media, "variancia_k": float(lista.var()),
        "razao_kmax_kmedio": float(lista.max() / media),
    }


# ---------------------------------------------------------------- centralidade

def betweenness(g):
    """Centralidade de intermediação (Brandes, não ponderada). O(V E)."""
    cb = {v: 0.0 for v in g.adj}
    for s in g.adj:
        pilha = []
        pred = {v: [] for v in g.adj}
        sigma = dict.fromkeys(g.adj, 0)
        sigma[s] = 1
        dist = {s: 0}
        fila = deque([s])
        while fila:
            v = fila.popleft()
            pilha.append(v)
            for w in g.adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1
                    fila.append(w)
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = dict.fromkeys(g.adj, 0.0)
        while pilha:
            w = pilha.pop()
            for v in pred[w]:
                delta[v] += sigma[v] / sigma[w] * (1 + delta[w])
            if w != s:
                cb[w] += delta[w]
    return {v: c / 2 for v, c in cb.items()}  # cada par contado duas vezes


# ---------------------------------------------------------------- robustez

def eficiencia_global(g):
    """Média de 1/d(u, v) sobre todos os pares (0 para pares desconectados)."""
    n = g.num_vertices()
    if n < 2:
        return 0.0
    soma = 0.0
    for v in g.adj:
        for u, d in bfs(g, v).items():
            if u != v:
                soma += 1 / d
    return soma / (n * (n - 1))


def metricas_robustez(g, n_original=None):
    """A: maior componente; A_partidas: fração dos vértices originais que permanecem
    no maior pedaço conexo da sua partida; B: componentes; C: distâncias; D: isolados."""
    comps = componentes_conexas(g)
    n = g.num_vertices()
    n_original = n_original or n
    maior_por_partida = {}
    for c in comps:
        partida = str(c[0]).split("_")[0]
        maior_por_partida[partida] = max(maior_por_partida.get(partida, 0), len(c))
    return {
        "A_lcc": len(comps[0]) if comps else 0,
        "A_partidas": sum(maior_por_partida.values()) / n_original,
        "B_componentes": len(comps),
        "C_caminho_medio": caminhos_medios(g),
        "C_eficiencia": eficiencia_global(g),
        "D_isolados": sum(1 for v in g.adj if not g.adj[v]) / n if n else 0.0,
    }


def robustez(g, fracao=0.05, repeticoes=30, semente=42):
    """Compara remoção aleatória (várias repetições) com ataques direcionados
    por grau e por betweenness (ranking fixo, calculado no grafo original)."""
    rng = random.Random(semente)
    r = int(math.floor(fracao * g.num_vertices()))
    vs = g.vertices()
    n = g.num_vertices()
    base = metricas_robustez(g)

    aleatorias = [metricas_robustez(g.remover_vertices(rng.sample(vs, r)), n)
                  for _ in range(repeticoes)]

    gr = graus(g)
    top_grau = sorted(vs, key=lambda v: (-gr[v], v))[:r]
    bc = betweenness(g)
    top_bc = sorted(vs, key=lambda v: (-bc[v], v))[:r]
    return {
        "r": r,
        "original": base,
        "aleatoria": aleatorias,
        "grau": metricas_robustez(g.remover_vertices(top_grau), n),
        "betweenness": metricas_robustez(g.remover_vertices(top_bc), n),
        "partidas_afetadas_grau": len({v.split("_")[0] for v in top_grau}),
        "partidas_afetadas_bc": len({v.split("_")[0] for v in top_bc}),
    }


# ---------------------------------------------------------------- modelos clássicos

def comparar_modelos(g, repeticoes=10, semente=42):
    """Gera ER, BA e WS com mesmo n e (aproximadamente) mesmo m e compara métricas."""
    rng = random.Random(semente)
    n, m = g.num_vertices(), g.num_arestas()
    k_medio = 2 * m / n
    m_ba = max(1, round(m / n))
    k_ws = max(2, 2 * round(k_medio / 2))
    geradores = {
        "Erdos-Renyi": lambda: gerar_erdos_renyi(n, m, rng),
        "Barabasi-Albert": lambda: gerar_barabasi_albert(n, m_ba, rng),
        "Watts-Strogatz": lambda: gerar_watts_strogatz(n, k_ws, 0.1, rng),
    }

    def resumo(h):
        gr = list(graus(h).values())
        return {"arestas": h.num_arestas(), "grau_max": max(gr), "var_grau": float(np.var(gr)),
                "C": clusterizacao_media(h), "L": caminhos_medios(h),
                "componentes": len(componentes_conexas(h))}

    saida = {"Observado": resumo(g)}
    for nome, gerar in geradores.items():
        amostras = [resumo(gerar()) for _ in range(repeticoes)]
        saida[nome] = {k: float(np.mean([a[k] for a in amostras])) for k in amostras[0]}
    saida["parametros"] = {"m_ba": m_ba, "k_ws": k_ws, "p_ws": 0.1}
    return saida
