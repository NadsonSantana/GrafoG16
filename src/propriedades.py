"""Análise de propriedades da rede (Parte III): small-world, lei de potência,
robustez e comparação com modelos clássicos. Implementações próprias."""

import math
import random
from collections import deque

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.special import zeta

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
    k_medio = 2 * m / n
    resultado = {"C": c_obs, "L": l_obs,
                 "teorico": {"L_rand": math.log(n) / math.log(k_medio), "C_rand": k_medio / n}}
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
    estimador de máxima verossimilhança discreto exato (Clauset et al.) com kmin por KS."""
    dist = distribuicao_graus(g)
    n = g.num_vertices()
    ks = np.array([k for k in dist if k > 0], dtype=float)
    pk = np.array([dist[int(k)] / n for k in ks])
    inclinacao, intercepto = np.polyfit(np.log10(ks), np.log10(pk), 1)
    previsto = inclinacao * np.log10(ks) + intercepto
    ss_res = float(np.sum((np.log10(pk) - previsto) ** 2))
    ss_tot = float(np.sum((np.log10(pk) - np.log10(pk).mean()) ** 2))

    lista = np.array([k for k in graus(g).values() if k > 0], dtype=float)
    mv = melhor_kmin(lista)

    ccdf = {int(k): float((lista >= k).mean()) for k in ks}
    media = float(lista.mean())
    poisson = {int(k): math.exp(-media) * media ** k / math.factorial(int(k)) for k in ks}
    return {
        "ks": ks.tolist(), "pk": pk.tolist(), "ccdf": ccdf, "poisson": poisson,
        "gamma_mq": float(-inclinacao), "constante": float(10 ** intercepto),
        "r2": 1 - ss_res / ss_tot if ss_tot else 0.0, "mv": mv,
        "k_medio": media, "variancia_k": float(lista.var()),
        "razao_kmax_kmedio": float(lista.max() / media),
    }


def _qui_quadrado(observados, esperados, parametros):
    """Teste qui-quadrado de aderência, juntando classes vizinhas até a frequência
    esperada atingir 5. Graus de liberdade: classes - 1 - parâmetros estimados."""
    classes = [[observados[0], esperados[0]]]
    for o, e in zip(observados[1:], esperados[1:]):
        if classes[-1][1] < 5:
            classes[-1][0] += o
            classes[-1][1] += e
        else:
            classes.append([o, e])
    while len(classes) > 1 and classes[-1][1] < 5:
        o, e = classes.pop()
        classes[-1][0] += o
        classes[-1][1] += e
    obs = np.array([c[0] for c in classes], dtype=float)
    esp = np.array([c[1] for c in classes])
    estatistica = float(((obs - esp) ** 2 / esp).sum())
    gl = len(classes) - 1 - parametros
    return {"qui2": estatistica, "gl": gl, "p": float(stats.chi2.sf(estatistica, gl)),
            "classes": len(classes)}


def aderencia_graus(g, blocos):
    """Compara a distribuição de graus com uma Poisson de mesma média e com a mistura de
    binomiais Bin(n_p - 1, p_p) de cada partida (n_p vértices, densidade p_p), que respeita
    o limite k <= n_p - 1. Inclui o teste de dispersão (Var/média = 1 sob Poisson)."""
    ks = np.array([len(g.adj[v]) for v in g.adj])
    n, media, var = len(ks), float(ks.mean()), float(ks.var(ddof=1))
    kmax = int(ks.max()) + 1
    observados = np.bincount(ks, minlength=kmax + 1)[:kmax + 1]
    poisson = stats.poisson.pmf(np.arange(kmax + 1), media)
    poisson[-1] = 1 - poisson[:-1].sum()  # última classe: k >= kmax
    mistura = np.zeros(kmax + 1)
    for bloco in blocos:
        nb = len(bloco)
        if nb < 2:
            mistura[0] += 1
            continue
        pb = 2 * g.subgrafo(bloco).num_arestas() / (nb * (nb - 1))
        pmf = nb * stats.binom.pmf(np.arange(nb), nb - 1, pb)
        mistura[:min(nb, kmax)] += pmf[:kmax]
        mistura[kmax] += pmf[kmax:].sum()  # última classe: k >= kmax
    dispersao = (n - 1) * var / media
    return {
        "indice_dispersao": var / media,
        "dispersao_qui2": dispersao,
        "dispersao_p_sub": float(stats.chi2.cdf(dispersao, n - 1)),  # unilateral: subdispersão
        "poisson": _qui_quadrado(observados, poisson * n, 1),
        "binomial_partidas": _qui_quadrado(observados, mistura, 0),
        "observados": observados.tolist(),
        "esperados_poisson": (poisson * n).tolist(),
        "esperados_binomial": mistura.tolist(),
    }


def log_verossimilhanca(gamma, cauda, kmin):
    """Log-verossimilhança da lei de potência discreta P(k) = k^-gamma / zeta(gamma, kmin)."""
    return -gamma * np.sum(np.log(cauda)) - len(cauda) * math.log(zeta(gamma, kmin))


def ajuste_mv(graus_lista, kmin):
    """MV discreto exato (Clauset, Shalizi e Newman, 2009) na cauda k >= kmin: maximiza
    numericamente a log-verossimilhança com a função zeta de Hurwitz. O erro padrão vem da
    curvatura da log-verossimilhança no ótimo, e a distância de Kolmogorov-Smirnov compara a
    CCDF empírica com P(K >= k) = zeta(gamma, k) / zeta(gamma, kmin). A aproximação contínua
    gamma = 1 + n / sum ln(k / (kmin - 1/2)) é mantida para comparação."""
    cauda = np.sort(graus_lista[graus_lista >= kmin])
    nt = len(cauda)
    otimo = minimize_scalar(lambda g: -log_verossimilhanca(g, cauda, kmin),
                            bounds=(1.0001, 50.0), method="bounded", options={"xatol": 1e-8})
    gamma = float(otimo.x)
    h = 1e-4
    curvatura = (log_verossimilhanca(gamma + h, cauda, kmin) - 2 * log_verossimilhanca(gamma, cauda, kmin)
                 + log_verossimilhanca(gamma - h, cauda, kmin)) / h ** 2
    valores = np.unique(cauda)
    empirica = 1 - np.searchsorted(cauda, valores, side="left") / nt
    modelo = zeta(gamma, valores) / zeta(gamma, kmin)
    return {"kmin": int(kmin), "gamma": gamma, "ks": float(np.max(np.abs(empirica - modelo))),
            "n_cauda": int(nt), "erro_padrao": float(1 / math.sqrt(-curvatura)),
            "gamma_aprox": float(1 + nt / np.sum(np.log(cauda / (kmin - 0.5))))}


def melhor_kmin(graus_lista, cauda_minima=50, valores_minimos=3):
    """Escolhe o kmin de menor distância KS entre os que deixam na cauda ao menos
    cauda_minima vértices e valores_minimos graus distintos (evita caudas degeneradas)."""
    candidatos = [ajuste_mv(graus_lista, k) for k in np.unique(graus_lista)
                  if (graus_lista >= k).sum() >= cauda_minima
                  and len(np.unique(graus_lista[graus_lista >= k])) >= valores_minimos]
    return min(candidatos, key=lambda r: r["ks"])


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


def pagerank(g, amortecimento=0.85, tolerancia=1e-10, max_iter=1000):
    """PageRank por iteração de potência (arestas sem peso). Vértices sem arestas de
    saída distribuem seu valor uniformemente entre todos os vértices."""
    vs = list(g.adj)
    n = len(vs)
    pr = dict.fromkeys(vs, 1 / n)
    for _ in range(max_iter):
        sem_saida = sum(pr[v] for v in vs if not g.adj[v])
        novo = dict.fromkeys(vs, (1 - amortecimento) / n + amortecimento * sem_saida / n)
        for u in vs:
            if g.adj[u]:
                parcela = amortecimento * pr[u] / len(g.adj[u])
                for v in g.adj[u]:
                    novo[v] += parcela
        erro = sum(abs(novo[v] - pr[v]) for v in vs)
        pr = novo
        if erro < tolerancia:
            break
    return pr


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


def robustez(g, fracao=0.05, repeticoes=30, semente=42, dirigido=None,
             fracoes_curva=(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)):
    """Compara remoção aleatória (várias repetições) com ataques direcionados por grau,
    betweenness e PageRank (ranking fixo, calculado no grafo original). O PageRank é
    calculado na versão dirigida, quando fornecida (quem recebe atenção).

    r = ceil(fracao * n), como no Apêndice C. A curva S(f) repete a comparação para
    frações crescentes, com os mesmos rankings fixos."""
    rng = random.Random(semente)
    vs = g.vertices()
    n = g.num_vertices()
    r = math.ceil(fracao * n)
    base = metricas_robustez(g)

    aleatorias = [metricas_robustez(g.remover_vertices(rng.sample(vs, r)), n)
                  for _ in range(repeticoes)]

    gr = graus(g)
    bc = betweenness(g)
    pr = pagerank(dirigido if dirigido is not None else g)
    rankings = {
        "grau": sorted(vs, key=lambda v: (-gr[v], v)),
        "betweenness": sorted(vs, key=lambda v: (-bc[v], v)),
        "pagerank": sorted(vs, key=lambda v: (-pr[v], v)),
    }
    top = {nome: ordem[:r] for nome, ordem in rankings.items()}
    return {
        "r": r,
        "original": base,
        "aleatoria": aleatorias,
        **{nome: metricas_robustez(g.remover_vertices(alvos), n) for nome, alvos in top.items()},
        "partidas_afetadas_pr": len({v.split("_")[0] for v in top["pagerank"]}),
        "sobreposicao_grau_pr": len(set(top["grau"]) & set(top["pagerank"])),
        "partidas_afetadas_grau": len({v.split("_")[0] for v in top["grau"]}),
        "partidas_afetadas_bc": len({v.split("_")[0] for v in top["betweenness"]}),
        "curva": curva_remocao(g, rankings, fracoes_curva, repeticoes, rng),
    }


def curva_remocao(g, rankings, fracoes, repeticoes, rng):
    """S(f): fração no núcleo da partida (A) e de isolados (D) para cada fração removida f,
    na remoção aleatória (média e desvio de várias repetições) e nos ataques de ranking fixo."""
    vs = g.vertices()
    n = g.num_vertices()
    curva = {"fracoes": [0.0] + list(fracoes), "aleatoria": [], "ataques": {k: [] for k in rankings}}
    base = metricas_robustez(g)
    curva["aleatoria"].append({"A_media": base["A_partidas"], "A_dp": 0.0,
                               "D_media": base["D_isolados"], "D_dp": 0.0})
    for nome in rankings:
        curva["ataques"][nome].append({"A": base["A_partidas"], "D": base["D_isolados"]})
    for f in fracoes:
        r = math.ceil(f * n)
        amostras = [metricas_robustez(g.remover_vertices(rng.sample(vs, r)), n)
                    for _ in range(repeticoes)]
        a = [x["A_partidas"] for x in amostras]
        d = [x["D_isolados"] for x in amostras]
        curva["aleatoria"].append({"A_media": float(np.mean(a)), "A_dp": float(np.std(a, ddof=1)),
                                   "D_media": float(np.mean(d)), "D_dp": float(np.std(d, ddof=1))})
        for nome, ordem in rankings.items():
            m = metricas_robustez(g.remover_vertices(ordem[:r]), n)
            curva["ataques"][nome].append({"A": m["A_partidas"], "D": m["D_isolados"]})
    return curva


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

    def ccdf(graus_lista):
        """P(K >= k) para cada grau k presente (graus de todas as amostras juntas)."""
        lista = np.array(graus_lista)
        return {int(k): float((lista >= k).mean()) for k in np.unique(lista) if k > 0}

    saida = {"Observado": resumo(g)}
    ccdfs = {"Observado": ccdf(list(graus(g).values()))}
    for nome, gerar in geradores.items():
        amostras, todos_graus = [], []
        for _ in range(repeticoes):
            h = gerar()
            amostras.append(resumo(h))
            todos_graus.extend(graus(h).values())
        saida[nome] = {k: float(np.mean([a[k] for a in amostras])) for k in amostras[0]}
        ccdfs[nome] = ccdf(todos_graus)
    saida["ccdf_graus"] = ccdfs
    saida["parametros"] = {"m_ba": m_ba, "k_ws": k_ws, "p_ws": 0.1}
    return saida
