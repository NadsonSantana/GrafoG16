"""Implementacoes proprias dos algoritmos da Parte II.

Grafos sao representados por listas de adjacencia (dict):
  nao dirigido: adj[u] = {v: peso, ...} com adj[u][v] == adj[v][u]
  dirigido:     adj[u] = {v: peso, ...} apenas arcos de saida
Nenhuma funcao usa NetworkX; ele e usado apenas para validar os resultados.
"""
import heapq
from collections import deque

import numpy as np

INF = float("inf")


# ---------------------------------------------------------------- buscas
def bfs(adj, origem):
    """Busca em largura a partir de origem. O(V + E).
    Retorna (ordem de visita, distancia em saltos, pai)."""
    dist = {origem: 0}
    pai = {origem: None}
    ordem = []
    fila = deque([origem])
    while fila:
        u = fila.popleft()
        ordem.append(u)
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                pai[v] = u
                fila.append(v)
    return ordem, dist, pai


def bfs_completa(adj):
    """BFS em todas as componentes (floresta de busca). O(V + E)."""
    visitado = set()
    componentes = []
    for s in adj:
        if s not in visitado:
            ordem, _, _ = bfs(adj, s)
            visitado.update(ordem)
            componentes.append(ordem)
    return componentes


def dfs(adj, origem, visitado=None):
    """Busca em profundidade iterativa (pilha explicita). O(V + E).
    Retorna (ordem de descoberta, ordem de termino)."""
    if visitado is None:
        visitado = set()
    descoberta, termino = [], []
    visitado.add(origem)
    descoberta.append(origem)
    pilha = [(origem, iter(adj[origem]))]
    while pilha:
        u, it = pilha[-1]
        avancou = False
        for v in it:
            if v not in visitado:
                visitado.add(v)
                descoberta.append(v)
                pilha.append((v, iter(adj[v])))
                avancou = True
                break
        if not avancou:
            pilha.pop()
            termino.append(u)
    return descoberta, termino


def dfs_completa(adj):
    """DFS em todas as componentes. O(V + E)."""
    visitado = set()
    componentes = []
    for s in adj:
        if s not in visitado:
            desc, _ = dfs(adj, s, visitado)
            componentes.append(desc)
    return componentes


# ---------------------------------------------------------- eulerianidade
def eulerianidade(adj):
    """Classifica um grafo nao dirigido. O(V + E).
    Euleriano: todas as arestas numa unica componente e 0 vertices de grau impar.
    Semi-euleriano: idem, com exatamente 2 vertices de grau impar."""
    impares = [u for u in adj if len(adj[u]) % 2 == 1]
    com_aresta = [u for u in adj if adj[u]]
    if not com_aresta:
        return {"classe": "euleriano (sem arestas)", "impares": 0, "conexo_nas_arestas": True}
    alcancados, _ = dfs(adj, com_aresta[0])
    conexo = len(set(alcancados)) == len(com_aresta)
    if conexo and not impares:
        classe = "euleriano"
    elif conexo and len(impares) == 2:
        classe = "semi-euleriano"
    else:
        classe = "nao euleriano"
    return {"classe": classe, "impares": len(impares), "conexo_nas_arestas": conexo}


# ---------------------------------------------------------- caminhos
def dijkstra(adj, origem):
    """Dijkstra com heap binario (lazy deletion). O((V + E) log V).
    Requer pesos nao negativos."""
    dist = {origem: 0.0}
    pai = {origem: None}
    heap = [(0.0, origem)]
    fechado = set()
    while heap:
        d, u = heapq.heappop(heap)
        if u in fechado:
            continue
        fechado.add(u)
        for v, w in adj[u].items():
            nd = d + w
            if nd < dist.get(v, INF):
                dist[v] = nd
                pai[v] = u
                heapq.heappush(heap, (nd, v))
    return dist, pai


def bellman_ford(vertices, arcos, origem, parada_antecipada=False):
    """Bellman-Ford classico: V-1 rodadas relaxando todos os arcos. O(V * E).
    arcos: lista (u, v, w); grafo nao dirigido entra com os dois sentidos.
    Retorna (dist, pai, tem_ciclo_negativo)."""
    dist = {v: INF for v in vertices}
    pai = {v: None for v in vertices}
    dist[origem] = 0.0
    for _ in range(len(dist) - 1):
        mudou = False
        for u, v, w in arcos:
            du = dist[u]
            if du != INF and du + w < dist[v]:
                dist[v] = du + w
                pai[v] = u
                mudou = True
        if parada_antecipada and not mudou:
            break
    ciclo_negativo = any(dist[u] != INF and dist[u] + w < dist[v] for u, v, w in arcos)
    return dist, pai, ciclo_negativo


def floyd_warshall(vertices, adj):
    """Floyd-Warshall. O(V^3). O laco em k e explicito; a relaxacao de cada
    rodada (V^2 pares) e vetorizada com NumPy. Retorna (indice, matriz D)."""
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    D = np.full((n, n), INF)
    np.fill_diagonal(D, 0.0)
    for u in adj:
        for v, w in adj[u].items():
            D[idx[u], idx[v]] = min(D[idx[u], idx[v]], w)
    for k in range(n):
        np.minimum(D, D[:, k:k + 1] + D[k:k + 1, :], out=D)
    return idx, D


# ---------------------------------------------------------- Tarjan
def tarjan_scc(adj):
    """Componentes fortemente conexas (Tarjan), versao iterativa. O(V + E)."""
    indice, low = {}, {}
    na_pilha, pilha, sccs = set(), [], []
    contador = 0
    for s in adj:
        if s in indice:
            continue
        trabalho = [(s, iter(adj[s]))]
        indice[s] = low[s] = contador
        contador += 1
        pilha.append(s)
        na_pilha.add(s)
        while trabalho:
            u, it = trabalho[-1]
            avancou = False
            for v in it:
                if v not in indice:
                    indice[v] = low[v] = contador
                    contador += 1
                    pilha.append(v)
                    na_pilha.add(v)
                    trabalho.append((v, iter(adj[v])))
                    avancou = True
                    break
                if v in na_pilha:
                    low[u] = min(low[u], indice[v])
            if avancou:
                continue
            trabalho.pop()
            if trabalho:
                p = trabalho[-1][0]
                low[p] = min(low[p], low[u])
            if low[u] == indice[u]:
                comp = []
                while True:
                    w = pilha.pop()
                    na_pilha.discard(w)
                    comp.append(w)
                    if w == u:
                        break
                sccs.append(comp)
    return sccs


def tarjan_articulacoes_pontes(adj):
    """Pontos de articulacao e pontes (Tarjan, lowlink) em grafo nao dirigido. O(V + E)."""
    disc, low = {}, {}
    articulacoes, pontes = set(), []
    tempo = 0
    for s in adj:
        if s in disc:
            continue
        disc[s] = low[s] = tempo
        tempo += 1
        filhos_raiz = 0
        trabalho = [(s, None, iter(adj[s]))]
        while trabalho:
            u, pai_u, it = trabalho[-1]
            avancou = False
            for v in it:
                if v == pai_u:
                    continue
                if v not in disc:
                    disc[v] = low[v] = tempo
                    tempo += 1
                    if u == s:
                        filhos_raiz += 1
                    trabalho.append((v, u, iter(adj[v])))
                    avancou = True
                    break
                low[u] = min(low[u], disc[v])
            if avancou:
                continue
            trabalho.pop()
            if pai_u is not None:
                low[pai_u] = min(low[pai_u], low[u])
                if low[u] > disc[pai_u]:
                    pontes.append((pai_u, u))
                if pai_u != s and low[u] >= disc[pai_u]:
                    articulacoes.add(pai_u)
        if filhos_raiz > 1:
            articulacoes.add(s)
    return articulacoes, pontes


# ---------------------------------------------------------- AGM
class UniaoBusca:
    def __init__(self, elementos):
        self.pai = {e: e for e in elementos}
        self.rank = {e: 0 for e in elementos}

    def busca(self, x):
        while self.pai[x] != x:
            self.pai[x] = self.pai[self.pai[x]]
            x = self.pai[x]
        return x

    def uniao(self, a, b):
        ra, rb = self.busca(a), self.busca(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.pai[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True


def kruskal(vertices, arestas):
    """Floresta geradora minima (Kruskal + uniao-busca). O(E log E).
    arestas: lista (u, v, w) sem duplicatas."""
    ub = UniaoBusca(vertices)
    floresta = []
    for u, v, w in sorted(arestas, key=lambda a: a[2]):
        if ub.uniao(u, v):
            floresta.append((u, v, w))
    return floresta


def prim(adj):
    """Floresta geradora minima (Prim com heap, reiniciado em cada componente).
    O(E log V)."""
    na_arvore = set()
    floresta = []
    for s in adj:
        if s in na_arvore:
            continue
        na_arvore.add(s)
        heap = [(w, s, v) for v, w in adj[s].items()]
        heapq.heapify(heap)
        while heap:
            w, u, v = heapq.heappop(heap)
            if v in na_arvore:
                continue
            na_arvore.add(v)
            floresta.append((u, v, w))
            for x, wx in adj[v].items():
                if x not in na_arvore:
                    heapq.heappush(heap, (wx, v, x))
    return floresta
