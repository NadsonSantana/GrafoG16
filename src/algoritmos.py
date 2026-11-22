"""Algoritmos clássicos implementados sobre a classe Grafo (sem bibliotecas prontas)."""

import heapq
from collections import deque

INF = float("inf")


def bfs(g, origem):
    """Busca em largura. Retorna {vértice: distância em saltos} dos alcançáveis. O(V + E)."""
    dist = {origem: 0}
    fila = deque([origem])
    while fila:
        u = fila.popleft()
        for v in g.adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                fila.append(v)
    return dist


def dfs(g, origem):
    """Busca em profundidade iterativa (pilha explícita). Retorna a ordem de visita. O(V + E)."""
    visitados = set()
    ordem = []
    pilha = [origem]
    while pilha:
        u = pilha.pop()
        if u in visitados:
            continue
        visitados.add(u)
        ordem.append(u)
        for v in reversed(list(g.adj[u])):
            if v not in visitados:
                pilha.append(v)
    return ordem


def componentes_conexas(g):
    """Componentes de um grafo não-direcionado via BFS, da maior para a menor. O(V + E)."""
    vistos = set()
    comps = []
    for s in g.adj:
        if s not in vistos:
            comp = list(bfs(g, s))
            vistos.update(comp)
            comps.append(comp)
    comps.sort(key=len, reverse=True)
    return comps


def verificar_euleriano(g):
    """Classifica o grafo não-direcionado como 'circuito', 'caminho' ou 'nao_euleriano'. O(V + E).

    Condições: as arestas devem estar em uma única componente e o número de
    vértices de grau ímpar deve ser 0 (circuito) ou 2 (caminho).
    """
    impares = [v for v in g.adj if len(g.adj[v]) % 2 == 1]
    com_arestas = [v for v in g.adj if g.adj[v]]
    conexo = True
    if com_arestas:
        alcancados = bfs(g, com_arestas[0])
        conexo = all(v in alcancados for v in com_arestas)
    if not conexo:
        tipo = "nao_euleriano"
    elif len(impares) == 0:
        tipo = "circuito"
    elif len(impares) == 2:
        tipo = "caminho"
    else:
        tipo = "nao_euleriano"
    return {"tipo": tipo, "conexo": conexo, "vertices_grau_impar": len(impares)}


def dijkstra(g, origem):
    """Caminhos mínimos de fonte única com heap binário. O((V + E) log V)."""
    dist = {origem: 0.0}
    heap = [(0.0, origem)]
    finalizados = set()
    while heap:
        d, u = heapq.heappop(heap)
        if u in finalizados:
            continue
        finalizados.add(u)
        for v, w in g.adj[u].items():
            nd = d + w
            if nd < dist.get(v, INF):
                dist[v] = nd
                heapq.heappush(heap, (nd, v))
    return dist


def bellman_ford(g, origem):
    """Caminhos mínimos com suporte a pesos negativos. O(V * E).

    Retorna (dist, tem_ciclo_negativo). Encerra cedo se uma rodada não relaxar nada.
    """
    arestas = [(u, v, w) for u in g.adj for v, w in g.adj[u].items()]
    dist = {v: INF for v in g.adj}
    dist[origem] = 0.0
    for _ in range(len(g.adj) - 1):
        mudou = False
        for u, v, w in arestas:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                mudou = True
        if not mudou:
            break
    ciclo_negativo = any(dist[u] + w < dist[v] for u, v, w in arestas)
    return {v: d for v, d in dist.items() if d < INF}, ciclo_negativo


def floyd_warshall(g):
    """Caminhos mínimos entre todos os pares. O(V^3). Retorna (vértices, matriz)."""
    vs = list(g.adj)
    idx = {v: i for i, v in enumerate(vs)}
    n = len(vs)
    d = [[INF] * n for _ in range(n)]
    for i in range(n):
        d[i][i] = 0.0
    for u in vs:
        for v, w in g.adj[u].items():
            d[idx[u]][idx[v]] = min(d[idx[u]][idx[v]], w)
    for k in range(n):
        dk = d[k]
        for i in range(n):
            dik = d[i][k]
            if dik == INF:
                continue
            di = d[i]
            for j in range(n):
                nd = dik + dk[j]
                if nd < di[j]:
                    di[j] = nd
    return vs, d


def floyd_warshall_numpy(g):
    """Floyd-Warshall vetorizado com numpy: cada iteração de k atualiza a matriz inteira
    de uma vez. Mesma complexidade O(V^3), sem poda. Retorna (vértices, matriz numpy)."""
    import numpy as np

    vs = list(g.adj)
    idx = {v: i for i, v in enumerate(vs)}
    n = len(vs)
    d = np.full((n, n), np.inf)
    np.fill_diagonal(d, 0.0)
    for u in vs:
        for v, w in g.adj[u].items():
            d[idx[u], idx[v]] = min(d[idx[u], idx[v]], w)
    for k in range(n):
        np.minimum(d, d[:, k, None] + d[None, k, :], out=d)
    return vs, d


def tarjan(g):
    """Componentes fortemente conexas (Tarjan), versão iterativa. O(V + E)."""
    indice = {}
    baixo = {}
    na_pilha = set()
    pilha = []
    sccs = []
    contador = 0
    for raiz in g.adj:
        if raiz in indice:
            continue
        trabalho = [(raiz, iter(g.adj[raiz]))]
        indice[raiz] = baixo[raiz] = contador
        contador += 1
        pilha.append(raiz)
        na_pilha.add(raiz)
        while trabalho:
            u, it = trabalho[-1]
            avancou = False
            for v in it:
                if v not in indice:
                    indice[v] = baixo[v] = contador
                    contador += 1
                    pilha.append(v)
                    na_pilha.add(v)
                    trabalho.append((v, iter(g.adj[v])))
                    avancou = True
                    break
                if v in na_pilha:
                    baixo[u] = min(baixo[u], indice[v])
            if avancou:
                continue
            trabalho.pop()
            if trabalho:
                pai = trabalho[-1][0]
                baixo[pai] = min(baixo[pai], baixo[u])
            if baixo[u] == indice[u]:
                comp = []
                while True:
                    w = pilha.pop()
                    na_pilha.discard(w)
                    comp.append(w)
                    if w == u:
                        break
                sccs.append(comp)
    return sccs


def prim(g):
    """Floresta geradora mínima (Prim com heap, reiniciado em cada componente). O(E log V)."""
    visitados = set()
    arvore = []
    for raiz in g.adj:
        if raiz in visitados:
            continue
        visitados.add(raiz)
        heap = [(w, raiz, v) for v, w in g.adj[raiz].items()]
        heapq.heapify(heap)
        while heap:
            w, u, v = heapq.heappop(heap)
            if v in visitados:
                continue
            visitados.add(v)
            arvore.append((u, v, w))
            for x, wx in g.adj[v].items():
                if x not in visitados:
                    heapq.heappush(heap, (wx, v, x))
    return arvore


class UnionFind:
    """Union-find com compressão de caminho e união por posto."""

    def __init__(self, elementos):
        self.pai = {e: e for e in elementos}
        self.posto = {e: 0 for e in elementos}

    def achar(self, x):
        raiz = x
        while self.pai[raiz] != raiz:
            raiz = self.pai[raiz]
        while self.pai[x] != raiz:
            self.pai[x], x = raiz, self.pai[x]
        return raiz

    def unir(self, a, b):
        ra, rb = self.achar(a), self.achar(b)
        if ra == rb:
            return False
        if self.posto[ra] < self.posto[rb]:
            ra, rb = rb, ra
        self.pai[rb] = ra
        if self.posto[ra] == self.posto[rb]:
            self.posto[ra] += 1
        return True


def kruskal(g):
    """Floresta geradora mínima (Kruskal + union-find). O(E log E) = O(E log V)."""
    uf = UnionFind(g.adj)
    arvore = []
    for u, v, w in sorted(g.arestas(), key=lambda a: a[2]):
        if uf.unir(u, v):
            arvore.append((u, v, w))
    return arvore


def bfs_completa(g):
    """BFS reiniciada em cada vértice ainda não visitado (cobre todas as componentes)."""
    visitados = {}
    for s in g.adj:
        if s not in visitados:
            visitados.update(bfs(g, s))
    return visitados


def dfs_completa(g):
    """DFS reiniciada em cada vértice ainda não visitado (cobre todas as componentes)."""
    visitados = set()
    ordem = []
    for s in g.adj:
        if s not in visitados:
            parcial = dfs(g, s)
            visitados.update(parcial)
            ordem.extend(parcial)
    return ordem


def dijkstra_todos(g):
    """Dijkstra a partir de cada vértice (caminhos mínimos entre todos os pares)."""
    return {s: dijkstra(g, s) for s in g.adj}


def bellman_ford_todos(g):
    """Bellman-Ford a partir de cada vértice."""
    return {s: bellman_ford(g, s)[0] for s in g.adj}


def pontes_articulacoes(g):
    """Pontes e pontos de articulação de um grafo não-direcionado (low-link de Tarjan),
    versão iterativa. Retorna (pontes, articulações). O(V + E)."""
    descoberta = {}
    baixo = {}
    pai = {}
    pontes = []
    articulacoes = set()
    tempo = 0
    for raiz in g.adj:
        if raiz in descoberta:
            continue
        descoberta[raiz] = baixo[raiz] = tempo
        tempo += 1
        pai[raiz] = None
        filhos_raiz = 0
        pilha = [(raiz, iter(g.adj[raiz]))]
        while pilha:
            u, it = pilha[-1]
            avancou = False
            for v in it:
                if v not in descoberta:
                    pai[v] = u
                    descoberta[v] = baixo[v] = tempo
                    tempo += 1
                    if u == raiz:
                        filhos_raiz += 1
                    pilha.append((v, iter(g.adj[v])))
                    avancou = True
                    break
                if v != pai[u]:
                    baixo[u] = min(baixo[u], descoberta[v])
            if avancou:
                continue
            pilha.pop()
            p = pai[u]
            if p is not None:
                baixo[p] = min(baixo[p], baixo[u])
                if baixo[u] > descoberta[p]:
                    pontes.append((p, u))
                if p != raiz and baixo[u] >= descoberta[p]:
                    articulacoes.add(p)
        if filhos_raiz > 1:
            articulacoes.add(raiz)
    return pontes, sorted(articulacoes)
