"""Métricas estruturais (Parte I). Cálculo exato: o grafo tem apenas 451 vértices."""

from collections import Counter

from .algoritmos import bfs, componentes_conexas


def graus(g):
    return {v: len(g.adj[v]) for v in g.adj}


def distribuicao_graus(g):
    """Retorna {k: número de vértices com grau k}, ordenado por k."""
    return dict(sorted(Counter(graus(g).values()).items()))


def densidade(g):
    n, m = g.num_vertices(), g.num_arestas()
    if n < 2:
        return 0.0
    return m / (n * (n - 1)) if g.dirigido else 2 * m / (n * (n - 1))


def excentricidades(g, vertices=None):
    """Excentricidade (em saltos) de cada vértice dentro da sua componente."""
    vs = vertices if vertices is not None else g.adj
    return {v: max(bfs(g, v).values()) for v in vs}


def caminhos_medios(g):
    """Comprimento médio dos caminhos mínimos sobre todos os pares alcançáveis (u != v).

    Em um grafo desconexo, pares de componentes distintas são ignorados (distância infinita).
    """
    soma = 0
    pares = 0
    for v in g.adj:
        for u, d in bfs(g, v).items():
            if u != v:
                soma += d
                pares += 1
    return soma / pares if pares else 0.0


def triangulos_por_vertice(g):
    """Número de triângulos que contêm cada vértice (interseção de vizinhanças)."""
    t = {}
    for v in g.adj:
        viz = g.adj[v]
        cont = 0
        for u in viz:
            cont += sum(1 for w in g.adj[u] if w in viz)
        t[v] = cont // 2
    return t


def numero_triangulos(g):
    return sum(triangulos_por_vertice(g).values()) // 3


def clusterizacao_local(g):
    """C_i = 2 T_i / (k_i (k_i - 1)); vale 0 quando k_i < 2."""
    t = triangulos_por_vertice(g)
    c = {}
    for v in g.adj:
        k = len(g.adj[v])
        c[v] = 2 * t[v] / (k * (k - 1)) if k >= 2 else 0.0
    return c


def clusterizacao_media(g):
    c = clusterizacao_local(g)
    return sum(c.values()) / len(c) if c else 0.0


def metricas_estruturais(g):
    """Calcula todas as métricas da Parte I e devolve um dicionário."""
    gr = graus(g)
    comps = componentes_conexas(g)
    tamanhos = [len(c) for c in comps]
    lcc = g.subgrafo(comps[0])
    exc_lcc = excentricidades(lcc)
    exc_todos = excentricidades(g)
    # vértices isolados têm excentricidade 0 e não definem raio de forma útil
    exc_nao_triviais = [e for v, e in exc_todos.items() if gr[v] > 0]
    return {
        "vertices": g.num_vertices(),
        "arestas": g.num_arestas(),
        "grau_min": min(gr.values()),
        "grau_max": max(gr.values()),
        "grau_medio": sum(gr.values()) / len(gr),
        "densidade": densidade(g),
        "num_componentes": len(comps),
        "tamanhos_componentes": tamanhos,
        "tamanho_lcc": tamanhos[0],
        "vertices_isolados": sum(1 for k in gr.values() if k == 0),
        "diametro_lcc": max(exc_lcc.values()),
        "raio_lcc": min(exc_lcc.values()),
        "caminho_medio_lcc": caminhos_medios(lcc),
        "diametro_max_componentes": max(exc_todos.values()),
        "raio_min_componentes": min(exc_nao_triviais) if exc_nao_triviais else 0,
        "caminho_medio_global": caminhos_medios(g),
        "clusterizacao_media": clusterizacao_media(g),
        "triangulos": numero_triangulos(g),
        "distribuicao_graus": distribuicao_graus(g),
    }
