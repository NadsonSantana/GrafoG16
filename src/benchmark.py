"""Medição de tempo real com média, desvio padrão e intervalo de confiança de 95%."""

import math
import statistics
import time

from scipy import stats


def intervalo_confianca(amostras, confianca=0.95):
    """IC da média: Normal (Z) se n >= 30; t-Student com n-1 g.l. se n < 30."""
    n = len(amostras)
    media = statistics.mean(amostras)
    dp = statistics.stdev(amostras)
    alfa = 1 - confianca
    if n >= 30:
        critico, dist = stats.norm.ppf(1 - alfa / 2), "Z"
    else:
        critico, dist = stats.t.ppf(1 - alfa / 2, n - 1), "t"
    margem = critico * dp / math.sqrt(n)
    return {"n": n, "media": media, "desvio": dp, "ic_inf": media - margem,
            "ic_sup": media + margem, "distribuicao": dist, "valor_critico": float(critico)}


def medir(funcao, repeticoes, aquecimento=1):
    """Executa a função 'repeticoes' vezes (após aquecimento) e devolve (resultado, estatísticas)."""
    resultado = None
    for _ in range(aquecimento):
        resultado = funcao()
    tempos = []
    for _ in range(repeticoes):
        inicio = time.perf_counter()
        resultado = funcao()
        tempos.append(time.perf_counter() - inicio)
    return resultado, intervalo_confianca(tempos)
