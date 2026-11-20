"""Etapa 5b - Parte III, questao 2: a distribuicao de graus segue lei de potencia?

Procedimento (Apendice B + Clauset, Shalizi & Newman, 2009), implementado aqui:
  1. Distribuicao de graus, histograma e CCDF em escala log-log.
  2. Regressao linear em log-log (P(k) e CCDF) - o teste "visual" do enunciado.
  3. Ajuste por maxima verossimilhanca da lei de potencia discreta
     p(k) = k^-alpha / zeta(alpha, k_min), com k_min escolhido pela menor
     distancia de Kolmogorov-Smirnov (KS).
  4. p-valor de aderencia por bootstrap semiparametrico.
  5. Razao de verossimilhanca (teste de Vuong) contra alternativas com a mesma
     cauda k >= k_min: exponencial discreta e Poisson truncada (modelo ER).
"""
import math
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, special, stats

from comum import (AMARELO, AQUA, AZUL, DIR_FIG, DIR_PROC, LARANJA, TEXTO, TEXTO_2,
                   carregar_grafo, estilo)

SEMENTE = 224215945
BOOTSTRAPS = 1000


# ------------------------------------------------------- lei de potencia
def ll_potencia(alpha, x, kmin):
    return -alpha * np.log(x).sum() - len(x) * np.log(special.zeta(alpha, kmin))


def ajustar_potencia(x, kmin):
    res = optimize.minimize_scalar(lambda a: -ll_potencia(a, x, kmin),
                                   bounds=(1.01, 15), method="bounded")
    return res.x


def cdf_potencia(k, alpha, kmin):
    """P(X <= k) para a lei de potencia discreta com suporte k >= kmin."""
    return 1 - special.zeta(alpha, k + 1) / special.zeta(alpha, kmin)


def ks(x, cdf_teorica):
    ks_ = np.arange(x.min(), x.max() + 1)
    emp = np.array([(x <= k).mean() for k in ks_])
    return np.abs(emp - cdf_teorica(ks_)).max()


def melhor_kmin(graus):
    melhor = None
    for kmin in sorted(set(graus[graus >= 1])):
        x = graus[graus >= kmin]
        if len(np.unique(x)) < 3:
            break
        a = ajustar_potencia(x, kmin)
        d = ks(x, lambda k: cdf_potencia(k, a, kmin))
        if melhor is None or d < melhor[2]:
            melhor = (kmin, a, d)
    return melhor


def amostrar_potencia(alpha, kmin, n, rng, kmax=10_000):
    ks_ = np.arange(kmin, kmax)
    p = ks_ ** -alpha
    return rng.choice(ks_, size=n, p=p / p.sum())


def bootstrap_pvalor(graus, kmin, alpha, d_obs, rng):
    """Bootstrap semiparametrico (Clauset et al.): reamostra o corpo k < kmin,
    gera a cauda pelo modelo ajustado e refaz todo o ajuste (inclusive k_min)."""
    n = len(graus)
    corpo = graus[(graus >= 1) & (graus < kmin)]
    n_cauda = int((graus >= kmin).sum())
    maiores = 0
    for _ in range(BOOTSTRAPS):
        nc = rng.binomial(n, n_cauda / n)
        sint = np.concatenate([rng.choice(corpo, n - nc) if len(corpo) else np.array([], int),
                               amostrar_potencia(alpha, kmin, nc, rng)])
        r = melhor_kmin(sint)
        if r is not None and r[2] >= d_obs:
            maiores += 1
    return maiores / BOOTSTRAPS


# ------------------------------------------------------- alternativas
def ll_exponencial_pontos(x, kmin):
    """Exponencial (geometrica) discreta truncada em kmin: p(k) = (1-e^-l) e^-l(k-kmin)."""
    lam = optimize.minimize_scalar(
        lambda l: -np.sum(np.log1p(-np.exp(-l)) - l * (x - kmin)),
        bounds=(1e-6, 20), method="bounded").x
    return np.log1p(-np.exp(-lam)) - lam * (x - kmin), lam


def ll_poisson_pontos(x, kmin):
    """Poisson truncada (k >= kmin): p(k) = Pois(k; mu) / P(K >= kmin)."""
    def ll(mu):
        return stats.poisson.logpmf(x, mu) - np.log(stats.poisson.sf(kmin - 1, mu))
    mu = optimize.minimize_scalar(lambda m: -ll(m).sum(), bounds=(1e-6, 50),
                                  method="bounded").x
    return ll(mu), mu


def vuong(ll1, ll2):
    """R = soma das diferencas de log-verossimilhanca; p bilateral (Vuong, 1989).
    R > 0 favorece o modelo 1 (lei de potencia)."""
    d = ll1 - ll2
    R = d.sum()
    sd = d.std(ddof=0)
    p = math.erfc(abs(R) / (math.sqrt(2 * len(d)) * sd)) if sd > 0 else float("nan")
    return R, p


def regressao(x, y):
    r = stats.linregress(x, y)
    return r.slope, r.intercept, r.rvalue ** 2


# ------------------------------------------------------- figura
def figura(graus, kmin, alpha, mu_ger, lam, mu_trunc):
    n = len(graus)
    pos = graus[graus >= 1]
    cont = Counter(pos)
    ks_ = np.array(sorted(cont))
    pk = np.array([cont[k] / n for k in ks_])
    ccdf = np.array([(graus >= k).mean() for k in ks_])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), dpi=150)
    ax = axes[0]
    ax.scatter(ks_, pk, s=40, color=AZUL, zorder=3, label="observado")
    s, b, r2 = regressao(np.log10(ks_), np.log10(pk))
    xx = np.linspace(ks_.min(), ks_.max(), 50)
    ax.plot(xx, 10 ** b * xx ** s, color=LARANJA, linewidth=2,
            label=f"reta log-log: inclinação {s:.2f}, R² = {r2:.2f}")
    ax.plot(ks_, stats.poisson.pmf(ks_, mu_ger), color=AQUA, linewidth=2, marker="s", ms=4,
            label=f"Poisson (<k> = {mu_ger:.2f}, modelo ER)")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xticks(ks_), ax.set_xticklabels(ks_)
    ax.set_xlabel("grau k (escala log)", color=TEXTO)
    ax.set_ylabel("P(k) (escala log)", color=TEXTO)
    ax.set_title("Distribuição de graus em log-log (k ≥ 1)", color=TEXTO, fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    estilo(ax, grade="both")

    ax = axes[1]
    ax.scatter(ks_, ccdf, s=40, color=AZUL, zorder=3, label="observado P(K ≥ k)")
    cauda = ks_[ks_ >= kmin]
    frac = (graus >= kmin).mean()
    ax.plot(cauda, frac * (1 - cdf_potencia(cauda - 1, alpha, kmin)), color=LARANJA,
            linewidth=2, label=f"lei de potência MLE (α = {alpha:.2f}, k_min = {kmin})")
    ax.plot(cauda, frac * np.exp(-lam * (cauda - kmin)), color=AMARELO, linewidth=2,
            linestyle="--", label="exponencial (k ≥ k_min)")
    ax.plot(cauda, frac * stats.poisson.sf(cauda - 1, mu_trunc)
            / stats.poisson.sf(kmin - 1, mu_trunc), color=AQUA, linewidth=2, linestyle=":",
            label="Poisson truncada (k ≥ k_min)")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xticks(ks_), ax.set_xticklabels(ks_)
    ax.set_xlabel("grau k (escala log)", color=TEXTO)
    ax.set_ylabel("CCDF  P(K ≥ k) (escala log)", color=TEXTO)
    ax.set_title("Distribuição acumulada complementar (CCDF)", color=TEXTO, fontsize=10,
                 loc="left")
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    estilo(ax, grade="both")
    fig.tight_layout()
    fig.savefig(DIR_FIG / "parte3_lei_potencia.png")
    plt.close(fig)


def main():
    rng = np.random.default_rng(SEMENTE)
    G = carregar_grafo()
    graus = np.array([d for _, d in G.degree()])
    n = len(graus)
    pos = graus[graus >= 1]

    # 2. Regressoes em log-log (k >= 1) e log-linear
    cont = Counter(pos)
    ks_ = np.array(sorted(cont))
    pk = np.array([cont[k] / n for k in ks_])
    ccdf = np.array([(graus >= k).mean() for k in ks_])
    s_pk, _, r2_pk = regressao(np.log10(ks_), np.log10(pk))
    s_cc, _, r2_cc = regressao(np.log10(ks_), np.log10(ccdf))
    _, _, r2_cc_exp = regressao(ks_, np.log10(ccdf))
    # P(k) e crescente ate k = 2: regressao so na parte decrescente
    dec = ks_ >= ks_[np.argmax(pk)]
    s_dec, _, r2_dec = regressao(np.log10(ks_[dec]), np.log10(pk[dec]))
    _, _, r2_dec_exp = regressao(ks_[dec], np.log10(pk[dec]))

    # 3. MLE + k_min por KS
    kmin, alpha, d_ks = melhor_kmin(graus)
    x = graus[graus >= kmin]
    erro_alpha = (alpha - 1) / math.sqrt(len(x))

    # 4. p-valor por bootstrap
    p_boot = bootstrap_pvalor(graus, kmin, alpha, d_ks, rng)

    # 5. Alternativas na mesma cauda
    ll_pot = -alpha * np.log(x) - np.log(special.zeta(alpha, kmin))
    ll_exp, lam = ll_exponencial_pontos(x, kmin)
    ll_poi, mu_trunc = ll_poisson_pontos(x, kmin)
    R_exp, p_exp = vuong(ll_pot, ll_exp)
    R_poi, p_poi = vuong(ll_pot, ll_poi)

    resumo = {
        "vertices": n, "grau_min": int(graus.min()), "grau_max": int(graus.max()),
        "grau_medio": graus.mean(), "variancia_grau": graus.var(ddof=1),
        "razao_variancia_media": graus.var(ddof=1) / graus.mean(),
        "decadas_cobertas_k>=1": math.log10(pos.max() / pos.min()),
        "inclinacao_loglog_Pk": s_pk, "R2_loglog_Pk": r2_pk,
        "inclinacao_loglog_Pk_decrescente": s_dec, "R2_loglog_Pk_decrescente": r2_dec,
        "R2_loglinear_Pk_decrescente": r2_dec_exp,
        "inclinacao_loglog_CCDF": s_cc, "R2_loglog_CCDF": r2_cc,
        "R2_loglinear_CCDF": r2_cc_exp,
        "k_min": kmin, "alpha_MLE": alpha, "erro_padrao_alpha": erro_alpha,
        "n_cauda": len(x), "KS": d_ks, "p_bootstrap": p_boot,
        "lambda_exponencial": lam, "mu_poisson_truncada": mu_trunc,
        "R_vuong_pot_vs_exp": R_exp, "p_vuong_pot_vs_exp": p_exp,
        "R_vuong_pot_vs_poisson": R_poi, "p_vuong_pot_vs_poisson": p_poi,
    }
    pd.DataFrame(list(resumo.items()), columns=["medida", "valor"]) \
        .to_csv(DIR_PROC / "parte3_lei_potencia.csv", index=False)
    figura(graus, kmin, alpha, graus.mean(), lam, mu_trunc)

    for k, v in resumo.items():
        print(f"  {k:36s} {v:.4f}" if isinstance(v, float) else f"  {k:36s} {v}")


if __name__ == "__main__":
    main()
