"""Intervalo de confiança e medição de tempo."""

import pytest

from src.benchmark import intervalo_confianca, medir


def test_intervalo_confianca_t():
    ic = intervalo_confianca([1.0, 2.0, 3.0, 4.0, 5.0])
    assert ic["distribuicao"] == "t" and ic["media"] == 3.0
    assert ic["valor_critico"] == pytest.approx(2.776, abs=1e-3)
    assert ic["ic_inf"] < 3.0 < ic["ic_sup"]


def test_intervalo_confianca_z():
    ic = intervalo_confianca([1.0, 2.0] * 15)
    assert ic["distribuicao"] == "Z" and abs(ic["valor_critico"] - 1.96) < 0.001


def test_medir():
    resultado, est = medir(lambda: 42, 5)
    assert resultado == 42 and est["n"] == 5 and est["mediana"] >= 0

