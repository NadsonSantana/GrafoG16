"""Executa todas as análises (Partes I, II e III) e gera as tabelas do relatório.

Uso: python -m src.executar [--rapido]
Cada parte também pode ser executada isoladamente: python -m src.parte1 | parte2 | parte3
"""

from . import parte1, parte2, parte3, tabelas
from .comum import RESULTADOS, FIGURAS, carregar, preparar


def main():
    rapido = preparar()
    g, d = carregar()
    print("Parte I - métricas estruturais")
    parte1.parte1(g)
    print("Parte II - algoritmos")
    parte2.parte2(g, d, rapido)
    parte2.parte2_escala(g, d, rapido)
    print("Parte III - propriedades da rede")
    parte3.parte3(g, d, rapido)
    tabelas.gerar()
    print(f"Resultados e tabelas em {RESULTADOS}; figuras em {FIGURAS}")


if __name__ == "__main__":
    main()
