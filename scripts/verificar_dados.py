"""Confere os dados brutos baixados com as estatísticas publicadas no SNAP.

Uso: python3 scripts/verificar_dados.py   (só usa a biblioteca padrão)
"""

import csv
import os
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORIGEM = os.path.join(RAIZ, "data", "raw", "comm-f2f-Resistance")
ESPERADO = {"partidas": 62, "participantes": 451, "arestas_temporais": 3_126_993}


def contar():
    with open(os.path.join(ORIGEM, "network_list.csv"), newline="") as f:
        partidas = [(int(l["NETWORK"]), int(l["NUMBER_OF_PARTICIPANTS"])) for l in csv.DictReader(f)]
    arestas = 0
    for partida, _ in partidas:
        with open(os.path.join(ORIGEM, "network", f"network{partida}.csv"), newline="") as f:
            leitor = csv.reader(f)
            next(leitor)
            for linha in leitor:
                arestas += sum(int(float(x)) for x in linha[1:])
    return {"partidas": len(partidas), "participantes": sum(n for _, n in partidas),
            "arestas_temporais": arestas}


def main():
    obtido = contar()
    ok = True
    for chave, valor in ESPERADO.items():
        confere = obtido[chave] == valor
        ok &= confere
        print(f"[{'OK' if confere else 'ERRO'}] {chave}: {obtido[chave]:,} (SNAP: {valor:,})")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
