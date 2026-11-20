"""Executa todas as etapas em sequencia (dados brutos -> tabelas e figuras)."""
import subprocess
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
SCRIPTS = [
    "01_inspecao.py",
    "02_tratamento.py",
    "03_parte1_estrutural.py",
    "04_parte2_algoritmos.py",
    "05a_small_world.py",
    "05b_lei_potencia.py",
    "05c_robustez.py",
    "05d_descoberta.py",
    "05e_modelos.py",
]

for s in SCRIPTS:
    print(f"\n===== {s} =====", flush=True)
    subprocess.run([sys.executable, str(RAIZ / "src" / s)], check=True)
