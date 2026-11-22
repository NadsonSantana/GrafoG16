#!/usr/bin/env bash
# Baixa e extrai o dataset comm-f2f-Resistance (SNAP) em data/raw/.
set -euo pipefail

RAIZ="$(cd "$(dirname "$0")/.." && pwd)"
DESTINO="$RAIZ/data/raw"
URL="https://snap.stanford.edu/data/comm-f2f-Resistance-network.tar.gz"

mkdir -p "$DESTINO"
if [ ! -d "$DESTINO/comm-f2f-Resistance" ]; then
    curl -L "$URL" -o "$DESTINO/comm-f2f-Resistance-network.tar.gz"
    tar xzf "$DESTINO/comm-f2f-Resistance-network.tar.gz" -C "$DESTINO"
fi
python3 "$RAIZ/scripts/verificar_dados.py"
echo "Dados disponíveis em $DESTINO/comm-f2f-Resistance"
