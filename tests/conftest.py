import pytest

from src.grafo import Grafo
from tests.auxiliares import grafo_aleatorio


@pytest.fixture(params=[1, 2, 3])
def grafo(request):
    """Grafo aleatório ponderado e desconexo (três sementes)."""
    return grafo_aleatorio(40, 0.15, request.param)


@pytest.fixture
def real():
    """Grafo tratado do trabalho (data/processed)."""
    return Grafo.de_csv("data/processed/arestas.csv", "data/processed/vertices.csv")
