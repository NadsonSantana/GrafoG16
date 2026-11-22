# Analisando Grafos do SNAP — G16 `comm-f2f-Resistance`

Trabalho prático de Teoria dos Grafos (UFBA — Instituto de Computação).
**Autor:** Nadson de Santana Sousa — Matrícula 224215945

## Determinação do grafo

Dada a matrícula `M = 224215945`:

| Passo | Valor |
|---|---|
| Soma dos dígitos `S` | 2+2+4+2+1+5+9+4+5 = **34** |
| Dois últimos dígitos `D` | **45** |
| `(S × (D + 1)) mod 129` | (34 × 46) mod 129 = 1564 mod 129 = **16** |

Grafo atribuído: **G16 — [comm-f2f-Resistance](https://snap.stanford.edu/data/comm-f2f-Resistance.html)**
(redes dinâmicas de interação face a face em 62 partidas do jogo *The Resistance*).

## Estrutura

```
data/processed/        dados tratados (CSV) usados nas análises
  vertices.csv           451 participantes (id, partida, participante)
  arestas_dirigidas.csv  1.241 arcos u->v (atenção e distância)
  arestas.csv            951 arestas não-direcionadas (atenção e distância)
  partidas.csv           metadados das 62 partidas
  etapas_tratamento.json evolução do grafo ao longo do tratamento
src/
  tratamento.py        dados brutos -> grafo estático tratado
  grafo.py             lista de adjacência
  algoritmos.py        BFS, DFS, Eulerianidade, Dijkstra, Bellman-Ford, Floyd-Warshall,
                       Tarjan (SCC e pontes/articulações), Prim, Kruskal
  metricas.py          métricas estruturais
  propriedades.py      small-world, lei de potência, centralidades, robustez e modelos
  benchmark.py         tempo real: média, desvio e IC 95% (Z ou t-Student)
  comum.py             caminhos, estilo dos gráficos e carga do grafo tratado
  parte1.py            Parte I  (python -m src.parte1)
  parte2.py            Parte II (python -m src.parte2), inclui o experimento de escala
  parte3.py            Parte III (python -m src.parte3)
  tabelas.py           tabelas LaTeX e macros do relatório a partir dos JSONs
  executar.py          executa as três partes e gera as tabelas
tests/                 testes (NetworkX usado apenas como oráculo)
results/               resultados em JSON, tabelas tab_*.tex e macros valores.tex
figures/               figuras do relatório
relatorio/             relatório em LaTeX (template WebMedia) e PDF
scripts/baixar_dados.sh  download do dataset original
```

## Tratamento realizado

1. **Temporal → estático:** para cada par ordenado (u, v) da mesma partida calcula-se
   a atenção `f(u, v)`, a fração dos instantes (1/3 s) em que u olha para v.
2. **Remoção do alvo LAPTOP** (não é participante) e de auto-laços.
3. **Limiar de atenção:** mantém-se u → v apenas se `f(u, v) ≥ 1/N`, sendo N o número de alvos
   possíveis de u (N−1 colegas + laptop), ou seja, se u olha para v mais do que sob atenção uniforme.
4. **Dirigido → não-direcionado:** {u, v} existe se u → v ou v → u; peso
   `w = (f(u,v) + f(v,u)) / 2` e distância `d = −ln(w)` (usada nos algoritmos ponderados).
5. **União das 62 partidas** em um único grafo (vértice `gXX_pY`). A versão dirigida é usada no Tarjan.

## Como executar

Requisitos: Python 3.10+.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
./scripts/baixar_dados.sh            # baixa ~44 MB do SNAP para data/raw/
.venv/bin/python -m src.tratamento   # regenera data/processed/
.venv/bin/python -m pytest           # testes
.venv/bin/python -m src.executar     # análises completas (~1 min); --rapido para teste
.venv/bin/python -m src.parte2       # ou apenas uma parte (parte1, parte2, parte3)
.venv/bin/python -m src.tabelas      # só regenera as tabelas do relatório a partir dos JSONs
```

O relatório inclui as tabelas de `results/` com `\input` e os valores citados no texto
por macros (`results/valores.tex`). Assim, após uma nova execução, basta recompilar o PDF.

Os dados tratados já estão versionados; o download só é necessário para refazer o tratamento.

Para compilar o relatório:

```bash
cd relatorio && pdflatex main.tex && pdflatex main.tex
```

## Onde mexer

| Quero alterar... | Arquivo |
|---|---|
| leitura dos dados brutos, limiar de atenção, pesos | `src/tratamento.py` |
| estrutura de grafo (lista de adjacência) | `src/grafo.py` |
| um algoritmo da Parte II | `src/algoritmos.py` |
| repetições, tamanhos do experimento de escala | `src/parte2.py` (`parte2`, `parte2_escala`) |
| uma métrica da Parte I | `src/metricas.py` |
| small-world, lei de potência, centralidades, robustez, modelos | `src/propriedades.py` |
| cálculo do IC 95% | `src/benchmark.py` |
| aparência de uma figura | função `figura_*` da respectiva `src/parteN.py`; cores em `src/comum.py` |
| formato de uma tabela ou valor citado no relatório | `src/tabelas.py` |
| texto do relatório | `relatorio/main.tex` |

## Referências

- S. Kumar, C. Bai, V. S. Subrahmanian, J. Leskovec. *Deception Detection in Group Video
  Conversations using Dynamic Interaction Networks*. ICWSM, 2021.
- C. Bai et al. *Predicting the Visual Focus of Attention in Multi-Person Discussion Videos*. IJCAI, 2019.
- J. Leskovec, A. Krevl. *SNAP Datasets: Stanford Large Network Dataset Collection*, 2014.
