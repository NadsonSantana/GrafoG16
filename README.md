# GrafoG16 – Analisando Grafos do SNAP

Trabalho prático da disciplina de Grafos (UFBA/IC/DCC): análise do grafo **G16**,
[comm-f2f-Resistance](https://snap.stanford.edu/data/comm-f2f-Resistance.html) do SNAP.

## Determinação do grafo

Matrícula `224215945`:

- soma dos dígitos = 2+2+4+2+1+5+9+4+5 = **34**
- últimos dois dígitos = **45**
- (34 × (45 + 1)) mod 129 = 1564 mod 129 = **16** → **G16**

## O grafo

Redes dinâmicas de interação face a face entre jogadores de *The Resistance*: 62 partidas,
5 a 8 jogadores por partida, 451 vértices e 3.126.993 arestas temporais (ponderadas e dirigidas).
A cada 1/3 s, a aresta u → v indica que o participante u olha para v (ou para o laptop).

Referência: S. Kumar, C. Bai, V. S. Subrahmanian, J. Leskovec. *Deception Detection in Group Video
Conversations using Dynamic Interaction Networks*. ICWSM, 2021.

## Como executar

Requisitos: Python 3.11+.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Baixe os dados brutos (não versionados, ~44 MB):

```bash
mkdir -p data/raw
curl -L -o data/raw/comm-f2f-Resistance-network.tar.gz \
  https://snap.stanford.edu/data/comm-f2f-Resistance-network.tar.gz
tar -xzf data/raw/comm-f2f-Resistance-network.tar.gz -C data/raw
```

SHA-256 do arquivo usado: `1155481efbd883804aabf4f31bb43188f0ab5875f70dd6a16a1177d81e7d10af`.

Execute tudo (cerca de 3 min):

```bash
.venv/bin/python executar_tudo.py
```

Ou cada etapa separadamente:

| Etapa | Comando | Saídas |
|---|---|---|
| 1. Inspeção dos dados brutos | `.venv/bin/python src/01_inspecao.py` | `data/processed/inspecao_bruta.csv` |
| 2. Tratamento | `.venv/bin/python src/02_tratamento.py` | `data/processed/g16_*.csv`, `pares_atencao.csv`, `sensibilidade_limiar.csv`, `figures/distribuicao_atencao.png` |
| 3. Parte I – análise estrutural | `.venv/bin/python src/03_parte1_estrutural.py` | `parte1_*.csv`, `figures/parte1_*.png` |
| 4. Parte II – algoritmos | `.venv/bin/python src/04_parte2_algoritmos.py` | `parte2_*.csv`, `figures/parte2_*.png` |
| 5a. Parte III – small-world | `.venv/bin/python src/05a_small_world.py` | `parte3_small_world_*.csv`, `figures/parte3_small_world.png` |
| 5b. Parte III – lei de potência | `.venv/bin/python src/05b_lei_potencia.py` | `parte3_lei_potencia.csv`, `figures/parte3_lei_potencia.png` |
| 5c. Parte III – robustez | `.venv/bin/python src/05c_robustez.py` | `parte3_robustez_*.csv`, `figures/parte3_robustez_*.png` |
| 5d. Parte III – descoberta | `.venv/bin/python src/05d_descoberta.py` | `parte3_descoberta_*.csv`, `figures/parte3_descoberta.png` |
| 5e. Parte III – extra (ER/BA/WS) | `.venv/bin/python src/05e_modelos.py` | `parte3_modelos.csv`, `figures/parte3_modelos.png` |

Etapas 3 a 5 dependem apenas de `data/processed/g16_*.csv` (versionados); a 5d também lê os dados brutos.
Todos os sorteios usam a semente `224215945`. Os tempos da Parte II dependem da máquina e da carga:
para medições estáveis, execute com a máquina ociosa.

## Tratamento dos dados

1. **Temporal → estático:** para cada partida e par dirigido (u, v), atenção `a_uv` = fração dos
   quadros em que u olha para v (versão binária; a probabilidade média da versão ponderada fica no atributo `prob`).
2. **Remoções:** arestas para o `LAPTOP` (não é participante) e auto-loops (inexistentes).
3. **Limiar relativo:** a ligação é mantida se a atenção ≥ `1/(n−1)`, isto é, acima da atenção uniforme
   entre os n−1 colegas da partida (sem limiar, 2.867 dos 2.870 pares dirigidos possíveis estariam presentes).
4. **Dirigido → não dirigido:** atenção mútua `m_uv = (a_uv + a_vu)/2` com o mesmo limiar;
   peso `distancia = 1 − m_uv` (usado em caminhos mínimos e árvore geradora).
5. **Vértices:** rótulo `R<partida>_P<participante>`; os 451 são mantidos, inclusive isolados.

| Arquivo | Conteúdo |
|---|---|
| `data/processed/g16_nos.csv` | 451 vértices (partida, participante) |
| `data/processed/g16_nao_dirigido.csv` | 548 arestas: `atencao`, `prob`, `distancia` |
| `data/processed/g16_dirigido.csv` | 1.045 arcos (usado no Tarjan de componentes fortemente conexas e na descoberta) |

O grafo é **desconexo por construção** (partidas não compartilham jogadores): 85 componentes,
a maior com 8 vértices. Por isso, medidas de distância são calculadas por componente e a robustez é
avaliada no grafo inteiro.

## Principais resultados

| Pergunta | Resposta |
|---|---|
| Small-world? | Indícios fracos: C/C_rand = 1,32 e L/L_rand = 0,95 por componente (σ mediano 1,56), mas C só é significativamente maior em 9 de 60 componentes |
| Lei de potência? | Não: graus de 0 a 7, variância ≈ média (Poisson); bootstrap rejeita a lei de potência (p = 0,00) |
| Robustez a 5% aleatório? | Sim: 88% dos pares conectados preservados |
| Robustez a 5% mais centrais? | Não: 73% (grau) e 67% (betweenness), 13 a 18 desvios padrão abaixo do aleatório |
| Descoberta | "Jogador focal": em 32/62 partidas todos olham um mesmo jogador acima do uniforme, e em 42/62 o mais olhado é o mesmo nos três terços da partida (1,2 esperado ao acaso) |
| Modelo clássico mais próximo (extra) | Barabási–Albert (distância de variação total dos graus = 0,04) |

## Estrutura

```
executar_tudo.py        executa todas as etapas
requirements.txt        dependências com versões fixas
src/
  algoritmos.py         implementações próprias (BFS, DFS, Euler, Dijkstra, Bellman-Ford,
                        Floyd-Warshall, Tarjan, Kruskal, Prim)
  comum.py              utilitários da Parte III
  01_ ... 05e_*.py      um script por etapa
data/processed/         dados tratados e tabelas de resultados
figures/                gráficos usados no relatório
relatorio/              relatório em LaTeX (formato WebMedia/ACM)
```
