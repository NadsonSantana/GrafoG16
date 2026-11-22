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

> O mesmo dataset aparece em quatro entradas da planilha ANALISANDO-SNAP-DATASET
> (G16, G37, G103 e G121). A regra da matrícula aponta para a entrada **G16**, a adotada aqui.

## Checklist da especificação

| Requisito | Onde está |
|---|---|
| Cálculo de determinação do grafo | Relatório, Seção 2 (Tabela 1); README, acima |
| Link para o repositório | Relatório: resumo (1ª página) e Seção 3 |
| Códigos e scripts | `src/`, `scripts/` |
| Conjunto de dados tratado | `data/processed/` |
| README explicando a execução | Este arquivo, "Como executar" |
| Descrição do grafo original | Relatório, Seção 4 (tabela do dataset e figura do dado bruto ao longo da partida) |
| Tratamento realizado | Relatório, Seção 5 (Tabela de evolução e figura do limiar de atenção); `src/tratamento.py` |
| Metodologia | Relatório, Seção 6 |
| Parte I: todas as medidas + tabela "computadas ou não" | Relatório, Seção 7.1 (tabela de métricas com colunas "Computada?" e "Método") |
| Parte I: visualização do grafo | Relatório, Figuras 1 e 3; `figures/grafo_completo.png` |
| Parte II: os 8 grupos de algoritmos | Relatório, Seção 7.2; `src/algoritmos.py` |
| Parte II: complexidade teórica e observada; média, desvio e IC (Z se n ≥ 30, t se n < 30) | Relatório, tabelas de tempos e de escala (V de 58 a ~29 mil, com réplicas das partidas); `src/benchmark.py`, `src/parte2.py` |
| Parte III.1: small-world (código próprio), faz sentido?, implicação prática | Relatório, Seção 7.3; `src/propriedades.py` |
| Parte III.2: lei de potência | Relatório, Seção 7.3 (regressão, CCDF, máxima verossimilhança e testes de aderência à Poisson e à binomial por partida) |
| Parte III.3: robustez aleatória e direcionada (boxplot, métricas A a D, r = ⌈0,05·n⌉) | Relatório, Seção 7.3 (figuras de robustez e da curva S(f) até 30%) |
| Parte III.4: descoberta mais interessante | Relatório, Seção 7.3 |
| Extra: modelos clássicos (ER, BA, WS) | Relatório, Seção 7.3 (tabela e CCDF dos graus) |
| Discussão crítica e conclusão | Relatório, Seções 8 e 9 (inclui a sensibilidade ao limiar de atenção) |
| Template WebMedia; PDF com o número da matrícula | `relatorio/main.tex`; `relatorio/224215945.pdf` |

## Estrutura

```
data/processed/        dados tratados (CSV) usados nas análises
  vertices.csv           451 participantes (id, partida, participante)
  arestas_dirigidas.csv  1.241 arcos u->v (atenção e distância)
  arestas.csv            951 arestas não-direcionadas (atenção e distância)
  partidas.csv           metadados das 62 partidas
  atencao_pares.csv      atenção f(u, v) de todos os 2.870 pares ordenados, antes do limiar
  serie_temporal.csv     laptop e foco da atenção ao longo de cada partida (20 faixas de tempo)
  etapas_tratamento.json evolução do grafo ao longo do tratamento
src/
  tratamento.py        dados brutos -> grafo estático tratado
  grafo.py             lista de adjacência
  algoritmos.py        BFS, DFS, Eulerianidade, Dijkstra, Bellman-Ford, Floyd-Warshall
                       (Python e numpy),
                       Tarjan (SCC e pontes/articulações), Prim, Kruskal
  metricas.py          métricas estruturais
  propriedades.py      small-world, lei de potência, centralidades, robustez e modelos
  sensibilidade.py     medidas principais refeitas para outros limiares de atenção
  benchmark.py         tempo real: média, desvio e IC 95% (Z ou t-Student)
  comum.py             caminhos, estilo dos gráficos e carga do grafo tratado
  parte1.py            Parte I  (python -m src.parte1)
  parte2.py            Parte II (python -m src.parte2), inclui o experimento de escala
  parte3.py            Parte III (python -m src.parte3)
  tabelas.py           tabelas LaTeX e macros do relatório a partir dos JSONs
  executar.py          executa as três partes e gera as tabelas
tests/                 testes por módulo (NetworkX usado apenas como oráculo)
results/               resultados em JSON, tabelas tab_*.tex e macros valores.tex
figures/               figuras em PDF (relatório) e PNG (README e notebook)
notebooks/analise.ipynb  reúne resultados, tabelas e figuras (não recalcula nada)
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

## Resultados principais

O relatório completo está em [`relatorio/224215945.pdf`](relatorio/224215945.pdf). Os números
abaixo vêm de `results/`; os tempos da Parte II variam a cada execução, os demais são
determinísticos (sementes fixas).

![Grafo completo: 62 partidas](figures/grafo_completo.png)

### Parte I: análise estrutural (tudo exato, sem amostragem)

| Medida | Valor |
|---|---:|
| Vértices / arestas | 451 / 951 |
| Grau mín. / máx. / médio / mediano | 1 / 7 / 4,22 / 4 |
| Densidade global / média intra-partida | 9,37×10⁻³ / 0,67 |
| Componentes (= partidas) | 62, de 5 a 8 vértices |
| Diâmetro / raio (maior componente) | 2 / 1 (diâmetro máximo entre componentes: 3) |
| Comprimento médio dos caminhos | 1,34 |
| Clusterização média / triângulos / transitividade | 0,754 / 791 / 0,676 |

### Parte II: algoritmos (grafo completo, IC 95%)

| Algoritmo | Complexidade | Média (ms) | IC 95% (ms) |
|---|---|---:|---|
| BFS | O(V+E) | 0,31 | [0,21; 0,42] |
| DFS | O(V+E) | 0,48 | [0,34; 0,62] |
| Eulerianidade | O(V+E) | 0,09 | [0,06; 0,12] |
| Dijkstra (todas as fontes) | O(V(V+E) log V) | 3,3 | [3,0; 3,5] |
| Bellman-Ford (todas as fontes) | O(V²E) | 384 | [328; 440] |
| Floyd-Warshall (Python, com poda) | O(V³) | 82 | [72; 93] (t-Student, n = 10) |
| Floyd-Warshall (numpy, sem poda) | O(V³) | 350 | [246; 453] (t-Student, n = 10) |
| Tarjan (SCC, dirigido) | O(V+E) | 0,71 | [0,69; 0,74] |
| Tarjan (pontes e articulações) | O(V+E) | 0,78 | [0,76; 0,81] |
| Prim | O(E log V) | 0,53 | [0,51; 0,55] |
| Kruskal | O(E log V) | 0,89 | [0,86; 0,93] |

- Dijkstra, Bellman-Ford e Floyd-Warshall produzem distâncias idênticas nos 2.870 pares alcançáveis.
- Floyd-Warshall vetorizado com numpy é ~4× **mais lento** que o puro: sem poda, faz sempre as
  V³ ≈ 9×10⁷ operações, enquanto a versão pura ignora intermediários inalcançáveis (cada vértice
  só alcança a própria partida).
- Tarjan: 131 SCCs (uma não-trivial por partida); apenas 3 pontes e 3 articulações em todo o grafo.
- Prim e Kruskal: floresta geradora mínima com 389 arestas (uma árvore por partida) e mesmo peso.
- **Complexidade observada** (ajuste log-log, V de 58 a 28.753): até 62 partidas, as partidas
  reais; acima disso, réplicas de partidas sorteadas com reposição, que estendem a faixa para
  mais de duas ordens de grandeza sem mudar o tipo de grafo. Expoente b entre 1,07 e 1,11 para
  os algoritmos lineares e de 0,97 a 1,00 para Prim e Kruskal; Dijkstra e Bellman-Ford de todas
  as fontes e Floyd-Warshall com poda crescem um grau de V abaixo do pior caso, porque cada fonte
  só alcança a própria partida. Medições com aquecimento e coletor de lixo desligado.

![Complexidade observada](figures/escala.png)

### Parte III: propriedades da rede

| Pergunta | Resposta |
|---|---|
| Small-world? | **Só aparentemente.** Contra um G(n,m) global, σ ≈ 301; contra um aleatório que respeita as partidas, σ ≈ 1,14 e L é o mesmo (1,34). Os caminhos curtos vêm do tamanho dos grupos (≤ 8 pessoas). |
| Lei de potência? | **Não.** A regressão log-log tem inclinação positiva (R² = 0,46); o ajuste por máxima verossimilhança discreto exato (zeta de Hurwitz) dá γ = 5,88 ± 0,38 numa cauda de só 3 valores de k, com KS = 0,13. A distribuição é subdispersa (Var(k) = 2,0 < ⟨k⟩ = 4,2): mais próxima de uma Poisson, mas os testes de dispersão e χ² rejeitam a Poisson (p < 10⁻¹⁵) e também a mistura de binomiais por partida (p ≈ 3×10⁻⁹); os graus são ainda mais homogêneos. |
| Robusta à remoção aleatória de 5%? | **Sim.** Com r = ⌈0,05·451⌉ = 23, 94,9% dos vértices continuam no núcleo conexo da sua partida (o máximo possível é 94,9%). |
| Robusta à remoção dos 5% mais centrais? | **Sim.** Por grau, betweenness e PageRank, a fração no núcleo cai só para 94,2%, 93,1% e 94,0%: não há hubs (k ≤ 7) e quase nenhuma partida tem vértice de corte. Até 30% removidos, a queda é gradual (54,5% no núcleo no pior ataque, contra 68,6% ao acaso), mas os isolados chegam a 18%. |
| Descoberta mais interessante | A estrutura global é inteiramente explicada pela divisão em partidas, para qualquer limiar testado (σ contra o aleatório por partidas entre 0,98 e 1,39); a informação está na direção das arestas (só 46,7% dos arcos são recíprocos e 13,7% dos participantes não recebem atenção acima do uniforme). |
| Extra: modelo clássico | Nenhum globalmente (BA cria hubs, ER e WS são conexos com caminhos longos); cada partida se aproxima de um Erdős–Rényi denso (p ≈ 0,67). |

![Robustez](figures/robustez.png)

![Curva de remoção](figures/robustez_curva.png)

## Como executar

Requisitos: Python 3.10+ (testado com Python 3.11.2 e as versões fixadas em `requirements.txt`).

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
./scripts/baixar_dados.sh            # baixa ~44 MB do SNAP e confere as contagens publicadas
.venv/bin/python -m src.tratamento   # regenera data/processed/
.venv/bin/python -m pytest           # testes
.venv/bin/python -m src.executar     # análises completas (~1 min); --rapido para teste
.venv/bin/python -m src.executar --so-figuras   # só redesenha as figuras a partir de results/
.venv/bin/python -m src.parte2       # ou apenas uma parte (parte1, parte2, parte3)
.venv/bin/python -m src.tabelas      # só regenera as tabelas do relatório a partir dos JSONs
.venv/bin/jupyter notebook notebooks/analise.ipynb   # consulta dos resultados
```

O relatório inclui as tabelas de `results/` com `\input` e os valores citados no texto
por macros (`results/valores.tex`). Assim, após uma nova execução, basta recompilar o PDF.

Os dados tratados já estão versionados; o download só é necessário para refazer o tratamento.

Para compilar o relatório:

```bash
cd relatorio && pdflatex main.tex && pdflatex main.tex
```

## Reprodutibilidade

- Sementes fixas (42) em todo sorteio: grafos aleatórios, remoções e modelos clássicos.
  Reexecutar reproduz todos os resultados, exceto os tempos da Parte II.
- Versões das dependências fixadas em `requirements.txt`.
- `scripts/verificar_dados.py` confere o download com o SNAP (62 partidas, 451 participantes,
  3.126.993 arestas temporais); `tests/` confere o tratamento e cada algoritmo próprio com o
  NetworkX (`python -m pytest`, 70 testes, um arquivo por módulo).
- Figuras gravadas sem data nos metadados e com ordem de vértices fixa: redesenhar não
  altera os arquivos versionados.

## Onde mexer

| Quero alterar... | Arquivo |
|---|---|
| leitura dos dados brutos, limiar de atenção, pesos | `src/tratamento.py` |
| estrutura de grafo (lista de adjacência) | `src/grafo.py` |
| um algoritmo da Parte II | `src/algoritmos.py` (testes em `tests/test_<família>.py`) |
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
