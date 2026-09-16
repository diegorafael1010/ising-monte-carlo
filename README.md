# ising-monte-carlo

Implementações em Python dos algoritmos de Monte Carlo (Metropolis, Wolff e
Swendsen-Wang) aplicados ao modelo de Ising bidimensional, desenvolvidas como
parte do Trabalho de Conclusão de Curso *"Comparação e Otimização de
Algoritmos de Monte Carlo Aplicados ao Modelo de Ising"* (Diego R. Oliveira,
Departamento de Física, UFPR).

O código é documentado com o objetivo de também servir como material
didático para outros alunos da graduação que queiram entender, na prática,
como uma simulação de Monte Carlo do modelo de Ising é construída do zero.
Cada função traz, em comentário, a equação correspondente do TCC, para que
o código possa ser lido lado a lado com o texto.

## Estrutura do repositório

```
ising-monte-carlo/
├── ising_utils.py                 # funções comuns aos três algoritmos (rede, vizinhos,
│                                   # energia, magnetização, solução exata de Onsager)
├── metropolis.ipynb                # algoritmo de Metropolis (dinâmica local)
├── wolff.ipynb                      # algoritmo de cluster único de Wolff
├── swendsen_wang.ipynb              # algoritmo de múltiplos clusters
├── consistencia_algoritmos.ipynb    # comparação cruzada entre os três algoritmos
├── autocorrelacao.ipynb             # tempo de autocorrelação e expoente dinâmico
├── desempenho.ipynb                 # Fator de Mérito e mapeamento de eficiência
├── vieses_metropolis.ipynb          # testes de sanidade e vieses propositais do Metropolis
└── controle_simulacoes_template.xlsx  # planilha modelo de controle das simulações
```

## Como usar no Google Colab

No início de cada notebook, clone este repositório e adicione o ao caminho
de importação do Python.

```python
!git clone https://github.com/diegorafael1010/ising-monte-carlo.git
import sys
sys.path.append('/content/ising-monte-carlo')

from ising_utils import (
    inicializar_rede,
    construir_tabela_vizinhos,
    energia_total,
    magnetizacao_total,
    energia_onsager,
    magnetizacao_onsager,
    temperatura_critica_onsager,
)
```

Se você editar `ising_utils.py` depois de já ter clonado o repositório numa
sessão do Colab, rode `!git -C /content/ising-monte-carlo pull` para
atualizar a cópia local, ou reinicie o runtime e clone novamente.

## Como usar localmente

```bash
git clone https://github.com/diegorafael1010/ising-monte-carlo.git
cd ising-monte-carlo
pip install numpy scipy matplotlib tqdm
python3 ising_utils.py   # roda a bateria de testes de sanidade do módulo
jupyter notebook          # abre os notebooks
```

## Planilha de controle

`controle_simulacoes_template.xlsx` lista todas as combinações de algoritmo,
tamanho de rede e temperatura usadas nas simulações oficiais do TCC (Seção
4.5 da Metodologia), com uma coluna de status para acompanhamento. É um
modelo em branco, faça uma cópia e preencha conforme for rodando suas
próprias simulações.

## Referência

Se este código for utilizado ou adaptado, por favor cite o TCC de origem.

> OLIVEIRA, D. R. *Comparação e Otimização de Algoritmos de Monte Carlo
> Aplicados ao Modelo de Ising*. Trabalho de Conclusão de Curso (Bacharelado
> em Física) -- Universidade Federal do Paraná, Curitiba, 2026.

## Licença

Este projeto está licenciado sob a licença MIT, que permite uso, cópia,
modificação e redistribuição livres, inclusive para fins comerciais, desde
que os créditos originais sejam mantidos. Veja o arquivo `LICENSE` para o
texto completo.
