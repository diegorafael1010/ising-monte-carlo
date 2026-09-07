"""
ising_utils.py
===============

Funções compartilhadas para simulações de Monte Carlo do modelo de Ising
bidimensional, usadas pelos três algoritmos comparados no TCC:
*Comparação e Otimização de Algoritmos de Monte Carlo Aplicados ao Modelo de
Ising* (Diego R. Oliveira, UFPR).

Este módulo reúne apenas o que é comum aos três métodos (Metropolis, Wolff e
Swendsen-Wang): a representação do sistema físico (rede, vizinhos, energia,
magnetização) e a solução analítica exata de Onsager, usada para validação.
A lógica específica de cada algoritmo de amostragem (critério de aceitação,
construção de aglomerados etc.) permanece nos respectivos notebooks
(`metropolis.ipynb`, `wolff.ipynb`, `swendsen_wang.ipynb`).

Referências às equações citadas nos comentários correspondem à numeração do
TCC (Capítulos 2 e 4).

Uso típico, a partir de um notebook Colab (após clonar o repositório):

    import sys
    sys.path.append('/content/ising-monte-carlo')
    from ising_utils import (
        inicializar_rede,
        construir_tabela_vizinhos,
        energia_total,
        magnetizacao_total,
        probabilidade_ligacao,
        temperatura_critica_onsager,
        energia_onsager,
        magnetizacao_onsager,
    )

Este arquivo também pode ser executado diretamente (`python ising_utils.py`)
para rodar sua própria bateria de testes de sanidade.
"""

import numpy as np
from scipy.special import ellipk

# Constante de acoplamento, em unidades reduzidas (J = 1, k_B = 1 -- Seção 4.2 do TCC)
J = 1.0


# ---------------------------------------------------------------------------
# Representação da rede e estado inicial (Seção 2.3 e 4.3.1 do TCC)
# ---------------------------------------------------------------------------
def inicializar_rede(
    L: int, modo: str = "quente", rng: np.random.Generator = None
) -> np.ndarray:
    """Cria a configuração inicial de spins da rede LxL.

    Parameters
    ----------
    L : int
        Tamanho linear da rede (a rede terá N = L**2 spins).
    modo : {"fria", "quente"}
        "fria"   -> todos os spins iguais a +1 (estado ordenado, T -> 0).
        "quente" -> spins sorteados aleatoriamente em {-1, +1} (T -> infinito).
    rng : np.random.Generator, opcional
        Gerador de números aleatórios. Se None, cria um novo gerador sem
        semente fixa (útil para rodadas de produção independentes).

    Returns
    -------
    np.ndarray
        Array 1D de inteiros (dtype=int8) de tamanho N = L**2, com os spins
        armazenados na ordem k = i*L + j (Equação 2.4 do TCC, adaptada para
        indexação 0-based, padrão em Python).
    """
    if rng is None:
        rng = np.random.default_rng()

    N = L * L
    if modo == "fria":
        estado = np.ones(N, dtype=np.int8)
    elif modo == "quente":
        estado = rng.choice(np.array([-1, 1], dtype=np.int8), size=N)
    else:
        raise ValueError("modo deve ser 'fria' ou 'quente'")
    return estado


# ---------------------------------------------------------------------------
# Tabela de vizinhos com condições de contorno periódicas (Seção 2.4 do TCC)
# ---------------------------------------------------------------------------
def construir_tabela_vizinhos(L: int) -> np.ndarray:
    """Pré-computa os índices dos 4 primeiros vizinhos de cada sítio da rede,
    sob condições de contorno periódicas (PBC).

    Parameters
    ----------
    L : int
        Tamanho linear da rede.

    Returns
    -------
    np.ndarray
        Array de forma (N, 4) e dtype int32. Para o sítio k, a linha
        vizinhos[k] = [norte, sul, leste, oeste] contém os índices lineares
        dos quatro vizinhos, já considerando o "wrap-around" periódico.
    """
    N = L * L
    vizinhos = np.zeros((N, 4), dtype=np.int32)
    for i in range(L):
        for j in range(L):
            k = i * L + j
            norte = ((i - 1) % L) * L + j
            sul = ((i + 1) % L) * L + j
            leste = i * L + ((j + 1) % L)
            oeste = i * L + ((j - 1) % L)
            vizinhos[k] = [norte, sul, leste, oeste]
    return vizinhos


# ---------------------------------------------------------------------------
# Energia e magnetização globais (Equações 2.6 e 2.9 do TCC)
# ---------------------------------------------------------------------------
def energia_total(estado: np.ndarray, vizinhos: np.ndarray) -> float:
    """Calcula a energia total da configuração (Equação 2.6 do TCC).

    Custo O(N); usada para inicializar a simulação, nos instantes de
    amostragem (bem mais raros que as tentativas de atualização) e em
    testes de consistência -- nunca dentro do laço de atualização local.
    """
    soma_vizinhos = estado[vizinhos].sum(axis=1)
    energia = -J * np.sum(estado * soma_vizinhos) / 2.0
    return float(energia)


def magnetizacao_total(estado: np.ndarray) -> float:
    """Calcula a magnetização total da configuração (Equação 2.9 do TCC)."""
    return float(estado.sum())


# ---------------------------------------------------------------------------
# Probabilidade de ligação para algoritmos de cluster (Equações 3.18 / 4.17)
# ---------------------------------------------------------------------------
def probabilidade_ligacao(beta: float) -> float:
    """Probabilidade de ativar uma ligação entre dois spins paralelos
    vizinhos, usada pelos algoritmos de Wolff e Swendsen-Wang:

        P_add = 1 - exp(-2 * beta * J)   (Equações 3.18 e 4.17 do TCC)

    Não é usada pelo Metropolis (que tem sua própria lookup table de
    aceitação para os cinco valores discretos de delta_E -- ver
    metropolis.ipynb).
    """
    return 1.0 - np.exp(-2.0 * beta * J)


# ---------------------------------------------------------------------------
# Solução analítica exata de Onsager, para validação (Seção 4.4.1 do TCC)
# ---------------------------------------------------------------------------
def temperatura_critica_onsager() -> float:
    """Temperatura crítica exata de Onsager (Equação 2.15 do TCC)."""
    return 2.0 * J / np.log(1.0 + np.sqrt(2.0))


def magnetizacao_onsager(T: float) -> float:
    """Magnetização espontânea exata de Onsager (Equação 4.21 do TCC),
    válida no limite termodinâmico (L -> infinito).
    """
    Tc = temperatura_critica_onsager()
    if T >= Tc:
        return 0.0
    beta = 1.0 / T
    return (1.0 - np.sinh(2 * beta * J) ** (-4)) ** (1.0 / 8.0)


def energia_onsager(T: float) -> float:
    """Energia por sítio exata de Onsager (Equação 4.22 do TCC), válida no
    limite termodinâmico (L -> infinito).

    Nota de implementação: `scipy.special.ellipk(m)` espera o parâmetro
    m = k1**2, não o módulo k1 diretamente (Equação 4.23 do TCC usa o
    módulo k1) -- por isso chamamos `ellipk(k1**2)` abaixo.
    """
    beta = 1.0 / T
    k1 = 2 * np.sinh(2 * beta * J) / np.cosh(2 * beta * J) ** 2  # Eq. 4.23
    K = ellipk(k1**2)
    termo = 1 + (2.0 / np.pi) * (2 * np.tanh(2 * beta * J) ** 2 - 1) * K
    return -J / np.tanh(2 * beta * J) * termo


# ---------------------------------------------------------------------------
# Testes de sanidade do próprio módulo
# ---------------------------------------------------------------------------
def _testes_sanidade() -> None:
    """Roda uma bateria mínima de testes de consistência das funções acima.
    Executado automaticamente quando este arquivo é rodado diretamente
    (`python ising_utils.py`).
    """
    rng = np.random.default_rng(42)

    # --- inicializar_rede ---
    fria = inicializar_rede(4, modo="fria")
    assert np.all(fria == 1), "Falha: rede fria deveria ter todos os spins +1"

    quente = inicializar_rede(64, modo="quente", rng=rng)
    assert abs(quente.mean()) < 0.2, "Magnetização da rede quente muito longe de 0"

    # --- construir_tabela_vizinhos ---
    viz = construir_tabela_vizinhos(3)
    assert viz[0, 0] == 6, "Falha na periodicidade do vizinho norte (L=3)"

    # --- energia_total / magnetizacao_total ---
    L_teste = 4
    estado = inicializar_rede(L_teste, modo="fria")
    vizinhos = construir_tabela_vizinhos(L_teste)
    N_teste = L_teste**2
    assert energia_total(estado, vizinhos) == -J * 2 * N_teste
    assert magnetizacao_total(estado) == N_teste

    # --- probabilidade_ligacao ---
    assert probabilidade_ligacao(beta=0.0) == 0.0, "P_add deveria ser 0 em T infinita"
    assert probabilidade_ligacao(beta=1e6) > 0.999, "P_add deveria ser ~1 em T baixa"

    # --- Onsager: Tc conhecido, magnetização e energia com sinais/limites corretos ---
    Tc = temperatura_critica_onsager()
    assert abs(Tc - 2.269185) < 1e-5, f"Tc calculado ({Tc}) difere do esperado"
    assert magnetizacao_onsager(Tc + 0.5) == 0.0, "m deveria ser 0 acima de Tc"
    assert magnetizacao_onsager(0.5) > 0.99, "m deveria ser ~1 em T baixa"
    assert energia_onsager(0.5) < energia_onsager(5.0), "energia deveria crescer com T"

    print("OK: todos os testes de sanidade de ising_utils.py passaram.")


if __name__ == "__main__":
    _testes_sanidade()
