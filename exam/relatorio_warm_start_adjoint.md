# Relatório de Análise: Aplicação de Warm Start no Solucionador Adjunto (`psi`)

Este relatório analisa a viabilidade e os impactos de aplicar a técnica de *Warm Start* no solucionador das equações adjuntas (responsáveis pelo cálculo dos gradientes) no arquivo [asalib.py](file:///home/widmark/github/aerostructural-optimization/exam/asalib.py).

---

## 1. Resumo dos Resultados dos Testes

Realizamos testes numéricos comparando o desempenho da otimização aeroestrutural com e sem a aplicação de *Warm Start* nas variáveis adjuntas ($\psi$ ou `resb`):

| Configuração | Warm Start Física? | Warm Start Adjunto? | Iterações SLSQP (`nit`) | Tempo total de execução |
| :--- | :--- | :--- | :--- | :--- |
| **Referência** | Sim | Não | **80** | **22.8 segundos** |
| **Teste Adjunto Warm Start** | Sim | Sim | **145** | **40.9 segundos** |

*Nota: Ambos os testes foram rodados sob as mesmas condições ideais de resíduos escalonados e tolerância de otimização de $10^{-8}$.*

Como demonstrado na tabela, **aplicar Warm Start no solucionador adjunto piorou o desempenho da otimização**, quase dobrando o número de iterações necessárias do SLSQP (de 80 para 145) e aumentando o tempo total de execução.

---

## 2. Por que o Warm Start do Adjunto prejudica a Otimização?

Esse fenômeno é muito conhecido em otimização numérica multidisciplinar (MDO). O motivo pelo qual o *Warm Start* funciona de forma excelente para o solver de física, mas falha para o solver adjunto, reside na **sensibilidade dos algoritmos de quase-Newton (como o BFGS usado pelo SLSQP)**:

1. **Introdução de Jitter (Ruído de Tolerância):**
   * Quando o solver `root` é iniciado com um chute inicial idêntico ao do passo anterior (muito próximo da raiz), ele atinge o critério de convergência rapidamente e para.
   * Porém, a posição exata onde ele para dentro da "bolha de tolerância" depende fortemente do histórico do chute inicial.
   * Isso introduz um pequeno erro de convergência de alta frequência (chamado de *jitter* ou ruído de histórico) no vetor adjoint $\psi$ e, consequentemente, nos gradientes calculados.

2. **Corrupção da Matriz Hessiana (BFGS):**
   * O otimizador SLSQP constrói uma aproximação da matriz Hessiana (curvatura do espaço de projeto) usando a fórmula BFGS. 
   * O BFGS baseia-se na diferença de gradientes entre iterações sucessivas:
     $$y_k = \nabla L(x_{k+1}) - \nabla L(x_k)$$
   * Como a diferença $x_{k+1} - x_k$ is muito pequena, a variação real do gradiente também é minúscula. 
   * Se o gradiente contiver mesmo que uma quantidade ínfima de *jitter* (ruído) decorrente do *warm start* do adjunto, a diferença $y_k$ será dominada pelo ruído, **corrompendo a atualização da matriz Hessiana**.
   * Sem uma estimativa de curvatura limpa, o SLSQP perde sua taxa de convergência quadrática/superlinear e passa a dar passos piores, exigindo muito mais iterações para convergir.

---

## 3. Conclusão e Recomendação

* **Física:** **Mantenha o Warm Start** (`self._last_stateVars = sol.x`). A física é resolvida de forma muito estável e não afeta a Hessiana de forma tão sensível.
* **Adjunto:** **NÃO use Warm Start** (sempre inicie o solver adjunto a partir do vetor fixo `resb0`). Iniciar a partir de um vetor fixo e limpo garante que o solucionador de sensibilidade siga sempre a mesma trajetória de convergência, gerando gradientes livres de ruído histórico e permitindo que o otimizador SLSQP convirja no menor número de passos possível.
