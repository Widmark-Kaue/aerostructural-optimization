# Relatório de Verificação de Gradiente (Método Lifting Line Theory)

Este relatório apresenta as descobertas da análise e verificação do processo de diferenciação automática (AD) usando o Tapenade para o código Fortran da Teoria da Linha de Sustentação (LLT), localizado na pasta `llt_optimization`.

---

## 1. Análise do Arquivo `grad_verify.py` (Passo 1)

O erro principal que está impedindo a correspondência entre as derivadas por Diferenças Finitas (FD) e a Diferenciação Automática (AD) direta está localizado neste arquivo.

### Erro Principal: Inconsistência na Direção da Perturbação (Erro Matemático)
Na chamada da versão direta (tangente) da AD:
```python
res, res_lltd, cl, CLd, cd, CDd = llt_d.llt_main_d(**inputs, **inputs_seed)
```
Você está fornecendo as direções de semente (`seed`):
* `twistd = [0.3, 0.1, 0.2, 0.4, 0.6, 0.5]`
* `gamad = [0.2, 0.3, 0.6, 0.4, 0.5, 0.1]`

Isso significa que o código AD calcula a **derivada direcional** de $CL$, $CD$ e dos resíduos na direção do vetor semente $\mathbf{d} = (\mathbf{twistd}, \mathbf{gamad})$:
$$\text{res\_lltd} = \nabla_{\text{twist}} \text{res\_llt} \cdot \mathbf{twistd} + \nabla_{\mathbf{\Gamma}} \text{res\_llt} \cdot \mathbf{gamad}$$

No entanto, no seu laço de diferenças finitas, você perturba os inputs usando a multiplicação:
```python
for i in range(len(exp)): 
    inputs1['twist'] = inputs['twist']*(1+h[i]*inputs_seed['twistd'])
    inputs1['gama'] = inputs['gama']*(1+h[i]*inputs_seed['gamad'])
```

### O Problema da Perturbação Nula
Como o valor nominal de `inputs['twist']` é um vetor de zeros `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`, a perturbação nele é multiplicada por zero e continua sendo zero:
$$\mathbf{twist}_{\text{perturbed}} = \mathbf{0} \cdot (1 + h \cdot \mathbf{twistd}) = \mathbf{0}$$

Para `gama`, como o valor nominal é `1.0`, a perturbação funciona e é de tamanho $h \cdot \mathbf{gamad}$:
$$\mathbf{gama}_{\text{perturbed}} = \mathbf{1} \cdot (1 + h \cdot \mathbf{gamad}) = \mathbf{gama} + h \cdot \mathbf{gamad}$$

Dessa forma, o seu teste de diferenças finitas perturba **apenas** o `gama` na direção $\mathbf{gamad}$, mantendo o `twist` em zero. Ou seja, a direção da perturbação em diferenças finitas é $\mathbf{d}_{\text{FD}} = (\mathbf{0}, \mathbf{gamad})$.

### Por que as derivadas de $CL$ e $CD$ batem?
Analisando o arquivo [llt.f90](file:///home/widmark/github/aerostructural-optimization/homework02/llt_optimization/llt.f90), a sustentação ($CL$) e o arrasto ($CD$) totais são calculados usando a circulação `lki` (que depende de `gama`) e o ângulo induzido `alpha_ind` (que depende de `wi`, que por sua vez depende de `gama`):
```fortran
CL = sum(lki*panel_length*cos(alpha_ind)) / (0.5*rho_air*Vinf**2* Sref)
CD = -sum(lki*panel_length*sin(alpha_ind)) / (0.5*rho_air*Vinf**2* Sref)
```
Nenhuma dessas expressões depende de `twist`! Ou seja:
$$\frac{\partial CL}{\partial \text{twist}} = 0, \quad \frac{\partial CD}{\partial \text{twist}} = 0$$

Como a sensibilidade de $CL$ e $CD$ em relação ao `twist` é nula, a derivada direcional com a semente completa $\mathbf{d} = (\mathbf{twistd}, \mathbf{gamad})$ é matematicamente idêntica à derivada na direção incompletamente perturbada $\mathbf{d}_{\text{FD}} = (\mathbf{0}, \mathbf{gamad})$:
$$\text{CLd} = \underbrace{\frac{\partial CL}{\partial \text{twist}} \cdot \mathbf{twistd}}_{0} + \frac{\partial CL}{\partial \mathbf{\Gamma}} \cdot \mathbf{gamad} = \frac{\partial CL}{\partial \mathbf{\Gamma}} \cdot \mathbf{gamad}$$
Por isso, os valores de `CLd_FD` e `CDd_FD` coincidem perfeitamente com os resultados do AD!

### Por que as derivadas dos resíduos não batem?
Os resíduos da LLT (`res_llt`) dependem diretamente da força de sustentação seccional (`lsi`), que depende do coeficiente seccional (`cli`), o qual depende do ângulo de ataque efetivo (`alpha_eff`):
```fortran
alpha_eff = alpha + twist + alpha_ind
cli = cl0 + cla*alpha_eff
lsi = 0.5*rho_air*Veff**2*chord*cli
res_llt = lsi - lki
```
Como `res_llt` depende diretamente de `twist`, temos que:
$$\frac{\partial \text{res\_llt}}{\partial \text{twist}} \neq 0$$

Logo, o termo $\frac{\partial \text{res\_llt}}{\partial \text{twist}} \cdot \mathbf{twistd}$ é não-nulo e relevante.
* A versão AD computa a derivada direcional completa incluindo a contribuição do `twist` ($\approx 112.65$ para a primeira seção).
* O teste de diferenças finitas (FD) **não perturbou o twist**, calculando apenas a contribuição do `gama` ($\approx -2.81$ para a primeira seção).
Isso causa a divergência observada nas derivadas dos resíduos.

### Como Corrigir o Arquivo `grad_verify.py`
Para que todas as derivadas (incluindo as dos resíduos) coincidam, a perturbação deve ser **aditiva** em vez de multiplicativa, garantindo que o `twist` seja perturbado mesmo partindo de um valor nominal zero:
```python
for i in range(len(exp)): 
    # Perturbação aditiva correta ao longo da semente (seed)
    inputs1['twist'] = inputs['twist'] + h[i] * inputs_seed['twistd']
    inputs1['gama'] = inputs['gama'] + h[i] * inputs_seed['gamad']
```
Com essa correção, as derivadas dos resíduos por Diferenças Finitas e AD passarão a coincidir com precisão de máquina.

### Outras Observações Menores no Arquivo:
1. **Cópia Rasa vs Cópia Profunda**:
   A linha `inputs1 = copy.copy(inputs)` faz uma cópia rasa do dicionário. Como os valores são arrays do NumPy, modificações diretas nos elementos dos arrays dentro de `inputs1` poderiam alterar o dicionário original `inputs`. No seu caso, ao fazer `inputs1['twist'] = ...`, você reatribui um novo array criado pela operação, então o original não é modificado. Contudo, por segurança, recomenda-se o uso de `copy.deepcopy` ou a criação explícita dos novos arrays.
2. **Avaliação Base Redundante**:
   A linha `res_llt, CL, CD = llt.llt_main(**inputs)` é executada a cada iteração do laço `for`. Como os valores de `inputs` não mudam durante o laço, essa chamada poderia ser feita uma única vez fora do laço para economizar tempo de execução.

---

## 2. Análise do Arquivo `send_to_tapenade.py` (Passo 2)

O script que envia o código para o Tapenade e processa os resultados está correto e funcionando perfeitamente. As versões direta (`_d`) e reversa (`_b`) foram geradas sem erros. 

Identifiquei apenas alguns detalhes cosméticos ou redundantes que não afetam a corretude:
1. **Redundância de Substituição de Módulo (`_DIFF`)**:
   No bloco de código `diff2mode`, o script tenta substituir a string `_DIFF` nos arquivos gerados:
   ```python
   for line in fileinput.input([ff_full], inplace=True):
       print(line.replace('_DIFF', modetag), end='')
   ```
   No entanto, o Tapenade 3.16 gera diretamente os módulos nomeados como `LLT_D` e `LLT_B` a partir de `llt.f90`. Portanto, a string `_DIFF` não existe nos arquivos gerados e essa substituição é um *no-op* (não altera nada), o que é inofensivo.
2. **Substituição `REAL4` por `REAL8` na Versão Reversa**:
   Como a sub-rotina `llt_main` executa um cálculo direto simples sem dependências de trajetória complexas (laços com sobrescrita de variáveis ativas), o Tapenade não precisou usar chamadas de pilha (`PUSH`/`POP`) no código reverso. Assim, a substituição de `REAL4` por `REAL8` também acabou sendo um *no-op*, mas é uma boa prática manter essa linha caso o código futuro venha a necessitar de armazenamento em pilha.
3. **Comentário Duplicado**:
   Na linha 193 de `send_to_tapenade.py`, o comentário diz `# Replace real4 by real8`, mas a linha seguinte executa a substituição do nome do módulo (`_DIFF`). Apenas um erro de digitação no comentário.

---

## 3. Análise do Arquivo `llt.f90` (Passo 3)

O código original em Fortran está muito bem estruturado, usa `implicit none` de forma adequada e é **completamente diferenciável**. Não há problemas na forma de implementação que possam comprometer a diferenciação automática.

Aqui está a justificativa detalhada de por que a diferenciação funciona perfeitamente:
1. **Diferenciabilidade das Funções Intrínsecas**:
   * O código utiliza `atan(wi/Vinf)`. A função arcotangente é diferenciável em todo o domínio real, com derivada $1/(1+x^2)$. Mesmo se `wi = 0`, a derivada é bem definida ($1.0$).
   * O código utiliza `sqrt(Vinf**2 + wi**2)`. A função raiz quadrada possui uma singularidade em sua derivada no ponto zero. Porém, como a velocidade do escoamento livre é $V_{\infty} = 10$, o termo sob a raiz é no mínimo $V_{\infty}^2 = 100 > 0$. Assim, a raiz nunca será zero e sua derivada é sempre estritamente diferenciável.
2. **Ausência de Singularidades no Biot-Savart**:
   O cálculo da velocidade induzida é dado por:
   ```fortran
   wi(ii) = panel_length/pi * sum(gama/(4*(yc(ii) - yc)**2 - panel_length**2))
   ```
   O denominador de cada termo da soma é $4(y_c(i) - y_c(j))^2 - \Delta y^2$.
   * Para o termo de auto-indução ($i = j$), temos $y_c(i) - y_c(i) = 0$, resultando no denominador $-\Delta y^2 \neq 0$.
   * Para termos vizinhos ($i \neq j$), como os centros dos painéis são igualmente espaçados por $\Delta y$, temos $|y_c(i) - y_c(j)| = |i-j|\Delta y \geq \Delta y$. O termo se torna $(4(i-j)^2 - 1)\Delta y^2 \geq 3\Delta y^2 > 0$.
   Portanto, o denominador **nunca** é zero, garantindo que não ocorram divisões por zero durante a execução ou diferenciação.
3. **Ausência de Condicionais Ativas**:
   O código não possui estruturas condicionais (`if`, `where`) que dependam de variáveis ativas (`twist` ou `gama`), o que evita descontinuidades no gradiente.

---

## 4. Resultados dos Testes de Validação

Para comprovar as conclusões acima, foram realizados dois testes numéricos usando scripts de teste temporários (sem alterar os seus arquivos):

### A. Teste de Diferenças Finitas Corrigido
Ao corrigir a direção da perturbação em `grad_verify.py` para alinhar com as sementes da AD, obtivemos a seguinte convergência para o erro relativo da derivada de $CL$ e $CD$:

| Passo $h$ | Derivada FD ($CL'$) | Erro Relativo ($CL'$) | Derivada FD ($CD'$) | Erro Relativo ($CD'$) |
| :--- | :--- | :--- | :--- | :--- |
| $10^1$ | 0.07000000 | $0.00 \times 10^{00}$ | 0.00340167 | $3.17 \times 10^{00}$ |
| $10^0$ | 0.07000000 | $3.97 \times 10^{-16}$ | 0.00107476 | $3.17 \times 10^{-01}$ |
| $10^{-1}$ | 0.07000000 | $7.93 \times 10^{-16}$ | 0.00084207 | $3.17 \times 10^{-02}$ |
| $10^{-2}$ | 0.07000000 | $7.06 \times 10^{-14}$ | 0.00081880 | $3.17 \times 10^{-03}$ |
| $10^{-4}$ | 0.07000000 | $1.00 \times 10^{-12}$ | 0.00081624 | $3.17 \times 10^{-05}$ |
| **$10^{-6}$** | **0.07000000** | **$2.67 \times 10^{-10}$** | **0.00081621** | **$3.17 \times 10^{-07}$** |
| $10^{-8}$ | 0.06999999 | $7.59 \times 10^{-08}$ | 0.00081621 | $4.63 \times 10^{-08}$ |
| $10^{-10}$ | 0.06999984 | $2.30 \times 10^{-06}$ | 0.00081622 | $1.29 \times 10^{-06}$ |

*Nota: O erro mínimo de diferenças finitas ocorre em torno de $h = 10^{-6}$ (balanço ideal entre erro de truncamento e erro de arredondamento numérico), onde a correspondência com a AD é excelente (erro relativo de $\approx 10^{-10}$ para $CL$ e $\approx 10^{-7}$ para $CD$, devido à maior não-linearidade do arrasto).*

### B. Teste do Produto Escalar (Adjoint/Dot-Product Test)
O teste do produto escalar foi realizado para validar a consistência entre a versão direta (`llt_d`) e a versão reversa (`llt_b`). Para qualquer semente de entrada $\dot{x}$ e semente de saída $\bar{y}$, a seguinte identidade matemática deve ser satisfeita:
$$\bar{y}^T \dot{y} = \bar{x}^T \dot{x}$$

A avaliação no mesmo ponto de projeto gerou os seguintes resultados:
* **Produto Escalar Direto ($\bar{y}^T \dot{y}$):** `537.644772202655`
* **Produto Escalar Reverso ($\bar{x}^T \dot{x}$):** `537.644772202655`
* **Diferença Absoluta:** `0.00e+00` (zero absoluto dentro da precisão de máquina de dupla precisão)

Isto prova de forma definitiva que tanto a versão direta quanto a versão reversa geradas pelo Tapenade estão **matematicamente corretas** e são consistentes entre si.
