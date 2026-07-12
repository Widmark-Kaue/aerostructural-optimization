# Relatório de Análise: Otimização de Peso (`Weight`) e Desempenho

Este relatório analisa os problemas de convergência encontrados ao definir o peso total (`Weight`) como função objetivo no arquivo [optimize.py](file:///home/widmark/github/aerostructural-optimization/exam/optimize.py), explica os motivos físicos/numéricos por trás desse comportamento e propõe melhorias para reduzir o tempo de execução das otimizações.

---

## 1. Resumo dos Diagnósticos e Soluções

1. **Typo no Nome do Parâmetro (`weithg`):** O parâmetro correto esperado pelo wrapper [asalib.py](file:///home/widmark/github/aerostructural-optimization/exam/asalib.py) é `'Weight'` (ou `'weight'`). Usar `'weithg'` gera um erro `KeyError` imediato no dicionário de resultados.
2. **Causa da Não-Convergência:** A tolerância do otimizador (`tol = 1e-8`) é muito mais rigorosa do que a precisão com que o solucionador não-linear de física resolve as equações aeroestruturais (`xtol = 1e-6`). 
3. **Comportamento de Minimização do Peso:** Para minimizar o peso total ($Weight = W_{\text{seco}} + FB$), o otimizador força a espessura da asa ($t$) ao limite mínimo permitido ($1\text{ mm}$) em todos os painéis. Isso ativa 20 restrições de limite inferior simultaneamente. Nesse regime altamente restrito, a busca linear do SLSQP sofre com o ruído numérico de $10^{-6}$ e falha com o erro `Positive directional derivative for linesearch`.
4. **Solução para Convergir:** Configurar a tolerância do otimizador para `tol = 1e-5` resolve o problema. O otimizador agora **converge com sucesso** (`success: True`) em 81 iterações para o mesmo peso ideal de **$9099.32\text{ N}$**.

---

## 2. Análise Detalhada da Não-Convergência

Fizemos testes executando a otimização de peso a partir da mesma condição inicial do seu script. 

Com a configuração padrão (`tol = 1e-8`), o otimizador avança até atingir o peso mínimo possível, mas falha nas iterações finais:
* **Iteração 78:** Peso = $9099.3163\text{ N}$ | `liftExcess` = $1.59 \times 10^{-9}$ (satisfeita) | `KSmargin` = $9.04 \times 10^{-9}$ (satisfeita)
* **Iterações 79-114:** O otimizador tenta dar passos extremamente curtos para satisfazer a tolerância de $10^{-8}$. Como a física interna só converge até $10^{-6}$, o otimizador encontra apenas ruído numérico e encerra sem sucesso (`success: False`).

Ao alterar a tolerância para **`tol = 1e-5`**, o critério de parada do SLSQP se torna compatível com a precisão física e a otimização converge perfeitamente:

```
Weight Optimization finished (tol=1e-5).
     message: Optimization terminated successfully
     success: True
      status: 0
         fun: 9099.316211
         nit: 81
```

---

## 3. Explicação Física: Por que a asa vai para a espessura mínima ($1\text{ mm}$)?

Diferente da otimização de combustível (`FB`), onde a espessura da asa varia ao longo da envergadura para dar rigidez e reduzir o arrasto induzido, a otimização de peso total ($Weight = W_{\text{seco}} + FB$) se comporta de forma diferente:
* O peso estrutural da asa ($W_{\text{seco}}$) é uma função diretamente linear da espessura $t$.
* Reduzir a espessura de $5\text{ mm}$ para $1\text{ mm}$ reduz drasticamente o peso estrutural da aeronave.
* Embora a asa mais fina sofra maior deformação flexional (o que aumenta o arrasto e, consequentemente, o consumo de combustível $FB$), **a redução no peso seco é muito maior do que o aumento no consumo de combustível**.
* Portanto, o peso total mínimo é obtido quando a asa é o mais fina e leve possível ($t = 1\text{ mm}$ em todos os painéis).

---

## 4. Como Acelerar as Otimizações (Melhorias de Desempenho)

Como as otimizações dependem da solução repetida de sistemas não-lineares, elas demoram um pouco. Você pode acelerá-las com as duas estratégias abaixo:

### A. Desativar verificação de limites (Bounds-Check) no compilador Fortran
No seu script de compilação [build_f2py_diff.sh](file:///home/widmark/github/aerostructural-optimization/exam/build_f2py_diff.sh), a flag de compilação `-fbounds-check` está activa na linha 9:
```bash
FFLAGS="-O2 -fPIC -fdefault-real-8 -fbounds-check -fallow-argument-mismatch"
```
A verificação de limites em tempo de execução adiciona um overhead considerável em loops internos no Fortran. 
* **Recomendação:** Remover `-fbounds-check` e usar `-O3` para compilar o código em produção. Isso reduzirá o tempo gasto em cada avaliação de resíduo.

### B. Implementar Warm Start no Solucionador de Física (Python)
No arquivo [asalib.py](file:///home/widmark/github/aerostructural-optimization/exam/asalib.py) (linha 306), a estimativa inicial (`gama0` e `d0`) para resolver as equações físicas do ponto atual é sempre reiniciada com vetores constantes:
```python
gama0 = np.ones(self.npanels)*8.0
d0 =  np.ones(2*(self.npanels+1))*0.1
```
Como o otimizador faz alterações muito pequenas nas variáveis de projeto a cada iteração, as distribuições de circulação ($\Gamma$) e deslocamentos ($d$) da iteração anterior são uma estimativa inicial **muito melhor** do que o vetor padrão.
* **Recomendação:** Guardar o último resultado convergente de `sol.x` em um atributo da classe (ex: `self.last_state`) e usá-lo como chute inicial em `solve_asa`. Isso fará com que o método `scipy.optimize.root` convirja em apenas 2 ou 3 iterações internas em vez de 10+.
