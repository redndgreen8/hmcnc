# HMCNC Mathematical Models

## Hidden Markov Model Framework

### State Space
- **S** = {s₀, s₁, s₂, s₃, s₄, s₅, s₆} representing copy numbers 0-6
- **Normal state**: s₂ (diploid, CN=2)

### Observations
At each genomic bin t:
- **cₜ** = coverage (read depth)
- **Lₜ** = left-clip count
- **Rₜ** = right-clip count
- **BAFₜ** = B-allele frequency (optional)

---

## Emission Distributions

### Coverage: Negative Binomial

For copy number state sᵢ with CN=i:

```
P(c | sᵢ) = NB(c; μᵢ, φ)
```

Where:
- **μᵢ** = μ_haploid × i (mean scales with copy number)
- **φ** = dispersion parameter (estimated from data)
- **μ_haploid** = genome-wide haploid mean coverage

**Log-likelihood:**
```
log P(c | sᵢ) = log Γ(c + φ) - log Γ(φ) - log(c!)
              + φ log(φ/(μᵢ + φ)) + c log(μᵢ/(μᵢ + φ))
```

### Clipping: Zero-Inflated Negative Binomial (ZINB)

Most genomic bins have zero clips; ZINB models this excess:

```
P(c | state) = π · I(c=0) + (1-π) · NB(c; μ_clip, φ_clip)
```

Where:
- **π** = zero-inflation probability (typically 0.8-0.95)
- **μ_clip** = mean clip count when clips present
- **φ_clip** = clip count dispersion

**Log-likelihood:**
```
If c = 0:
  log P(0) = log(π + (1-π) · NB(0; μ, φ))

If c > 0:
  log P(c) = log(1-π) + log NB(c; μ, φ)
```

---

## Directional Clipping Model

### Separate Left and Right Clips

Each bin has two clip observations:
- **P_L** = P(left clips | position)
- **P_R** = P(right clips | position)

Each modeled independently with ZINB:
```
P_L ~ ZINB(π_L, μ_L, φ_L)
P_R ~ ZINB(π_R, μ_R, φ_R)
```

### Biological Interpretation

| Event | Left Clips | Right Clips |
|-------|-----------|-------------|
| Deletion start | High | Low |
| Deletion end | Low | High |
| Duplication start | Low | High |
| Duplication end | High | Low |

---

## Transition Model

### Base Transitions

Standard HMM transition matrix **A_base** with:
- High self-transition: P(sᵢ → sᵢ) ≈ 0.999
- Small adjacent transitions: P(sᵢ → sᵢ±₁) ≈ 0.0005
- Very small non-adjacent: P(sᵢ → sⱼ) ≈ 10⁻⁶ for |i-j| > 1

### Clipping-Modified Transitions

At position t, modify transitions based on local clipping:

```
aₜ,ᵢⱼ = a_base,ᵢⱼ × (1 + α × mᵢⱼ(P_L, P_R))
```

Where:
- **α** = clipping influence weight (tunable, ~0.1-0.5)
- **mᵢⱼ(P_L, P_R)** = direction-aware modifier

### Modifier Function

```
mᵢⱼ(P_L, P_R) = {
  P_L - P_R   if transition i→j expects left clips
  P_R - P_L   if transition i→j expects right clips
  0           if transition is neutral
}
```

**Expected clip patterns:**
- s₂ → s₁ (entering deletion): expects P_L
- s₁ → s₂ (exiting deletion): expects P_R
- s₂ → s₃ (entering duplication): expects P_R
- s₃ → s₂ (exiting duplication): expects P_L

---

## Baum-Welch Algorithm

### E-Step: Forward-Backward

**Forward variable:**
```
αₜ(i) = P(o₁...oₜ, qₜ=sᵢ | λ)

α₁(i) = πᵢ · bᵢ(o₁)
αₜ(i) = [Σⱼ αₜ₋₁(j) · aⱼᵢ] · bᵢ(oₜ)
```

**Backward variable:**
```
βₜ(i) = P(oₜ₊₁...oₜ | qₜ=sᵢ, λ)

βₜ(i) = 1
βₜ(i) = Σⱼ aᵢⱼ · bⱼ(oₜ₊₁) · βₜ₊₁(j)
```

**Posterior:**
```
γₜ(i) = P(qₜ=sᵢ | O, λ) = αₜ(i) · βₜ(i) / P(O|λ)
```

### Dual Forward-Backward

Run two passes:
1. **Neutral**: Using A_base
2. **Clipped**: Using Aₜ (position-dependent)

Combine posteriors:
```
γₜ(i) = wₙ(t) · γₜⁿᵉᵘᵗʳᵃˡ(i) + wc(t) · γₜᶜˡⁱᵖᵖᵉᵈ(i)
```

Where weights wₙ, wc depend on local clipping evidence.

### M-Step: Parameter Updates

**Transition probabilities:**
```
âᵢⱼ = Σₜ ξₜ(i,j) / Σₜ γₜ(i)
```

**Emission parameters (coverage):**
```
μ̂ᵢ = Σₜ γₜ(i) · cₜ / Σₜ γₜ(i)
```

---

## Newton-Raphson for Dispersion

### Problem
Estimate φ in NB or ZINB by maximizing log-likelihood.

### Iteration
```
φₖ₊₁ = φₖ - f(φₖ) / f'(φₖ)
```

### For Negative Binomial

**Score function:**
```
f(φ) = ∂ℓ/∂φ = Σᵢ [ψ(cᵢ + φ) - ψ(φ) + log(φ/(μ+φ)) + (μ-cᵢ)/(μ+φ)]
```

**Hessian:**
```
f'(φ) = ∂²ℓ/∂φ² = Σᵢ [ψ'(cᵢ + φ) - ψ'(φ) + 1/φ - 2/(μ+φ) + (μ-cᵢ)/(μ+φ)²]
```

Where:
- ψ(x) = digamma function = d/dx log Γ(x)
- ψ'(x) = trigamma function

### For ZINB

More complex due to zero-inflation; typically use EM or profile likelihood.

### Convergence
- Tolerance: |φₖ₊₁ - φₖ| < ε (e.g., ε = 10⁻⁶)
- Max iterations: 50-100
- Initialize: φ₀ from method of moments

---

## Method of Moments Initialization

### Negative Binomial
Given sample mean x̄ and variance s²:
```
μ̂ = x̄
φ̂ = x̄² / (s² - x̄)   if s² > x̄
```

### ZINB
```
π̂ = (observed zeros - expected NB zeros) / n
μ̂ = x̄ / (1 - π̂)
φ̂ = from NB method of moments on non-zero data
```

---

## Log-Space Computations

All probabilities computed in log space for numerical stability:

**Log-sum-exp:**
```
log(Σᵢ exp(aᵢ)) = max(a) + log(Σᵢ exp(aᵢ - max(a)))
```

**Pair sum:**
```cpp
double PairSumOfLogP(double a, double b) {
    if (a > b) return a + log1p(exp(b - a));
    else return b + log1p(exp(a - b));
}
```

---

## References

1. Rabiner, L.R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition.
2. Yau, C. et al. (2011). A statistical approach for detecting genomic aberrations in heterogeneous tumor samples from single nucleotide polymorphism genotyping data.
3. Mathematical specifications in `definitions/*.pdf`
