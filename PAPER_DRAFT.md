# 후속 논문 초안

작업 기록은 `RESEARCH_LOG.md`. 이 문서는 논문으로 쓸 형태만 정리한다.

---

## 제목 후보

1. **Strain route to the quantum-paramagnetic transition in Ba₂FeSi₂O₇ reversed by the exchange coupling**
2. Exchange, not anisotropy, controls the quantum critical point in Ba₂FeSi₂O₇
3. First-principles reassessment of the D/J phase boundary in the S = 2 melilite Ba₂FeSi₂O₇

1번이 결론을 가장 직접적으로 담는다.

---

## Abstract (초안)

> Ba₂FeSi₂O₇ is a quasi-two-dimensional S = 2 antiferromagnet that sits close to
> the boundary between an antiferromagnetic and a quantum-paramagnetic ground
> state, the two being separated by a critical ratio of the single-ion anisotropy
> D to the exchange J. A recent X-ray linear dichroism study determined D from
> full-multiplet cluster calculations and, holding J at its measured value,
> predicted that compressing the FeO₄ tetrahedron further would drive the system
> across that boundary. We compute both quantities from first principles for the
> same deformation: D from a charge-transfer cluster model and J from
> broken-symmetry density-functional theory. The two reproduce the measured
> parameters, giving D/J = 15.9 against 16.0 from inelastic neutron scattering,
> and place the material at 84% of the critical ratio. Under uniaxial strain
> along c, however, J varies by 68% across a four-degree range of the distortion
> angle while D varies by only 2%, so D/J *decreases* with compression. The
> transition is therefore reached by tensile rather than compressive strain, at
> about +4.7% along c. We trace the discrepancy to a geometric inconsistency: the
> earlier D(Δθ) was evaluated along the structural trend between Ba₂FeSi₂O₇ and
> Sr₂FeSi₂O₇, which reflects A-site substitution and moves the Fe–O distance in
> the opposite direction to uniaxial strain. We further show that Kramers
> degeneracy restricts the transition to integer-spin members of the melilite
> family, leaving Fe²⁺ as effectively the only candidate.

---

## 핵심 주장

| # | 주장 | 근거 | 그림/표 |
|---|---|---|---|
| **1** | 변형 방향이 반대다. QPM 도달에는 c축 **인장** +4.7%가 필요하다 | J +68% vs D −2% | Fig 3 |
| **2** | 원 논문에 기하 불일치가 있다 | Poisson 기하와 c축 기하가 Fe–O를 반대로 움직임 | Fig 2 |
| **3** | D와 J를 모두 제일원리로 얻어 실측을 재현했다 | D/J = 15.9 vs 16.0 | Table 1 |
| **4** | 이 계열에서 Fe²⁺가 유일한 후보다 | Kramers 정리 + Ni²⁺ 궤도 삼중항 | Table 2 |

보조 기여: 결정장 정규화 규약 의존성, SFSO THz 재현 개선, BFSO THz 예측,
4차 이방성과 터널링 분열.

---

## 섹션 구성

### 1. Introduction

- 멜리라이트 A₂MX₂O₇, 준2차원 정방격자 반강자성
- 정수 스핀 + 큰 easy-plane 이방성 → Sz = 0 단일항 바닥상태 가능
- D/J가 임계비를 넘으면 QPM. Do et al.이 α_c ≈ 0.158 (유효 S=1)로 확정
- BFSO는 임계의 84% 지점 — 변형으로 넘길 수 있다는 제안 (Choi et al.)
- 그러나 그 예측은 J를 고정한 채 D만 계산한 것. **본 연구는 둘 다 계산한다**

### 2. Methods

**2.1 Charge-transfer cluster model**
- FeO₄ + 리간드 전하이동, Fock 기저 완전 대각화
- 파라미터: Δ = 5.5, U_dd = 4.25, U_dc = 5.75 eV, Slater 82.5%
- **결정장 정규화 규약을 명시** — t₂/e 그룹을 대칭성 라벨로 나눈다.
  준위 부호로 나누는 관행은 왜곡이 커지면 불연속을 만든다 (보조자료)
- 왜곡 스캔에서는 기준 구조 한 점에서만 보정하고 점전하 모델의 절대 변화를 따른다

**2.2 Broken-symmetry DFT**
- QE 7.5, PBE + U(Fe), PAW, 60/480 Ry, 4×4×6, 24원자 셀
- E_FM − E_AFM = 8 J S(S+1)
- U = 5 eV가 실측 J를 2.6% 이내로 재현 (§3.2)

**2.3 규약 통일**
- Do et al.은 유효 S = 1 모델. 그들 식 (2) 아래에 J̃ = 3J, D̃ = D 명시
- 모든 값을 S = 2 규약으로 환산. (D/J)_c = 3/α_c = 19.0

**2.4 변형 기하**
- c축 스케일, 분수좌표 고정 → 모든 다면체가 함께 반응
- O만 옮기면 Si–O가 1.62 → 1.76 Å로 깨져 SCF가 수렴하지 않는다 (보조자료)

### 3. Results

**3.1 검증**
- 우리 구현이 원조 XTLS 스펙트럼을 5.8% RMS로 재현
- D(BFSO) = 1.45 meV (중성자 1.42)
- D(SFSO) = 1.32 meV 예측 → THz 실측(β = 1.0, γ = 1.4 THz)과 −4.4%/−9.0%.
  원 논문의 1.22 meV는 −11.5%/−15.7%

**3.2 Exchange coupling**
- J(U): U = 5.07 eV가 중성자값 재현
- D/J = 15.9 vs 실측 16.0, 임계의 84%
- **논문이 계산하지 않은 양을 채운 것이 핵심**

**3.3 D와 J의 왜곡 의존성** ← 핵심
- c축 변형 6°~10°: J +67.8%, D −2.5%
- D/J 감소 → 임계 교차 6.95°, +4.7% 인장
- Poisson 기하에서는 D +33%로 부호와 크기가 모두 다름

**3.4 계열 스크리닝**
- Mn²⁺ D ≈ 0 (궤도 단일항, 모델 검증)
- Co²⁺ D = 3.5 meV, 왜곡 의존 +594% — 그러나 Kramers라 QPM 불가
- Ni²⁺ ³T₁ 궤도 삼중항, 스핀 해밀토니안 부적절
- → Fe²⁺가 유일

### 4. Discussion

- 왜 J가 D보다 민감한가: c축 압축이 Fe–O를 줄여 Fe–O–O–Fe 초교환을 강화.
  D는 이 특정 변형에 거의 반응하지 않는다
- 원 논문의 정성적 상쇄 논증("J decreases by the former, increases by the latter")이
  깨지는 지점
- 실험 제안: +4.7% 인장 에피택시. 기판 후보 검토 필요
- 한계: D는 클러스터, J는 DFT. U 의존성. c축 균일 스케일 근사

### 5. Conclusion

---

## 그림 계획

| Fig | 내용 | 상태 |
|---|---|---|
| 1 | 결정구조 + 준위 도식 (원 논문 Fig 1 갱신) | 미작성 |
| **2** | **기하 의존성**: Fe–O와 D를 두 변형에 대해 | `outputs/paper/fig_geometry.png` ✓ |
| **3** | **J(Δθ)와 D/J 상태도** | `outputs/paper/fig_exchange.png` ✓ |
| 4 | 계열 스크리닝 D vs d 전자수 | `outputs/screening/family_scan.png` (수정 필요) |
| S1 | 정규화 규약에 따른 D(Δθ) 불연속 | 미작성 |
| S2 | XTLS 대조 스펙트럼 | `examples/` 활용 |

---

## 표 계획

**Table 1 — 파라미터 비교**

| | D (meV) | J (meV) | D/J | 출처 |
|---|---|---|---|---|
| 중성자 (Do et al.) | 1.42 | 0.0887 | 16.0 | INS + one-loop |
| **이 작업** | **1.4445** | **0.0910** | **15.88** | 클러스터 + DFT |
| 임계 | — | — | 19.0 | Do et al. α_c = 0.158 |

**Table 2 — 계열 스크리닝** (BFSO 기하, 공통 Δ, U_dd)

| 이온 | d | S | 항 | D (meV) | QPM 가능? |
|---|---|---|---|---|---|
| Mn²⁺ | 5 | 5/2 | ⁶A₁ | 0.0018 | 불가 (Kramers) |
| **Fe²⁺** | 6 | 2 | ⁵E | 1.539 | **가능** |
| Co²⁺ | 7 | 3/2 | ⁴A₂ | 3.496 | 불가 (Kramers) |
| Ni²⁺ | 8 | 1 | ³T₁ | 정의 불가 | 궤도 삼중항 |

**Table 3 — J(Δθ)와 상태도**

| Δθ | c 스케일 | Fe–O (Å) | J (meV) | D (meV) | D/J |
|---|---|---|---|---|---|
| 6.000° | 1.0979 | 2.0292 | 0.0662 | 1.4730 | 22.25 |
| 8.225° | 1.0000 | 1.9875 | 0.0910 | 1.4445 | 15.88 |
| 10.000° | 0.9247 | 1.9575 | 0.1111 | 1.4554 | 13.10 |

---

## 투고처 후보

| 저널 | 적합성 |
|---|---|
| **Phys. Rev. B** | 가장 자연스럽다. 계산 중심, 분량 여유 |
| npj Quantum Materials | 임팩트 높지만 실험 없이는 어려울 수 있음 |
| J. Korean Phys. Soc. | 원 논문과 같은 저널, 후속으로 자연스러움 |

PRB를 1순위로 본다.

---

## 남은 작업

- [ ] Fig 1 (구조 + 준위 도식) 작성
- [ ] Fig 4 계열 스크리닝 그림 정리 (Kramers 구분 표시)
- [ ] S1 정규화 규약 그림
- [ ] +4.7% 인장을 실현할 기판 후보 조사
- [ ] Ba₂CoGe₂O₇ 실제 구조로 Co²⁺ 절대값 검증 (보조 주장 강화)
- [ ] 본문 작성

## 검토가 필요한 점

1. **α_c = 0.158의 적용 범위.** Do et al.은 J̃′ = 0.1J̃, Δ̃ = 1/3 조건에서 얻었다.
   변형 하에서 이 값들이 변하면 임계비도 변한다. 현재는 고정 가정.
2. **c축 균일 스케일**이 실제 에피택셜 변형을 얼마나 재현하는가.
   면내 격자가 기판에 고정되고 c가 반응하는 것이 실제 상황이므로 방향은 맞지만,
   내부 좌표 이완을 하지 않았다.
3. **D는 클러스터 모델, J는 DFT.** 서로 다른 방법이며 규약만 통일했다.
