# xtls-py

전이금속 화합물의 **X선 흡수 분광(XAS)·선형이색성(LD)** 과 **광전자분광(XPS)** 스펙트럼을
**전하이동 다중항 클러스터 모델(charge-transfer multiplet cluster model)** 로 계산하는 Python 코드입니다.

Arata Tanaka의 Fortran 프로그램 **XTLS**와 같은 물리 모델을 독립적으로 Python으로 구현했습니다.
저장소에 함께 둔 `Xtls900_man_Jpn.pdf`가 원조 프로그램의 매뉴얼입니다.

| 모드 | 분광법 | XTLS 대응 |
|---|---|---|
| `xas` | L-edge XAS, 선형이색성 | `Mode=XAS` |
| `xps` | 코어준위 XPS, 원자가띠 XPS | `Mode=XPS` (`PES`) |
| `both` | 둘 다 연속 실행 | — |

두 계산은 **같은 파라미터 세트**(Δ, U_dd, U_dc, 10Dq, 혼성화)를 공유합니다.
XAS와 XPS를 동시에 맞추는 것이 이 모델의 표준 사용법입니다 —
2p XPS의 satellite 구조가 XAS 선형보다 Δ에 훨씬 민감하기 때문입니다.

### 전하이동 배치

d 전자수 `n`인 이온에 대해:

| | 초기상태 | 종상태 |
|---|---|---|
| XAS | `2p⁶3d ⁿ`, `2p⁶3dⁿ⁺¹L̲`, … | `2p⁵3dⁿ⁺¹`, `2p⁵3dⁿ⁺²L̲`, … |
| 코어준위 XPS | 〃 | `2p⁵3d ⁿ`, `2p⁵3dⁿ⁺¹L̲`, … |
| 원자가띠 XPS | 〃 | `2p⁶3dⁿ⁻¹`, `2p⁶3d ⁿL̲`, … |

XAS는 코어 전자가 3d로 **올라가고**, XPS는 광전자가 클러스터를 **떠납니다**.
그래서 XPS 종상태는 전자가 하나 적고, d 전자수는 초기와 같습니다.

---

## 설치

```bash
pip install -e .
```

의존성은 `numpy`, `scipy`, `matplotlib` 뿐입니다.
(실행 스크립트가 `src/`를 직접 `sys.path`에 추가하므로 설치 없이도 실행됩니다.)

## 실행

진입점은 `run.py` 하나입니다.

```bash
python run.py                                  # 상단 설정대로
python run.py inputs/Fe_Ba2FeSi2O7.py xas
python run.py inputs/Fe_Ba2FeSi2O7.py xps
python run.py inputs/NiO.py both               # 두 분광법 연속 실행
```

Spyder에서는 `run.py`를 열어 상단 두 줄만 고치고 **F5**를 누르면 됩니다.

```python
INPUT_FILE = ROOT / "inputs" / "Fe_Ba2FeSi2O7.py"
SPECTROSCOPY = "xas"        # "xas" | "xps" | "both"
```

`run_xas.py`, `run_xps.py`를 직접 실행해도 됩니다. 같은 입력 파일을 읽습니다.

### 입력 파일 = 물질 하나

한 파일이 두 분광법을 모두 기술합니다.

```python
# 공통 — 클러스터 자체. 한 번만 씁니다.
ion = "Ni2+"                # element, d전자수, 출력경로, 플롯 오프셋 자동 결정
delta = 4.7
u_charge_transfer = 7.3
coordination_geometry = "octahedral"

# 분광법별로 값이 달라야 하는 것만 접두사를 붙입니다.
xas_energy_min = -15.0
xps_energy_min = -3.0

# 한쪽만 아는 설정은 접두사가 필요 없습니다 (상대는 무시).
grazing_angle_deg = 23.5        # XAS 전용
photoemission_shell = "2p"      # XPS 전용
```

접두사 붙은 값이 없는 값을 이깁니다. 따라서 Δ를 한 번 고치면 두 스펙트럼이 함께 움직입니다 —
XAS와 XPS를 한 파라미터 세트로 동시에 맞추는 것이 이 모델의 표준 사용법이므로, 이 구조가 불일치를 원천 차단합니다.

결과는 `outputs/<case_name>/`에 저장됩니다 (git 추적 제외).
스펙트럼 `.txt`, 그림 `.png`, 사용 파라미터 `.json`, 배치 에너지, 고유상태 분석이 함께 나옵니다.
XAS와 XPS 결과는 파일명이 달라 같은 폴더에 공존합니다.

### 수록된 입력 예제

| 파일 | 내용 |
|---|---|
| `inputs/Fe_Ba2FeSi2O7.py` | Ba₂FeSi₂O₇의 Fe²⁺ (d⁶, 사면체) — XAS/LD + 2p XPS |
| `inputs/NiO.py` | NiO의 Ni²⁺ (d⁸, 팔면체) — **XPS는 XTLS 매뉴얼 예제 재현** |

`NiO.py`의 XPS 블록은 매뉴얼 20–21쪽 X-card를 그대로 옮긴 것으로, 구현 검증용 기준입니다.
(매뉴얼에 NiO XAS 카드는 없어서, 그쪽은 미검증입니다.)

---

## 구조

```
run.py               진입점 — 입력 파일과 분광법 선택
run_xas.py           XAS/LD 계산 (단독 실행도 가능)
run_xps.py           XPS 계산 (단독 실행도 가능)
inputs/              계산 조건 입력 파일 (물질 하나당 하나)
src/xtls_py/
  engine/
    basis.py           Fock 기저 (정수 비트패턴)
    operators.py       생성·소멸 연산자, 조밀 1체/2체 행렬
    sparse.py          벡터화된 희소행렬 커널 — 해밀토니안과 전이 연산자
    appendix_slater.py Slater 적분 표 (Appendix A, 3d/4d)
    shells.py          d·p 껍질 행렬, 스핀궤도결합, p→d 쌍극자
    configuration.py   전자 배치 정의와 배치 에너지
    cluster.py         전하이동 클러스터 (d + ligand + core hole)
    lanczos.py         Lanczos 삼중대각화 + 연분수 Green 함수
  geometry.py        배위 기하 → 결정장 / pd 혼성화
  spectrum.py        에너지 의존 브로드닝, 파일 저장
data/                원조 XTLS 참조 스펙트럼 (오버레이 비교용)
examples/            Fe/Ba₂FeSi₂O₇ XAS 계산 결과 예시
Xtls900_man_Jpn.pdf  원조 XTLS 9.0 매뉴얼 (일본어)
```

## 계산 흐름

1. 입력 파일을 읽고 검증 → 기저 크기 추정
2. 배위 기하로부터 결정장·혼성화 행렬 생성 (`geometry.py`)
3. 초기·종상태 전하이동 해밀토니안 구성 (`cluster.py`)
4. 바닥상태 근처 고유상태를 구하고 온도에 따른 Boltzmann 가중치 계산
5. 전이 연산자를 적용한 뒤 **Lanczos + 연분수**로 스펙트럼 계산
   (`spectrum_method = "exact"`로 정확대각화도 가능)
   - XAS: 쌍극자 연산자, 편광 x/y/z 3채널 → 실험 기하에 맞춰 LD 조합
   - XPS: 소멸 연산자, 코어 스핀궤도 6채널(또는 원자가 10채널)을 비간섭 합산
6. 브로드닝 후 저장·플롯

## 주요 입력 파라미터

### 공통

| 파라미터 | 의미 |
|---|---|
| `ion`, `case_name` | 예: `"Fe2+"` — element와 d 전자수를 자동 결정 |
| `ten_dq` | 결정장 분리 10Dq (eV) |
| `delta` | 전하이동 에너지 Δ (eV) |
| `u_charge_transfer` | d-d 쿨롱 상호작용 U_dd (eV) |
| `core_hole_potential` | core hole 인력 U_dc (eV) |
| `coordination_geometry` | `tetrahedral`, `octahedral`, `square_planar`, `custom_xyz` 등 |
| `hybridization_mode` | `"symmetry"` (V(e_g), V(t₂g) 직접 지정), `"geometry"` (기하에서 계산), `"scalar"` |
| `max_ligand_holes` | 전하이동 섹터에 포함할 최대 ligand hole 수 |
| `*_scale` | Slater 적분·스핀궤도 축소 인자 (XTLS 기본 `Red = 0.8`) |
| `temperature_kelvin` | 초기상태 열적 분포 온도 |
| `lorentzian_hwhm`, `gaussian_hwhm` | 수명·장치 브로드닝 (에너지 의존 지정 가능) |

### XAS 전용

| 파라미터 | 의미 |
|---|---|
| `grazing_angle_deg`, `inplane_curve` | 실험 LD 기하 |
| `pd_sigma`, `pd_ratio` | pd 혼성화 세기와 σ/π 비 (`hybridization_mode="geometry"`) |

### XPS 전용

| 파라미터 | 의미 |
|---|---|
| `photoemission_shell` | `"2p"` (코어준위), `"3d"` (원자가띠), `"ligand"` |
| `spin_resolved` | 스핀 분해 스펙트럼 함께 출력 |
| `ligand_ten_dq` | 리간드 분자궤도 분열 — XTLS의 `10Dq(Ld) = 2·T_pp` |
| `v_eg`, `v_t2g` | 대칭성별 혼성화 — XTLS의 `VOh(Ld 3d) = {V(e_g), V(t₂g)}` |
| `plot_binding_energy_axis` | 결합에너지 축 반전 (XPS 관례) |

---

## 검증

- **NiO 2p XPS**: 매뉴얼 예제 파라미터로 2p₃/₂ 주선 – satellite 간격 6.3 eV,
  스핀궤도 분리 17.5 eV를 재현합니다 (ζ₂p·3/2 = 17.3 eV).
- **합 규칙**: 코어준위 XPS의 총 스펙트럼 세기가 껍질 점유수(2p⁶ → 6.000)와 일치합니다.
- **Appendix-A 표**: 매뉴얼이 인쇄한 NiO Slater 적분값과 소수점 셋째 자리까지 일치합니다.
- **교차 검증**: 같은 입력의 XAS와 XPS가 동일한 초기 바닥상태 에너지를 줍니다.

### 축퇴 상태에 대한 주의

`lowest_eigenpairs`는 ARPACK 시작 벡터를 고정해 실행 간 재현성을 보장합니다.
다만 이것은 선택을 **고정**할 뿐이며, 축퇴 부공간 안에서는 어떤 기저도 똑같이 유효합니다.

따라서 부공간 불변량이 아닌 값 — `state_analysis.txt`의 개별 상태 Sz·Lz,
스핀 분해 XPS 채널, 정방정계가 아닌 경우의 x/y 편광 — 은 물리적으로 유일하게 정해지지 않습니다.
합성량(`iso`, `ld`, `total`)은 영향받지 않습니다.
개별 상태량이 필요하시면 축퇴 부공간 안에서 Sz를 대각화해 기저를 고정하는 후처리가 필요합니다.

## 참고

- A. Tanaka, *Xtls version 9.0 매뉴얼* (2009) — 저장소의 `Xtls900_man_Jpn.pdf`
- A. Tanaka and T. Jo, *J. Phys. Soc. Jpn.* **63**, 2788 (1994) — XTLS 클러스터 모델
- F. de Groot and A. Kotani, *Core Level Spectroscopy of Solids* (CRC Press, 2008)
