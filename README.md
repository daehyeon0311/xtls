# xtls-py

전이금속 화합물의 **X선 흡수 분광(XAS)·선형이색성(LD)** 과 **광전자분광(XPS)** 스펙트럼을
**전하이동 다중항 클러스터 모델(charge-transfer multiplet cluster model)** 로 계산하는 Python 코드입니다.

Arata Tanaka의 Fortran 프로그램 **XTLS**와 같은 물리 모델을 독립적으로 Python으로 구현했습니다.
저장소에 함께 둔 `Xtls900_man_Jpn.pdf`가 원조 프로그램의 매뉴얼입니다.

| 실행 스크립트 | 분광법 | XTLS 대응 |
|---|---|---|
| `run_xas.py` | L-edge XAS, 선형이색성 | `Mode=XAS` |
| `run_xps.py` | 코어준위 XPS, 원자가띠 XPS | `Mode=XPS` (`PES`) |

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

입력 파일에서 파라미터를 수정한 뒤 해당 스크립트를 실행합니다.

```bash
python run_xas.py    # inputs/Fe_Ba2FeSi2O7.py
python run_xps.py    # inputs/xps_NiO.py
```

Spyder에서는 스크립트를 열고 **F5**를 누르면 됩니다.
다른 입력 파일을 쓰려면 스크립트 상단의 `INPUT_FILE`을 바꿉니다.

원소·원자가는 `ion` 한 줄로 지정하면 나머지가 자동으로 채워집니다.

```python
ion = "Ni2+"        # element, d 전자수, output 경로, 플롯 오프셋 자동 결정
```

결과는 `outputs/<case_name>/`에 저장됩니다 (git 추적 제외).
스펙트럼 `.txt`, 그림 `.png`, 사용 파라미터 `.json`, 배치 에너지, 고유상태 분석이 함께 나옵니다.

### 수록된 입력 예제

| 파일 | 내용 |
|---|---|
| `inputs/Fe_Ba2FeSi2O7.py` | Ba₂FeSi₂O₇의 Fe²⁺ (d⁶, 사면체) L-edge XAS/LD |
| `inputs/xps_Fe_Ba2FeSi2O7.py` | 같은 클러스터의 Fe 2p XPS |
| `inputs/xps_NiO.py` | NiO의 Ni 2p XPS — **XTLS 매뉴얼 예제 재현** |

`xps_NiO.py`는 매뉴얼 20–21쪽의 X-card를 그대로 옮긴 것으로, 구현 검증용 기준입니다.

---

## 구조

```
run_xas.py           XAS/LD 실행 스크립트
run_xps.py           XPS 실행 스크립트
inputs/              계산 조건 입력 파일 (케이스 하나당 하나)
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

## 참고

- A. Tanaka, *Xtls version 9.0 매뉴얼* (2009) — 저장소의 `Xtls900_man_Jpn.pdf`
- A. Tanaka and T. Jo, *J. Phys. Soc. Jpn.* **63**, 2788 (1994) — XTLS 클러스터 모델
- F. de Groot and A. Kotani, *Core Level Spectroscopy of Solids* (CRC Press, 2008)
