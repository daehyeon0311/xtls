# xtls-py

전이금속 **L-edge XAS(X선 흡수 분광)** 및 **선형이색성(LD)** 스펙트럼을
**전하이동 다중항 클러스터 모델(charge-transfer multiplet cluster model)** 로 계산하는 Python 코드입니다.

Arata Tanaka의 Fortran 프로그램 **XTLS**와 같은 물리 모델을 독립적으로 Python으로 구현했고,
`data/`에 넣어 둔 원조 XTLS 결과를 그래프에 겹쳐서 검증할 수 있습니다.

수록된 예제는 **Ba₂FeSi₂O₇의 Fe²⁺ (3d⁶), 사면체 배위**입니다.

| | 전하이동 배치(configuration) |
|---|---|
| 초기상태 | `2p⁶3d⁶`, `2p⁶3d⁷L̲`, `2p⁶3d⁸L̲²` |
| 종상태 | `2p⁵3d⁷`, `2p⁵3d⁸L̲`, `2p⁵3d⁹L̲²` |

---

## 설치

```bash
pip install -e .
```

의존성은 `numpy`, `scipy`, `matplotlib` 뿐입니다.
(`run_xas.py`는 `src/`를 직접 `sys.path`에 추가하므로 설치 없이도 실행됩니다.)

## 실행

1. `inputs/Fe_Ba2FeSi2O7.py`에서 파라미터를 수정합니다.
2. `run_xas.py`를 실행합니다.

```bash
python run_xas.py
```

Spyder에서는 `run_xas.py`를 열고 **F5**를 누르면 됩니다.

다른 물질을 계산하려면 입력 파일을 복사해서 고친 뒤,
`run_xas.py` 상단의 `INPUT_FILE` 을 그 파일로 바꿉니다.

```python
case_name = "Ni_NiO"
ion = "Ni2+"        # element, d 전자수, output 경로, L3 오프셋이 자동으로 채워짐
```

결과는 `outputs/<case_name>/` 에 저장됩니다 (git 추적 제외).
스펙트럼 `.txt`, 그림 `.png`, 사용된 파라미터 `.json`,
배치 에너지, 고유상태 분석 결과가 함께 나옵니다.

---

## 구조

```
run_xas.py           실행 스크립트 — 입력 검증, 계산, 상태 분석, 플롯
inputs/              계산 조건 입력 파일 (물질 하나당 하나)
src/xtls_py/
  engine/
    basis.py           Fock 기저
    operators.py       생성·소멸 연산자, 1체/2체 행렬
    sparse.py          scipy 희소행렬 버전 (큰 기저용)
    appendix_slater.py Slater 적분 표 (Appendix A, 3d/4d)
    shells.py          d·p 껍질 행렬, 스핀궤도결합, p→d 쌍극자 전이
    configuration.py   전자 배치 정의와 배치 에너지
    cluster.py         전하이동 클러스터 해밀토니안 (d + ligand + core hole)
    lanczos.py         Lanczos 삼중대각화 + 연분수 Green 함수
    transitions.py     정확대각화 전이 강도, stick → 곡선 변환
  geometry.py        배위 기하 → 결정장 / pd 혼성화
  spectrum.py        에너지 의존 Lorentzian·Gaussian 브로드닝, 파일 저장
data/                원조 XTLS 참조 스펙트럼 (오버레이 비교용)
examples/            Fe/Ba₂FeSi₂O₇ 계산 결과 예시
```

## 계산 흐름

1. 입력 파일을 읽고 검증 → 기저 크기 추정
2. 배위 기하로부터 결정장·혼성화 행렬 생성 (`geometry.py`)
3. 초기·종상태 전하이동 해밀토니안 구성 (`cluster.py`)
4. 바닥상태 근처 고유상태를 구하고 (`lowest_eigenpairs`),
   온도에 따른 Boltzmann 가중치 계산
5. 쌍극자 연산자를 적용한 뒤 **Lanczos + 연분수**로 스펙트럼 계산
   (`spectrum_method = "exact"` 로 정확대각화도 가능)
6. 편광별 곡선(x, y, z)을 실험 기하(grazing angle)에 맞춰 조합 → LD 생성
7. 브로드닝 후 저장·플롯

## 주요 입력 파라미터

| 파라미터 | 의미 |
|---|---|
| `ion`, `case_name` | 예: `"Fe2+"` — element와 d 전자수를 자동 결정 |
| `ten_dq` | 결정장 분리 10Dq (eV) |
| `delta` | 전하이동 에너지 Δ (eV) |
| `u_charge_transfer` | d-d 쿨롱 상호작용 U (eV) |
| `core_hole_potential` | core hole 인력 Q (eV) |
| `pd_sigma`, `pd_ratio` | pd 혼성화 세기와 σ/π 비 |
| `coordination_geometry` | `tetrahedral`, `octahedral`, `square_planar`, `custom_xyz` 등 |
| `max_ligand_holes` | 전하이동 섹터에 포함할 최대 ligand hole 수 |
| `*_scale` | Slater 적분·스핀궤도 축소 인자 (보통 0.8 전후) |
| `temperature_kelvin` | 초기상태 열적 분포 온도 |
| `lorentzian_hwhm`, `gaussian_hwhm` | 수명·장치 브로드닝 (에너지 의존 지정 가능) |
| `grazing_angle_deg` | 실험 LD 기하의 입사각 |

---

## 참고

- A. Tanaka and T. Jo, *J. Phys. Soc. Jpn.* **63**, 2788 (1994) — XTLS 클러스터 모델
- F. de Groot and A. Kotani, *Core Level Spectroscopy of Solids* (CRC Press, 2008)
