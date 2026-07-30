"""Interactive front end for the charge-transfer cluster calculations.

    streamlit run app.py

Drag a parameter, watch the spectrum move. The point is fitting: pin curves to
compare, overlay measured data to see the residual, and show XAS and XPS side
by side so one parameter set has to satisfy both. When a set looks right,
export it back to an input file so `run.py` reproduces it at full quality.

Nothing here reimplements the physics -- it drives `run_xas.calculate_spectrum`
and `run_xps.calculate_spectrum` directly, so the curves shown are the same
ones the command-line runners produce.
"""

from __future__ import annotations

import io
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _configure_font() -> None:
    """Pick a font that can draw the Korean plot labels.

    matplotlib's default has no Hangul glyphs, so legends render as tofu boxes.
    Falling through the list leaves the default in place, which still draws the
    curves correctly.
    """
    from matplotlib import font_manager

    available = {font.name for font in font_manager.fontManager.ttflist}
    for name in ("Malgun Gothic", "AppleGothic", "NanumGothic", "Noto Sans CJK KR", "Noto Sans KR"):
        if name in available:
            plt.rcParams["font.family"] = name
            break
    # Hangul fonts often lack U+2212, so keep the ASCII hyphen for minus signs.
    plt.rcParams["axes.unicode_minus"] = False


_configure_font()


# ---------------------------------------------------------------------------
# Quality presets.
#
# Exploring at h=1 on a coarse grid takes a few seconds, which is fast enough
# to drag a slider against. The final pass uses whatever the input file says.

QUALITY_PRESETS = {
    "탐색 (h=1, 거친 격자)": {
        "max_ligand_holes": 1,
        "n_recursion": 200,
        "energy_step": 0.05,
        "n_initial_states": 3,
        "n_analyzed_states": 3,
    },
    "중간 (h=1, 세밀)": {
        "max_ligand_holes": 1,
        "n_recursion": 300,
        "energy_step": 0.02,
        "n_initial_states": 10,
        "n_analyzed_states": 5,
    },
    "최종 (입력 파일 그대로)": {},
}

# Sliders. Each entry is (label, minimum, maximum, step).
CLUSTER_SLIDERS = {
    "delta": ("Δ  전하이동 에너지", -5.0, 15.0, 0.05),
    "u_charge_transfer": ("U_dd  d-d 쿨롱", 0.0, 12.0, 0.05),
    "core_hole_potential": ("U_dc  core hole 인력", 0.0, 12.0, 0.05),
    "ten_dq": ("10Dq  결정장", -2.0, 3.0, 0.01),
    "ligand_ten_dq": ("10Dq(L)  리간드 분열", -3.0, 3.0, 0.01),
}

HYBRIDIZATION_SLIDERS = {
    "geometry": {
        "pd_sigma": ("pdσ", -4.0, 0.0, 0.05),
        "pd_ratio": ("pdπ/pdσ 비", -4.0, 0.0, 0.05),
    },
    "symmetry": {
        "v_eg": ("V(e_g)", 0.0, 4.0, 0.05),
        "v_t2g": ("V(t_2g)", 0.0, 4.0, 0.05),
    },
    "scalar": {
        "hopping": ("hopping", 0.0, 4.0, 0.05),
    },
}

BROADENING_SLIDERS = {
    "lorentzian_hwhm": ("Lorentzian HWHM", 0.02, 1.5, 0.01),
    "gaussian_hwhm": ("Gaussian HWHM", 0.0, 2.0, 0.01),
}

SLATER_KEYS = ("fdd2_scale", "fdd4_scale", "fpd2_scale", "gpd1_scale", "gpd3_scale")


# ---------------------------------------------------------------------------
# Calculation.


def _module(mode: str):
    if mode == "xas":
        import run_xas as module
    else:
        import run_xps as module
    return module


@st.cache_data(show_spinner=False, max_entries=64)
def compute(mode: str, input_file: str, overrides_json: str):
    """Run one spectrum. Cached on the full parameter set."""
    module = _module(mode)
    module._load_input_file(input_file)
    module.__dict__.update(json.loads(overrides_json))
    # Broadening point lists would override the scalar widths set above.
    module.__dict__["lorentzian_hwhm_points"] = []
    module.__dict__["gaussian_hwhm_points"] = []
    module.__dict__["gaussian_sigma"] = 0.0
    energy, curves, metadata = module.calculate_spectrum()
    return energy, {name: np.asarray(values) for name, values in curves.items()}, metadata


@st.cache_data(show_spinner=False)
def defaults_for(mode: str, input_file: str) -> dict:
    """Parameter values an input file starts from, after prefix resolution."""
    module = _module(mode)
    module._load_input_file(input_file)
    keys = (
        list(CLUSTER_SLIDERS)
        + list(BROADENING_SLIDERS)
        + list(SLATER_KEYS)
        + ["hybridization_mode", "energy_min", "energy_max", "energy_shift", "so3d_scale", "so2p_scale"]
        + [key for group in HYBRIDIZATION_SLIDERS.values() for key in group]
    )
    return {key: module.__dict__.get(key) for key in keys}


# ---------------------------------------------------------------------------
# Sidebar.


def sidebar() -> dict:
    st.sidebar.title("xtls-py")

    inputs = sorted((ROOT / "inputs").glob("*.py"))
    if not inputs:
        st.sidebar.error("inputs/ 에 입력 파일이 없습니다.")
        st.stop()
    input_file = st.sidebar.selectbox("물질", inputs, format_func=lambda path: path.stem)

    view = st.sidebar.radio("분광법", ["XAS", "XPS", "XAS + XPS"], horizontal=True)
    modes = {"XAS": ["xas"], "XPS": ["xps"], "XAS + XPS": ["xas", "xps"]}[view]

    quality = st.sidebar.selectbox("계산 품질", list(QUALITY_PRESETS), index=0)
    auto = st.sidebar.checkbox("슬라이더 조작 시 자동 계산", value=True)

    base = defaults_for(modes[0], str(input_file))

    st.sidebar.divider()
    st.sidebar.subheader("클러스터")
    overrides: dict = {}
    for key, (label, low, high, step) in CLUSTER_SLIDERS.items():
        value = base.get(key)
        if value is None:
            continue
        overrides[key] = st.sidebar.slider(label, low, high, float(value), step, key=f"s_{key}")

    hybridization_mode = st.sidebar.selectbox(
        "혼성화 방식",
        list(HYBRIDIZATION_SLIDERS),
        index=list(HYBRIDIZATION_SLIDERS).index(base.get("hybridization_mode", "geometry")),
    )
    overrides["hybridization_mode"] = hybridization_mode
    for key, (label, low, high, step) in HYBRIDIZATION_SLIDERS[hybridization_mode].items():
        value = base.get(key)
        if value is None:
            continue
        overrides[key] = st.sidebar.slider(label, low, high, float(value), step, key=f"s_{key}")

    st.sidebar.divider()
    st.sidebar.subheader("다중항·브로드닝")
    slater = st.sidebar.slider(
        "Slater 축소 인자",
        0.5,
        1.0,
        float(base.get("fdd2_scale") or 0.8),
        0.005,
        help="F²dd, F⁴dd, F²pd, G¹pd, G³pd 에 일괄 적용됩니다. XTLS 기본값은 0.8.",
    )
    overrides.update({key: slater for key in SLATER_KEYS})
    for key, (label, low, high, step) in BROADENING_SLIDERS.items():
        value = base.get(key)
        if value is None:
            continue
        overrides[key] = st.sidebar.slider(label, low, high, float(value), step, key=f"s_{key}")

    st.sidebar.divider()
    with st.sidebar.expander("에너지 범위"):
        overrides["energy_min"] = st.number_input("최소 (eV)", value=float(base.get("energy_min", -20.0)))
        overrides["energy_max"] = st.number_input("최대 (eV)", value=float(base.get("energy_max", 20.0)))
        overrides["energy_shift"] = st.number_input("이동 (eV)", value=float(base.get("energy_shift", 0.0)))

    return {
        "input_file": str(input_file),
        "modes": modes,
        "quality": quality,
        "auto": auto,
        "overrides": overrides,
    }


# ---------------------------------------------------------------------------
# Experimental overlay.


def load_experiment(upload) -> tuple[np.ndarray, np.ndarray, list[str]] | None:
    """Read a two-column or headed multi-column text file."""
    if upload is None:
        return None
    raw = upload.getvalue().decode("utf-8", errors="replace") if hasattr(upload, "getvalue") else Path(upload).read_text()
    lines = [line for line in raw.splitlines() if line.strip()]
    if not lines:
        return None
    header: list[str] = []
    try:
        float(lines[0].split()[0])
    except (ValueError, IndexError):
        header = lines[0].replace(",", " ").split()
        lines = lines[1:]
    data = np.loadtxt(io.StringIO("\n".join(lines)))
    if data.ndim == 1:
        data = data[:, None]
    if not header:
        header = [f"col{i}" for i in range(data.shape[1])]
    return data[:, 0], data[:, 1:], header[1:] if len(header) > 1 else header


def residual(energy, curve, exp_energy, exp_values) -> float:
    """RMS difference after matching peak heights, in percent of the maximum."""
    overlap = (exp_energy >= energy.min()) & (exp_energy <= energy.max())
    if overlap.sum() < 5:
        return float("nan")
    interpolated = np.interp(exp_energy[overlap], energy, curve)
    reference = exp_values[overlap]
    scale = np.max(np.abs(interpolated)) or 1.0
    exp_scale = np.max(np.abs(reference)) or 1.0
    return 100.0 * float(np.sqrt(np.mean((interpolated / scale - reference / exp_scale) ** 2)))


# ---------------------------------------------------------------------------
# Input-file export.


def rewrite_input(text: str, overrides: dict, prefix: str) -> str:
    """Substitute values into the original file, keeping comments intact."""
    updated = text
    for key, value in overrides.items():
        for name in (f"{prefix}_{key}", key):
            pattern = rf"^(\s*{re.escape(name)}\s*=\s*)([^\n#]*)"
            if not re.search(pattern, updated, flags=re.MULTILINE):
                continue
            replacement = repr(value) if isinstance(value, str) else f"{value}"

            def substitute(match, text=replacement):
                # Keep the gap before a trailing comment so the file stays tidy.
                original = match.group(2)
                padding = " " * (len(original) - len(original.rstrip()))
                return match.group(1) + text + padding

            updated = re.sub(pattern, substitute, updated, count=1, flags=re.MULTILINE)
            break
    return updated


# ---------------------------------------------------------------------------
# Main.


def main() -> None:
    st.set_page_config(page_title="xtls-py", layout="wide")
    settings = sidebar()
    overrides = dict(settings["overrides"])
    overrides.update(QUALITY_PRESETS[settings["quality"]])

    state = st.session_state
    state.setdefault("pinned", [])

    top = st.container()
    with top:
        columns = st.columns([1, 1, 1, 3])
        run_now = columns[0].button("계산", type="primary", width='stretch')
        pin = columns[1].button("현재 곡선 고정", width='stretch')
        clear = columns[2].button("고정 해제", width='stretch')
    if clear:
        state["pinned"] = []

    if not (settings["auto"] or run_now or state.get("last")):
        st.info("좌측에서 파라미터를 조정한 뒤 **계산**을 누르세요.")
        return

    results = {}
    overrides_json = json.dumps(overrides, sort_keys=True)
    for mode in settings["modes"]:
        label = mode.upper()
        with st.spinner(f"{label} 계산 중…"):
            try:
                results[mode] = compute(mode, settings["input_file"], overrides_json)
            except Exception as error:  # noqa: BLE001 -- surfaced to the user
                st.error(f"{label} 계산 실패: {error}")
    if not results:
        return
    state["last"] = True

    if pin:
        summary = ", ".join(
            f"{key}={overrides[key]:g}" for key in ("delta", "u_charge_transfer", "ten_dq") if key in overrides
        )
        state["pinned"].append(
            {
                "label": summary,
                "curves": {mode: (energy.copy(), {k: v.copy() for k, v in curves.items()}) for mode, (energy, curves, _) in results.items()},
            }
        )

    experiment = None
    with st.expander("실험 데이터 비교", expanded=False):
        upload_columns = st.columns([2, 1, 1])
        upload = upload_columns[0].file_uploader("측정 스펙트럼 (2열 또는 헤더 있는 다열 텍스트)", type=["txt", "csv", "dat"])
        if upload is None:
            bundled = sorted((ROOT / "data").glob("*.txt"))
            if bundled:
                choice = upload_columns[0].selectbox(
                    "또는 data/ 에서 선택", ["(없음)"] + [path.name for path in bundled]
                )
                if choice != "(없음)":
                    upload = ROOT / "data" / choice
        loaded = load_experiment(upload)
        if loaded is not None:
            exp_energy, exp_matrix, exp_names = loaded
            column = upload_columns[1].selectbox("열", exp_names, index=min(len(exp_names) - 1, 0))
            shift = upload_columns[2].slider("에너지 이동 (eV)", -20.0, 20.0, 0.0, 0.05)
            experiment = (exp_energy + shift, exp_matrix[:, exp_names.index(column)], column)

    for mode in settings["modes"]:
        energy, curves, metadata = results[mode]
        st.subheader(f"{mode.upper()}  ·  {Path(settings['input_file']).stem}")
        default_columns = ["iso", "ld"] if mode == "xas" else ["total"]
        available = list(curves)
        chosen = st.multiselect(
            "표시할 곡선",
            available,
            default=[name for name in default_columns if name in available] or available[:1],
            key=f"cols_{mode}",
        )
        figure = plot(mode, energy, curves, chosen, state["pinned"], experiment)
        st.pyplot(figure, width='stretch')
        plt.close(figure)

        info = st.columns(4)
        info[0].metric("초기 기저", int(metadata["initial_basis_size"]))
        info[1].metric("최종 기저", int(metadata["final_basis_size"]))
        info[2].metric("onset (eV)", f"{metadata['onset']:.3f}")
        if mode == "xps":
            info[3].metric("합 규칙", f"{metadata['sum_rule_weight']:.3f}")
        if experiment is not None and chosen:
            deviations = [
                f"{name}: {residual(energy, curves[name], experiment[0], experiment[1]):.2f}%" for name in chosen
            ]
            st.caption("실험 대비 RMS 편차 (최대값 정규화 후) — " + ",  ".join(deviations))

        with st.expander("배치 에너지 / 상태 조성"):
            table = st.columns(2)
            table[0].dataframe(
                [
                    {"label": row["label"], "배치": row["configuration"], "E (eV)": round(row["energy_eV"], 4)}
                    for row in metadata["configuration_energies"]
                ],
                width='stretch',
                hide_index=True,
            )
            analysis = metadata.get("state_analysis", {})
            rows = analysis.get("initial", [])[:5]
            if rows and "composition" in rows[0]:
                table[1].dataframe(
                    [
                        {"상태": row["index"], "dE (meV)": round(row["energy_relative_meV"], 3), "조성": row["composition"]}
                        for row in rows
                    ],
                    width='stretch',
                    hide_index=True,
                )

    export(settings, overrides)


def plot(mode, energy, curves, chosen, pinned, experiment):
    figure, axis = plt.subplots(figsize=(9.0, 4.2))

    for index, entry in enumerate(pinned):
        if mode not in entry["curves"]:
            continue
        pin_energy, pin_curves = entry["curves"][mode]
        for name in chosen:
            if name in pin_curves:
                axis.plot(
                    pin_energy,
                    pin_curves[name],
                    lw=0.9,
                    alpha=0.45,
                    color=f"C{index % 10}",
                    label=f"고정 {index + 1}: {entry['label']}" if name == chosen[0] else None,
                )

    for name in chosen:
        axis.plot(energy, curves[name], lw=1.7, label=name)

    if experiment is not None:
        exp_energy, exp_values, exp_name = experiment
        reference = curves[chosen[0]] if chosen else next(iter(curves.values()))
        scale = (np.max(np.abs(reference)) or 1.0) / (np.max(np.abs(exp_values)) or 1.0)
        axis.plot(exp_energy, exp_values * scale, "k--", lw=1.0, alpha=0.75, label=f"실험: {exp_name}")

    axis.set_xlabel("Binding energy (eV)" if mode == "xps" else "Relative energy (eV)")
    axis.set_ylabel("Intensity (arb.)")
    axis.axhline(0.0, color="0.7", lw=0.6)
    if mode == "xps":
        axis.invert_xaxis()
    axis.legend(frameon=False, fontsize=8, ncol=2)
    axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout()
    return figure


def export(settings, overrides) -> None:
    st.divider()
    st.subheader("입력 파일로 저장")
    st.caption(
        "지금 슬라이더 값을 원본 입력 파일에 반영해 내보냅니다. 주석은 그대로 유지되고, "
        "`run.py` 로 최종 품질 재계산이 가능합니다. 품질 프리셋은 반영하지 않습니다."
    )

    source = Path(settings["input_file"])
    text = source.read_text(encoding="utf-8")
    tunable = {key: value for key, value in overrides.items() if key not in QUALITY_PRESETS[settings["quality"]]}
    updated = rewrite_input(text, tunable, prefix=settings["modes"][0])

    columns = st.columns([2, 1])
    name = columns[0].text_input("파일 이름", value=f"{source.stem}_fit.py")
    columns[1].download_button(
        "내려받기",
        data=updated.encode("utf-8"),
        file_name=name,
        mime="text/x-python",
        width='stretch',
    )
    if columns[1].button("inputs/ 에 저장", width='stretch'):
        target = ROOT / "inputs" / name
        if target.exists():
            st.warning(f"이미 있는 파일입니다: {target.name}. 이름을 바꾸세요.")
        else:
            target.write_text(updated, encoding="utf-8")
            st.success(f"저장했습니다: inputs/{target.name}")

    with st.expander("변경 내용 미리보기"):
        st.code("\n".join(f"{key} = {value}" for key, value in sorted(tunable.items())), language="python")


if __name__ == "__main__":
    main()
