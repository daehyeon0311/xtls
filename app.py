"""Input-card builder for the charge-transfer cluster calculations.

    streamlit run app.py

This writes input files; it does not run spectra. Every number shown updates
instantly because nothing here diagonalizes anything -- the panels are basis
sizes from binomial coefficients, sector energies from a closed formula, a 5x5
crystal-field matrix, and a table lookup. Set the parameters, save the card,
then run the calculation from the command line:

    python run.py inputs/<name>.py both

Loading an existing file as the starting point keeps its comments: only the
values you changed are substituted back in.
"""

from __future__ import annotations

import re
import sys
from math import comb
from pathlib import Path

import numpy as np
import streamlit as st


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
for path in (ROOT, SRC):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from xtls_py.engine import (  # noqa: E402
    available_appendix_a_3d,
    get_appendix_a_3d,
    sector_energy,
)
from xtls_py.geometry import crystal_field  # noqa: E402


D_ELECTRON_GROUP = {
    "K": 1, "Ca": 2, "Sc": 3, "Ti": 4, "V": 5, "Cr": 6,
    "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
}

GEOMETRIES = (
    "tetrahedral",
    "octahedral",
    "square_planar",
    "square_pyramidal",
    "trigonal_bipyramidal",
    "custom_xyz",
    "custom_spherical",
)

D_ORBITAL_LABELS = ("xy", "yz", "zx", "x²-y²", "3z²-r²")


# ---------------------------------------------------------------------------
# Parameter definitions.
#
# Each entry is (key, widget, label, options). `options` carries the widget's
# bounds or choices. Everything is rendered from this table, and the resulting
# values are what get substituted into the input card.

FIELDS: dict[str, list[tuple]] = {
    "기본": [
        ("case_name", "text", "케이스 이름", {}),
        ("ion", "text", "이온 (예: Fe2+)", {"help": "element, d 전자수, 출력 경로, 플롯 오프셋이 자동으로 정해집니다."}),
        ("max_ligand_holes", "int", "최대 ligand hole 수", {"min": 0, "max": 4}),
    ],
    "클러스터": [
        ("delta", "number", "Δ  전하이동 에너지 (eV)", {"step": 0.1}),
        ("u_charge_transfer", "number", "U_dd  d-d 쿨롱 (eV)", {"step": 0.1}),
        ("core_hole_potential", "number", "U_dc  core hole 인력 (eV)", {"step": 0.1}),
        ("ten_dq", "number", "10Dq  결정장 (eV)", {"step": 0.01}),
        ("ligand_ten_dq", "number", "10Dq(L)  리간드 분열 (eV)", {"step": 0.01, "help": "XTLS `10Dq(Ld) = 2*Tpp`. 0이면 리간드 궤도가 축퇴됩니다."}),
    ],
    "기하": [
        ("coordination_geometry", "select", "배위 기하", {"choices": GEOMETRIES}),
        ("ligand_radius", "number", "금속-리간드 거리 (Å)", {"step": 0.01}),
        ("ligand_angle_offset_deg", "number", "각도 오프셋 (deg)", {"step": 0.1, "help": "사면체 일그러짐. 0이면 이상적인 사면체."}),
        ("r2", "number", "r₂", {"step": 0.001}),
        ("r4", "number", "r₄", {"step": 0.001}),
    ],
    "혼성화": [
        ("hybridization_mode", "select", "방식", {"choices": ("geometry", "symmetry", "scalar")}),
        ("pd_sigma", "number", "pdσ  (geometry)", {"step": 0.05}),
        ("pd_ratio", "number", "pdπ/pdσ  (geometry)", {"step": 0.05}),
        ("d_ref", "number", "기준 거리 (geometry)", {"step": 0.05}),
        ("v_eg", "number", "V(e_g)  (symmetry)", {"step": 0.05}),
        ("v_t2g", "number", "V(t₂g)  (symmetry)", {"step": 0.05}),
        ("hopping", "number", "hopping  (scalar)", {"step": 0.05}),
    ],
    "다중항 축소": [
        ("fdd2_scale", "number", "F²dd", {"step": 0.005}),
        ("fdd4_scale", "number", "F⁴dd", {"step": 0.005}),
        ("fpd2_scale", "number", "F²pd", {"step": 0.005}),
        ("gpd1_scale", "number", "G¹pd", {"step": 0.005}),
        ("gpd3_scale", "number", "G³pd", {"step": 0.005}),
        ("so3d_scale", "number", "ζ₃d", {"step": 0.005}),
        ("so2p_scale", "number", "ζ₂p", {"step": 0.005}),
    ],
    "솔버": [
        ("n_initial_states", "int", "초기 상태 수", {"min": 1, "max": 100}),
        ("temperature_kelvin", "number", "온도 (K)", {"step": 10.0}),
        ("spectrum_method", "select", "방법", {"choices": ("lanczos", "exact")}),
        ("n_recursion", "int", "Lanczos 반복", {"min": 20, "max": 3000}),
        ("n_analyzed_states", "int", "분석 상태 수", {"min": 1, "max": 100}),
        ("energy_step", "number", "에너지 격자 (eV)", {"step": 0.005}),
        ("normalize", "select", "정규화", {"choices": ("max", "area", "none")}),
    ],
}

# Settings that differ between the two spectroscopies carry a prefix.
PREFIXED_FIELDS = [
    ("energy_min", "number", "에너지 최소 (eV)", {"step": 0.5}),
    ("energy_max", "number", "에너지 최대 (eV)", {"step": 0.5}),
    ("lorentzian_hwhm", "number", "Lorentzian HWHM (eV)", {"step": 0.01}),
    ("gaussian_hwhm", "number", "Gaussian HWHM (eV)", {"step": 0.01}),
    ("plot_relative_energy_min", "number", "플롯 최소 (eV)", {"step": 0.5}),
    ("plot_relative_energy_max", "number", "플롯 최대 (eV)", {"step": 0.5}),
]

XAS_ONLY = [
    ("make_experimental_geometry_curves", "bool", "실험 기하 곡선 생성", {}),
    ("grazing_angle_deg", "number", "입사각 (deg)", {"step": 0.5}),
    ("inplane_curve", "select", "면내 곡선", {"choices": ("ab", "x", "y")}),
    ("overlay_xtls", "bool", "XTLS 참조 겹치기", {}),
]

XPS_ONLY = [
    ("photoemission_shell", "select", "방출 껍질", {"choices": ("2p", "3d", "ligand")}),
    ("spin_resolved", "bool", "스핀 분해 출력", {}),
    ("plot_binding_energy_axis", "bool", "결합에너지 축 반전", {}),
]


# ---------------------------------------------------------------------------
# Reading a template.


@st.cache_data(show_spinner=False)
def read_template(path: str) -> tuple[str, dict]:
    """Return the file text and the values it defines, without side effects."""
    import runpy

    text = Path(path).read_text(encoding="utf-8")
    values = runpy.run_path(path, init_globals={"ROOT": ROOT, "Path": Path})
    values = {key: value for key, value in values.items() if not key.startswith("__")}
    return text, values


def resolved(values: dict, key: str, prefix: str, fallback=None):
    """Prefixed key wins, matching how the runners read an input file."""
    if f"{prefix}_{key}" in values:
        return values[f"{prefix}_{key}"]
    return values.get(key, fallback)


# ---------------------------------------------------------------------------
# Instant panels. None of this diagonalizes anything.


def parse_ion(text: str) -> tuple[str, int] | None:
    match = re.fullmatch(r"\s*([A-Za-z]{1,2})\s*([0-9]+)\s*\+\s*", str(text or ""))
    if not match:
        return None
    element = match.group(1).capitalize()
    if element not in D_ELECTRON_GROUP:
        return None
    return element, int(match.group(2))


def d_count(element: str, charge: int) -> int:
    return D_ELECTRON_GROUP[element] - charge


def effective_holes(element: str, n_d: int, requested: int, mode: str) -> tuple[int, str]:
    """How many ligand-hole sectors actually survive the Appendix-A table."""
    offset = {"xas": 1, "xps_core": 0, "xps_valence": -1}[mode]
    p_electrons = 6 if mode == "xps_valence" else 5
    effective = -1
    reason = "none"
    for holes in range(requested + 1):
        initial_d = n_d + holes
        final_d = n_d + offset + holes
        if initial_d > 10 or not 0 <= final_d <= 10:
            reason = f"3d 껍질이 h={holes} 에서 참"
            break
        if not has_entry(element, 6, initial_d) or not has_entry(element, p_electrons, final_d):
            reason = f"h={holes} 에 Appendix-A 항목 없음 (3d{initial_d} / 3d{final_d})"
            break
        effective = holes
    return max(effective, 0), reason


def has_entry(element: str, p_electrons: int, d_electrons: int) -> bool:
    if d_electrons in (0, 10):
        return True  # closed shell: filled with zeros, as XTLS X-cards do
    try:
        get_appendix_a_3d(element, p_electrons, d_electrons)
        return True
    except KeyError:
        return False


def basis_sizes(n_d: int, holes: int, mode: str) -> tuple[int, int]:
    offset = {"xas": 1, "xps_core": 0, "xps_valence": -1}[mode]
    core_factor = 1 if mode == "xps_valence" else 6
    initial = final = 0
    for h in range(holes + 1):
        if n_d + h <= 10:
            initial += comb(10, n_d + h) * comb(10, h)
        final_d = n_d + offset + h
        if 0 <= final_d <= 10:
            final += core_factor * comb(10, final_d) * comb(10, h)
    return initial, final


def runtime_hint(final_size: int) -> str:
    """Rough wall-clock guide, anchored on measured Fe d6 runs."""
    if final_size < 2000:
        return "1초 안팎"
    if final_size < 5000:
        return "수 초"
    if final_size < 12000:
        return "10초 안팎"
    if final_size < 30000:
        return "30초~1분"
    return "수 분 이상"


def geometry_positions(values: dict) -> np.ndarray | None:
    """Ligand sites as [r, theta, phi]; None when a custom list is required."""
    geometry = values.get("coordination_geometry")
    radius = float(values.get("ligand_radius") or 2.0)
    if geometry == "octahedral":
        return np.array([
            [radius, np.pi / 2, 0.0], [radius, np.pi / 2, np.pi],
            [radius, np.pi / 2, np.pi / 2], [radius, np.pi / 2, 3 * np.pi / 2],
            [radius, 0.0, 0.0], [radius, np.pi, 0.0],
        ])
    if geometry == "tetrahedral":
        angle = np.deg2rad(54.7356 + float(values.get("ligand_angle_offset_deg") or 0.0))
        return np.array([
            [radius, angle, np.deg2rad(-45)], [radius, angle, np.deg2rad(135)],
            [radius, np.pi - angle, np.deg2rad(45)], [radius, np.pi - angle, np.deg2rad(-135)],
        ])
    if geometry == "square_planar":
        return np.array([[radius, np.pi / 2, phi] for phi in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2)])
    if geometry == "square_pyramidal":
        return np.array(
            [[radius, np.pi / 2, phi] for phi in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2)]
            + [[radius, 0.0, 0.0]]
        )
    if geometry == "trigonal_bipyramidal":
        return np.array(
            [[radius, np.pi / 2, phi] for phi in (0.0, 2 * np.pi / 3, 4 * np.pi / 3)]
            + [[radius, 0.0, 0.0], [radius, np.pi, 0.0]]
        )
    return None


# ---------------------------------------------------------------------------
# Writing the card.


def rewrite(text: str, overrides: dict, prefixes: dict[str, str]) -> str:
    """Substitute values into the template, leaving comments and layout alone."""
    updated = text
    for key, value in overrides.items():
        candidates = [key]
        if key in prefixes:
            candidates.insert(0, f"{prefixes[key]}_{key}")
        for name in candidates:
            pattern = rf"^(\s*{re.escape(name)}\s*=\s*)([^\n#]*)"
            if not re.search(pattern, updated, flags=re.MULTILINE):
                continue
            literal = repr(value) if isinstance(value, str) else str(value)

            def substitute(match, text=literal):
                original = match.group(2)
                padding = " " * (len(original) - len(original.rstrip()))
                return match.group(1) + text + padding

            updated = re.sub(pattern, substitute, updated, count=1, flags=re.MULTILINE)
            break
    return updated


def widget(container, key: str, kind: str, label: str, options: dict, current):
    handle = f"w_{key}"
    help_text = options.get("help")
    if kind == "text":
        return container.text_input(label, value=str(current or ""), key=handle, help=help_text)
    if kind == "bool":
        return container.checkbox(label, value=bool(current), key=handle, help=help_text)
    if kind == "select":
        choices = list(options["choices"])
        index = choices.index(current) if current in choices else 0
        return container.selectbox(label, choices, index=index, key=handle, help=help_text)
    if kind == "int":
        return container.number_input(
            label,
            value=int(current if current is not None else 1),
            min_value=options.get("min", 0),
            max_value=options.get("max", 1000),
            step=1,
            key=handle,
            help=help_text,
        )
    return container.number_input(
        label,
        value=float(current if current is not None else 0.0),
        step=options.get("step", 0.1),
        format="%.4g",
        key=handle,
        help=help_text,
    )


# ---------------------------------------------------------------------------
# Main.


def main() -> None:
    st.set_page_config(page_title="xtls-py 입력 카드", layout="wide")
    st.title("입력 카드 만들기")
    st.caption(
        "여기서는 계산하지 않습니다. 파라미터를 정해 입력 파일을 저장한 뒤, "
        "`python run.py inputs/<이름>.py both` 로 실행하세요."
    )

    templates = sorted((ROOT / "inputs").glob("*.py"))
    if not templates:
        st.error("inputs/ 에 템플릿으로 쓸 입력 파일이 없습니다.")
        st.stop()

    header = st.columns([2, 1])
    template_path = header[0].selectbox(
        "시작 템플릿", templates, format_func=lambda path: path.stem,
        help="이 파일의 값과 주석을 그대로 물려받고, 바꾼 값만 치환합니다.",
    )
    text, base = read_template(str(template_path))

    values: dict = {}
    prefixes: dict[str, str] = {}

    tabs = st.tabs(list(FIELDS) + ["XAS 설정", "XPS 설정", "미리보기 · 저장"])

    for tab, (section, fields) in zip(tabs, FIELDS.items()):
        with tab:
            columns = st.columns(2)
            for index, (key, kind, label, options) in enumerate(fields):
                values[key] = widget(columns[index % 2], key, kind, label, options, base.get(key))

    with tabs[len(FIELDS)]:
        st.caption("XAS 전용 설정과, 두 분광법이 다르게 쓰는 값의 XAS 쪽입니다.")
        columns = st.columns(2)
        for index, (key, kind, label, options) in enumerate(PREFIXED_FIELDS):
            values[f"xas__{key}"] = widget(
                columns[index % 2], f"xas_{key}", kind, label, options, resolved(base, key, "xas")
            )
        st.divider()
        columns = st.columns(2)
        for index, (key, kind, label, options) in enumerate(XAS_ONLY):
            values[key] = widget(columns[index % 2], key, kind, label, options, base.get(key))

    with tabs[len(FIELDS) + 1]:
        st.caption("XPS 전용 설정과, 두 분광법이 다르게 쓰는 값의 XPS 쪽입니다.")
        columns = st.columns(2)
        for index, (key, kind, label, options) in enumerate(PREFIXED_FIELDS):
            values[f"xps__{key}"] = widget(
                columns[index % 2], f"xps_{key}", kind, label, options, resolved(base, key, "xps")
            )
        st.divider()
        columns = st.columns(2)
        for index, (key, kind, label, options) in enumerate(XPS_ONLY):
            values[key] = widget(columns[index % 2], key, kind, label, options, base.get(key))

    # Flatten the prefixed entries back into real key names.
    overrides: dict = {}
    for key, value in values.items():
        if "__" in key:
            prefix, name = key.split("__", 1)
            overrides[name] = value
            prefixes[name] = prefix
        else:
            overrides[key] = value

    # The prefixed keys appear twice, once per spectroscopy, so they cannot go
    # through the single-substitution rewrite together. Handle them separately.
    shared = {key: value for key, value in overrides.items() if key not in prefixes}
    updated = rewrite(text, shared, {})
    for spectroscopy in ("xas", "xps"):
        block = {
            name: values[f"{spectroscopy}__{name}"]
            for name in prefixes
            if f"{spectroscopy}__{name}" in values
        }
        updated = rewrite(updated, block, {name: spectroscopy for name in block})

    st.divider()
    summary(values, base)

    with tabs[-1]:
        preview(updated, values, template_path)


def summary(values: dict, base: dict) -> None:
    """Instant feedback: basis sizes, sector energies, crystal field, Slater."""
    ion = parse_ion(values.get("ion", ""))
    st.subheader("즉시 확인")

    if ion is None:
        st.warning(f"이온 표기를 인식하지 못했습니다: {values.get('ion')!r}. 예: `Fe2+`, `Ni2+`")
        return
    element, charge = ion
    n_d = d_count(element, charge)
    if not 0 <= n_d <= 10:
        st.error(f"{element}{charge}+ 는 d{n_d} 로, 0~10 범위를 벗어납니다.")
        return

    requested = int(values.get("max_ligand_holes", 1))
    xps_mode = "xps_valence" if values.get("photoemission_shell") == "3d" else "xps_core"

    columns = st.columns(3)
    columns[0].markdown(f"**{element}{charge}+ → 3d{n_d}**")

    rows = []
    for label, mode in (("XAS", "xas"), ("XPS", xps_mode)):
        holes, reason = effective_holes(element, n_d, requested, mode)
        initial, final = basis_sizes(n_d, holes, mode)
        rows.append(
            {
                "분광법": label,
                "유효 h": holes,
                "초기 기저": f"{initial:,}",
                "최종 기저": f"{final:,}",
                "예상 시간": runtime_hint(final),
                "비고": "" if holes >= requested else reason,
            }
        )
    columns[0].dataframe(rows, hide_index=True, width="stretch")

    # Sector energies -- a closed formula, no diagonalization.
    delta = float(values.get("delta", 0.0))
    u_dd = float(values.get("u_charge_transfer", 0.0))
    u_dc = float(values.get("core_hole_potential", 0.0))
    holes_xas, _ = effective_holes(element, n_d, requested, "xas")
    energy_rows = []
    for h in range(holes_xas + 1):
        energy_rows.append({
            "배치": f"2p⁶ 3d{n_d + h}" + ("" if h == 0 else f" L{h if h > 1 else ''}"),
            "종류": "초기",
            "E (eV)": round(sector_energy(h, 0, delta, u_dd, u_dc, d_electron_offset=0), 3),
        })
    for h in range(holes_xas + 1):
        energy_rows.append({
            "배치": f"2p⁵ 3d{n_d + 1 + h}" + ("" if h == 0 else f" L{h if h > 1 else ''}"),
            "종류": "XAS 종",
            "E (eV)": round(sector_energy(h, 1, delta, u_dd, u_dc, d_electron_offset=1), 3),
        })
    for h in range(holes_xas + 1):
        energy_rows.append({
            "배치": f"2p⁵ 3d{n_d + h}" + ("" if h == 0 else f" L{h if h > 1 else ''}"),
            "종류": "XPS 종",
            "E (eV)": round(sector_energy(h, 1, delta, u_dd, u_dc, d_electron_offset=0), 3),
        })
    columns[1].markdown("**배치 에너지**")
    columns[1].dataframe(energy_rows, hide_index=True, width="stretch", height=260)

    # Crystal field -- a 5x5 eigenvalue problem, microseconds.
    positions = geometry_positions(values)
    columns[2].markdown("**결정장 준위**")
    if positions is None:
        columns[2].info("사용자 지정 좌표는 입력 파일에서 직접 적어주세요.")
    else:
        try:
            matrix, *_ = crystal_field(positions, float(values.get("ten_dq", 0.0)),
                                       float(values.get("r2", 0.4)), float(values.get("r4", 0.38)))
            levels = np.real(np.diag(matrix))
            columns[2].dataframe(
                [{"궤도": name, "E (eV)": round(float(level), 4)}
                 for name, level in zip(D_ORBITAL_LABELS, levels)],
                hide_index=True, width="stretch",
            )
            spread = float(np.max(levels) - np.min(levels))
            columns[2].caption(f"전체 분열폭 {spread:.3f} eV")
        except Exception as error:  # noqa: BLE001 -- shown to the user
            columns[2].warning(f"결정장 계산 불가: {error}")

    # Slater integrals actually available for these sectors.
    with st.expander("이 배치에서 쓰이는 Slater 적분 (Appendix A)"):
        table = []
        wanted = [(6, n_d + h) for h in range(holes_xas + 1)]
        wanted += [(5, n_d + 1 + h) for h in range(holes_xas + 1)]
        wanted += [(5, n_d + h) for h in range(holes_xas + 1)]
        for p_electrons, d_electrons in sorted(set(wanted)):
            if not 0 <= d_electrons <= 10:
                continue
            try:
                entry = get_appendix_a_3d(element, p_electrons, d_electrons)
            except KeyError:
                table.append({
                    "배치": f"2p{p_electrons} 3d{d_electrons}",
                    "비고": "표 없음 — 닫힌 껍질이면 0으로 채워짐" if d_electrons in (0, 10) else "표 없음 (계산 불가)",
                })
                continue
            scale = float(values.get("fdd2_scale", 1.0))
            table.append({
                "배치": f"2p{p_electrons} 3d{d_electrons}",
                "F²dd": round(entry.fdd2 * scale, 3),
                "F⁴dd": round(entry.fdd4 * scale, 3),
                "ζ₃d": round(entry.zeta_d * float(values.get("so3d_scale", 1.0)), 4),
                "ζ₂p": None if entry.zeta_2p is None else round(entry.zeta_2p * float(values.get("so2p_scale", 1.0)), 3),
                "비고": "",
            })
        st.dataframe(table, hide_index=True, width="stretch")
        st.caption(
            f"{element} 에 대해 표에 있는 배치: "
            + ", ".join(f"2p{p}3d{d}" for _e, p, d in sorted(available_appendix_a_3d(element)))
        )


def preview(updated: str, values: dict, template_path: Path) -> None:
    st.subheader("저장")
    name_default = f"{values.get('case_name') or template_path.stem}.py"
    columns = st.columns([2, 1, 1])
    name = columns[0].text_input("파일 이름", value=name_default)
    if not name.endswith(".py"):
        name += ".py"

    columns[1].download_button(
        "내려받기", data=updated.encode("utf-8"), file_name=name,
        mime="text/x-python", width="stretch",
    )
    target = ROOT / "inputs" / name
    overwrite = columns[2].checkbox("덮어쓰기 허용", value=False)
    if columns[2].button("inputs/ 에 저장", width="stretch", type="primary"):
        if target.exists() and not overwrite:
            st.warning(f"이미 있는 파일입니다: inputs/{name}. 덮어쓰려면 위 체크박스를 켜세요.")
        else:
            target.write_text(updated, encoding="utf-8")
            st.success(f"저장했습니다: inputs/{name}")
            st.code(f"python run.py inputs/{name} both", language="bash")

    st.divider()
    st.caption("생성될 파일 전체입니다. 템플릿의 주석과 구조는 그대로 유지됩니다.")
    st.code(updated, language="python")


if __name__ == "__main__":
    main()
