#!/usr/bin/env python3
"""
=============================================================
 배터리/에너지 소재 데모 (Battery & Energy Materials Demo)
=============================================================
 Pymatgen 기반 결정 구조 생성/분석, 배터리 전압 계산,
 상평형도 분석, ASE 원자 시뮬레이션, PyBaMM 배터리 모델링
=============================================================
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore')

print("=" * 65)
print(" 🔋 배터리/에너지 소재 데모 — Battery Materials Pipeline")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# 1. 결정 구조 생성 및 분석 (Pymatgen)
# ─────────────────────────────────────────────────────────────
print("\n📦 Step 1: 결정 구조 생성 및 분석 (Pymatgen)")
print("-" * 50)

from pymatgen.core import Structure, Lattice, Element, Composition

# LiFePO4 (리튬인산철) — 대표적 배터리 양극재
print("\n  [1-1] LiFePO4 올리빈 구조 생성")
lattice = Lattice.from_parameters(
    a=10.332, b=6.010, c=4.692,
    alpha=90, beta=90, gamma=90
)
# 단순화된 LiFePO4 구조 (Pnma 공간군)
species = ["Li", "Fe", "P", "O", "O", "O", "O"]
coords = [
    [0.0, 0.0, 0.0],     # Li
    [0.282, 0.25, 0.975], # Fe
    [0.094, 0.25, 0.418], # P
    [0.097, 0.25, 0.743], # O1
    [0.457, 0.25, 0.206], # O2
    [0.165, 0.046, 0.285],# O3
    [0.165, 0.454, 0.285],# O4
]
lifepo4 = Structure(lattice, species, coords)

print(f"  화학식: {lifepo4.composition.reduced_formula}")
print(f"  격자 상수: a={lattice.a:.3f}, b={lattice.b:.3f}, c={lattice.c:.3f} Å")
print(f"  격자 부피: {lattice.volume:.2f} ų")
print(f"  원자 수: {len(lifepo4)}")
print(f"  밀도: {lifepo4.density:.4f} g/cm³")

# ─────────────────────────────────────────────────────────────
# 2. 배터리 양극재 비교 분석
# ─────────────────────────────────────────────────────────────
print("\n\n📊 Step 2: 배터리 양극재 후보 비교 분석")
print("-" * 50)

cathode_materials = {
    "LiFePO4": {
        "formula": "LiFePO4",
        "voltage": 3.4,  # V vs Li/Li+
        "capacity": 170, # mAh/g (이론)
        "cycle_life": 3000,
        "safety": "매우 우수",
        "cost": "낮음",
        "energy_density": 578,  # Wh/kg
    },
    "LiCoO2": {
        "formula": "LiCoO2",
        "voltage": 3.9,
        "capacity": 274,
        "cycle_life": 500,
        "safety": "보통",
        "cost": "높음",
        "energy_density": 1069,
    },
    "LiMn2O4": {
        "formula": "LiMn2O4",
        "voltage": 4.1,
        "capacity": 148,
        "cycle_life": 1500,
        "safety": "우수",
        "cost": "낮음",
        "energy_density": 607,
    },
    "NMC811": {
        "formula": "LiNi0.8Mn0.1Co0.1O2",
        "voltage": 3.7,
        "capacity": 200,
        "cycle_life": 800,
        "safety": "보통",
        "cost": "중간",
        "energy_density": 740,
    },
    "NCA": {
        "formula": "LiNi0.8Co0.15Al0.05O2",
        "voltage": 3.7,
        "capacity": 200,
        "cycle_life": 1000,
        "safety": "보통",
        "cost": "중간",
        "energy_density": 740,
    },
    "LiFeSO4F": {
        "formula": "LiFeSO4F",
        "voltage": 3.6,
        "capacity": 151,
        "cycle_life": 2000,
        "safety": "우수",
        "cost": "낮음",
        "energy_density": 544,
    },
}

print(f"  {'소재':10s} {'전압(V)':>8s} {'용량':>8s} {'에너지':>8s} {'수명':>8s} {'안전성':>8s} {'비용':>6s}")
print(f"  {'─'*10} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*6}")

for name, props in cathode_materials.items():
    print(f"  {name:10s} {props['voltage']:>7.1f}V {props['capacity']:>6d}mAh {props['energy_density']:>6d}Wh {props['cycle_life']:>6d}회 {props['safety']:>8s} {props['cost']:>6s}")

# ─────────────────────────────────────────────────────────────
# 3. 원소 분석 및 자원 가용성 평가
# ─────────────────────────────────────────────────────────────
print("\n\n🌍 Step 3: 원소 분석 및 자원 가용성 평가")
print("-" * 50)

critical_elements = {
    "Li":  {"price_kg": 25.0,  "abundance": 20,    "criticality": "높음"},
    "Co":  {"price_kg": 30.0,  "abundance": 25,    "criticality": "매우 높음"},
    "Ni":  {"price_kg": 16.0,  "abundance": 84,    "criticality": "높음"},
    "Mn":  {"price_kg": 2.0,   "abundance": 950,   "criticality": "낮음"},
    "Fe":  {"price_kg": 0.08,  "abundance": 56300, "criticality": "매우 낮음"},
    "Al":  {"price_kg": 2.5,   "abundance": 81300, "criticality": "매우 낮음"},
    "P":   {"price_kg": 0.07,  "abundance": 1050,  "criticality": "낮음"},
}

print(f"  {'원소':>4s} {'가격($/kg)':>12s} {'지각존재량(ppm)':>16s} {'공급위험':>12s}")
print(f"  {'─'*4} {'─'*12} {'─'*16} {'─'*12}")

for elem, info in critical_elements.items():
    bar = "█" * min(int(np.log10(info["abundance"]+1) * 4), 20)
    print(f"  {elem:>4s} {info['price_kg']:>11.2f}$ {info['abundance']:>14d} {info['criticality']:>12s} {bar}")

# ─────────────────────────────────────────────────────────────
# 4. 조성 기반 이론 용량 계산
# ─────────────────────────────────────────────────────────────
print("\n\n⚡ Step 4: 이론 용량 계산 (Theoretical Capacity)")
print("-" * 50)

from pymatgen.core import Composition

def calc_theoretical_capacity(formula, n_electrons=1):
    """
    이론 비용량(mAh/g) = n * F / (3.6 * MW)
    n: 전달 전자 수, F: 패러데이 상수 (96485 C/mol)
    """
    comp = Composition(formula)
    mw = comp.weight  # g/mol
    F = 96485  # C/mol
    capacity = n_electrons * F / (3.6 * mw)
    return capacity, mw

print(f"  {'소재':20s} {'화학식':>20s} {'분자량':>10s} {'이론용량':>12s}")
print(f"  {'─'*20} {'─'*20} {'─'*10} {'─'*12}")

calc_targets = [
    ("LiFePO4", "LiFePO4", 1),
    ("LiCoO2", "LiCoO2", 1),
    ("LiMn2O4", "LiMn2O4", 1),
    ("LiFeSO4F", "LiFeSO4F", 1),
    ("Li2FeSiO4 (2e)", "Li2FeSiO4", 2),
    ("Na3V2(PO4)3", "Na3V2P3O12", 2),
    ("LiTiS2", "LiTiS2", 1),
]

for name, formula, n_e in calc_targets:
    cap, mw = calc_theoretical_capacity(formula, n_e)
    bar = "█" * int(cap / 15)
    print(f"  {name:20s} {formula:>20s} {mw:>9.2f}g {cap:>10.1f}mAh/g {bar}")

# ─────────────────────────────────────────────────────────────
# 5. ASE — 원자 시뮬레이션 (간단한 벌크 금속)
# ─────────────────────────────────────────────────────────────
print("\n\n🔬 Step 5: ASE 원자 시뮬레이션 (결정 구조)")
print("-" * 50)

from ase import Atoms
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.eos import EquationOfState
import io, sys

# 5-1: 벌크 금속 격자 상수 및 에너지 계산
print("\n  [5-1] 벌크 금속 에너지 계산 (EMT 포텐셜)")
metals = ["Cu", "Au", "Ag", "Pt", "Pd", "Ni", "Al"]
print(f"  {'금속':>4s} {'구조':>5s} {'격자상수(Å)':>12s} {'에너지(eV)':>12s} {'부피(ų)':>12s}")
print(f"  {'─'*4} {'─'*5} {'─'*12} {'─'*12} {'─'*12}")

for metal in metals:
    atoms = bulk(metal, 'fcc')
    atoms.calc = EMT()
    energy = atoms.get_potential_energy()
    vol = atoms.get_volume()
    a = atoms.cell.cellpar()[0]
    print(f"  {metal:>4s} {'FCC':>5s} {a:>11.4f} {energy:>11.4f} {vol:>11.4f}")

# 5-2: 상태방정식 (Equation of State)으로 평형 격자 상수 결정
print(f"\n  [5-2] Cu의 상태방정식 (Equation of State)")
cu = bulk('Cu', 'fcc')
cu.calc = EMT()

volumes = []
energies = []
a0 = cu.cell.cellpar()[0]

for scale in np.linspace(0.95, 1.05, 11):
    atoms = cu.copy()
    atoms.set_cell(cu.cell * scale, scale_atoms=True)
    atoms.calc = EMT()
    vol = atoms.get_volume()
    energy = atoms.get_potential_energy()
    volumes.append(vol)
    energies.append(energy)

# Birch-Murnaghan 상태방정식 피팅
eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()
a_eq = (4 * v0) ** (1/3)  # FCC: 4 atoms per unit cell → a = (4V)^(1/3)

print(f"  평형 부피: {v0:.4f} ų")
print(f"  평형 에너지: {e0:.6f} eV")
print(f"  체적 탄성률: {B / 1e9 * 160.2:.1f} GPa")  # eV/ų → GPa
print(f"  평형 격자 상수: {a_eq:.4f} Å")

# 5-3: 간단한 분자 에너지 계산
print(f"\n  [5-3] 간단한 분자 구조 생성")
molecules = ["H2", "H2O", "NH3", "CH4"]
for mol_name in molecules:
    try:
        mol = molecule(mol_name)
        positions = mol.get_positions()
        n_atoms = len(mol)
        symbols = mol.get_chemical_formula()
        print(f"  {mol_name:>5s} ({symbols:>5s}): 원자 {n_atoms}개, 위치 범위 [{positions.min():.2f}, {positions.max():.2f}] Å")
    except Exception as e:
        print(f"  {mol_name:>5s}: 생성 실패 ({e})")

# ─────────────────────────────────────────────────────────────
# 6. Pymatgen — 전기화학 분석
# ─────────────────────────────────────────────────────────────
print("\n\n⚡ Step 6: 전기화학 분석 (Pymatgen)")
print("-" * 50)

from pymatgen.core import Element

# 원소별 전기화학 특성
print(f"\n  [6-1] 배터리 관련 원소 특성")
battery_elements = ["Li", "Na", "K", "Mg", "Ca", "Al", "Zn"]
print(f"  {'원소':>4s} {'원자량':>8s} {'이온화E(eV)':>12s} {'전기음성도':>12s} {'이온반경':>10s}")
print(f"  {'─'*4} {'─'*8} {'─'*12} {'─'*12} {'─'*10}")

for elem_str in battery_elements:
    elem = Element(elem_str)
    ie = elem.ionization_energies[0] if elem.ionization_energies else 0
    en = elem.X if elem.X else 0
    ar = elem.atomic_radius if elem.atomic_radius else 0
    print(f"  {elem_str:>4s} {elem.atomic_mass:>7.2f}u {ie:>10.2f} {en:>11.2f} {ar:>9.2f}Å")

# ─────────────────────────────────────────────────────────────
# 7. PyBaMM — 리튬이온 배터리 시뮬레이션
# ─────────────────────────────────────────────────────────────
print("\n\n🔋 Step 7: PyBaMM 리튬이온 배터리 시뮬레이션")
print("-" * 50)

import pybamm

# SPM (Single Particle Model) — 단순하고 빠른 모델
print("\n  [7-1] SPM (Single Particle Model) — 1C 방전")
model_spm = pybamm.lithium_ion.SPM()
sim_spm = pybamm.Simulation(model_spm)

# 1C 방전 시뮬레이션 (3600초)
solution_spm = sim_spm.solve([0, 3600])

# 결과 추출
time_spm = solution_spm["Time [s]"].entries
voltage_spm = solution_spm["Voltage [V]"].entries
current_spm = solution_spm["Current [A]"].entries

print(f"  시뮬레이션 시간: {time_spm[-1]:.0f}초 ({time_spm[-1]/60:.1f}분)")
print(f"  초기 전압: {voltage_spm[0]:.4f} V")
print(f"  최종 전압: {voltage_spm[-1]:.4f} V")
print(f"  전압 강하: {voltage_spm[0] - voltage_spm[-1]:.4f} V")
print(f"  평균 전압: {np.mean(voltage_spm):.4f} V")
print(f"  방전 전류: {current_spm[0]:.4f} A")

# 전압 프로파일 텍스트 시각화
print(f"\n  [전압 프로파일 — 1C 방전]")
n_points = 20
indices = np.linspace(0, len(time_spm)-1, n_points, dtype=int)
v_min, v_max = voltage_spm.min(), voltage_spm.max()

for idx in indices:
    t = time_spm[idx]
    v = voltage_spm[idx]
    bar_len = int((v - v_min) / (v_max - v_min) * 40)
    bar = "█" * bar_len
    print(f"  {t/60:6.1f}min | {v:.3f}V | {bar}")

# 7-2: 다른 C-rate 비교
print(f"\n  [7-2] C-rate별 방전 비교")
c_rates = [0.5, 1.0, 2.0, 3.0]
print(f"  {'C-rate':>7s} {'초기V':>7s} {'최종V':>7s} {'평균V':>7s} {'용량유지':>10s}")
print(f"  {'─'*7} {'─'*7} {'─'*7} {'─'*7} {'─'*10}")

ref_capacity = None
for c_rate in c_rates:
    model = pybamm.lithium_ion.SPM()
    param = model.default_parameter_values
    param["Current function [A]"] = param["Nominal cell capacity [A.h]"] * c_rate
    
    sim = pybamm.Simulation(model, parameter_values=param)
    try:
        sol = sim.solve([0, 3600 / c_rate])
        v = sol["Voltage [V]"].entries
        t = sol["Time [s]"].entries
        
        # 방전 용량 추정
        capacity = np.trapz(np.abs(sol["Current [A]"].entries), t / 3600)
        if ref_capacity is None:
            ref_capacity = capacity
        retention = (capacity / ref_capacity) * 100
        
        print(f"  {c_rate:>6.1f}C {v[0]:>7.3f} {v[-1]:>7.3f} {np.mean(v):>7.3f} {retention:>8.1f}%")
    except Exception as e:
        print(f"  {c_rate:>6.1f}C  시뮬레이션 실패: {str(e)[:40]}")

# 7-3: DFN 모델 (더 정밀한 모델)
print(f"\n  [7-3] DFN (Doyle-Fuller-Newman) 모델 — 1C 방전")
model_dfn = pybamm.lithium_ion.DFN()
sim_dfn = pybamm.Simulation(model_dfn)
solution_dfn = sim_dfn.solve([0, 3600])

voltage_dfn = solution_dfn["Voltage [V]"].entries
time_dfn = solution_dfn["Time [s]"].entries

print(f"  초기 전압: {voltage_dfn[0]:.4f} V")
print(f"  최종 전압: {voltage_dfn[-1]:.4f} V")
print(f"  평균 전압: {np.mean(voltage_dfn):.4f} V")

# SPM vs DFN 비교
print(f"\n  [SPM vs DFN 비교]")
print(f"  {'모델':>6s} {'초기V':>7s} {'최종V':>7s} {'평균V':>7s} {'전압차':>7s}")
print(f"  {'─'*6} {'─'*7} {'─'*7} {'─'*7} {'─'*7}")
print(f"  {'SPM':>6s} {voltage_spm[0]:>7.3f} {voltage_spm[-1]:>7.3f} {np.mean(voltage_spm):>7.3f} {voltage_spm[0]-voltage_spm[-1]:>7.3f}")
print(f"  {'DFN':>6s} {voltage_dfn[0]:>7.3f} {voltage_dfn[-1]:>7.3f} {np.mean(voltage_dfn):>7.3f} {voltage_dfn[0]-voltage_dfn[-1]:>7.3f}")

# ─────────────────────────────────────────────────────────────
# 8. 열역학적 안정성 분석
# ─────────────────────────────────────────────────────────────
print("\n\n🌡️ Step 8: 열역학적 안정성 분석")
print("-" * 50)

# 간이 열안정성 평가 (분해 온도 데이터 기반)
thermal_data = {
    "LiFePO4":       {"decomp_temp": 400, "onset_exotherm": 310, "heat_release": 0.1},
    "LiCoO2":        {"decomp_temp": 200, "onset_exotherm": 150, "heat_release": 1.0},
    "LiMn2O4":       {"decomp_temp": 280, "onset_exotherm": 250, "heat_release": 0.3},
    "NMC811":        {"decomp_temp": 210, "onset_exotherm": 170, "heat_release": 0.9},
    "LiNi0.5Mn1.5O4":{"decomp_temp": 280, "onset_exotherm": 240, "heat_release": 0.4},
}

print(f"  {'소재':18s} {'분해온도(°C)':>12s} {'발열개시(°C)':>12s} {'발열량':>8s} {'안전등급':>10s}")
print(f"  {'─'*18} {'─'*12} {'─'*12} {'─'*8} {'─'*10}")

for name, data in thermal_data.items():
    safety = "🟢 안전" if data["decomp_temp"] > 300 else ("🟡 주의" if data["decomp_temp"] > 250 else "🔴 위험")
    bar = "█" * int(data["decomp_temp"] / 20)
    print(f"  {name:18s} {data['decomp_temp']:>10d}°C {data['onset_exotherm']:>10d}°C {data['heat_release']:>7.1f}J {safety:>10s} {bar}")

# ─────────────────────────────────────────────────────────────
# 9. 종합 평가 — 차세대 배터리 소재 랭킹
# ─────────────────────────────────────────────────────────────
print("\n\n📋 Step 9: 종합 평가 — 배터리 소재 랭킹")
print("=" * 65)

def score_material(name, props, thermal):
    """가중 점수 계산"""
    # 정규화 점수 (0~1)
    energy_score = min(props["energy_density"] / 1100, 1.0)
    cycle_score = min(props["cycle_life"] / 3000, 1.0)
    safety_score = min(thermal["decomp_temp"] / 400, 1.0) if thermal else 0.5
    cost_map = {"낮음": 1.0, "중간": 0.6, "높음": 0.3}
    cost_score = cost_map.get(props["cost"], 0.5)
    
    # 가중 합 (에너지 30%, 수명 25%, 안전성 25%, 비용 20%)
    total = energy_score * 0.30 + cycle_score * 0.25 + safety_score * 0.25 + cost_score * 0.20
    return total, energy_score, cycle_score, safety_score, cost_score

print(f"  {'소재':10s} {'에너지':>8s} {'수명':>6s} {'안전':>6s} {'비용':>6s} {'종합':>8s}")
print(f"  {'─'*10} {'─'*8} {'─'*6} {'─'*6} {'─'*6} {'─'*8}")

rankings = []
for name, props in cathode_materials.items():
    thermal = thermal_data.get(name, None)
    total, e_s, c_s, s_s, co_s = score_material(name, props, thermal)
    rankings.append((name, total, e_s, c_s, s_s, co_s))

rankings.sort(key=lambda x: x[1], reverse=True)
for rank, (name, total, e_s, c_s, s_s, co_s) in enumerate(rankings, 1):
    grade = "🥇" if rank == 1 else ("🥈" if rank == 2 else ("🥉" if rank == 3 else "  "))
    bar = "█" * int(total * 20)
    print(f"  {name:10s} {e_s:>7.2f} {c_s:>6.2f} {s_s:>6.2f} {co_s:>6.2f}  {grade} {total:.3f} {bar}")

print("\n" + "=" * 65)
print(" 🔋 배터리/에너지 소재 데모 완료!")
print("=" * 65)
