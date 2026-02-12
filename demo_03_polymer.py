#!/usr/bin/env python3
"""
=============================================================
 폴리머/고분자 데모 (Polymer Science Demo)
=============================================================
 RDKit 기반 단량체 분석, 폴리머 물성 예측,
 중합 반응 시뮬레이션, 고분자 특성 평가
=============================================================
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Crippen, rdMolDescriptors

print("=" * 65)
print(" 🧪 폴리머/고분자 데모 — Polymer Science Pipeline")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# 1. 단량체(Monomer) 라이브러리 구축
# ─────────────────────────────────────────────────────────────
print("\n📦 Step 1: 단량체(Monomer) 라이브러리 구축")
print("-" * 50)

monomer_library = {
    "Ethylene":        {"smiles": "C=C",            "polymer": "PE (폴리에틸렌)",      "type": "비닐"},
    "Propylene":       {"smiles": "CC=C",           "polymer": "PP (폴리프로필렌)",    "type": "비닐"},
    "Styrene":         {"smiles": "C=Cc1ccccc1",    "polymer": "PS (폴리스타이렌)",    "type": "비닐"},
    "Vinyl_Chloride":  {"smiles": "C=CCl",          "polymer": "PVC (폴리염화비닐)",   "type": "비닐"},
    "Methyl_Methacrylate": {"smiles": "CC(=C)C(=O)OC", "polymer": "PMMA (아크릴)",   "type": "아크릴"},
    "Acrylonitrile":   {"smiles": "C=CC#N",         "polymer": "PAN (폴리아크릴로니트릴)", "type": "비닐"},
    "Tetrafluoroethylene": {"smiles": "FC(F)=C(F)F", "polymer": "PTFE (테프론)",      "type": "불소"},
    "Caprolactam":     {"smiles": "O=C1CCCCCN1",    "polymer": "Nylon 6",            "type": "아미드"},
    "Ethylene_Oxide":  {"smiles": "C1CO1",          "polymer": "PEO (폴리에틸렌옥사이드)", "type": "에테르"},
    "Lactic_Acid":     {"smiles": "CC(O)C(=O)O",    "polymer": "PLA (폴리유산)",      "type": "에스터"},
    "Bisphenol_A":     {"smiles": "CC(C)(C1=CC=C(O)C=C1)C1=CC=C(O)C=C1", "polymer": "PC (폴리카보네이트)", "type": "카보네이트"},
    "Dimethylsiloxane": {"smiles": "C[Si](C)(O)O",  "polymer": "PDMS (실리콘)",       "type": "실록산"},
}

print(f"  {'단량체':20s} {'SMILES':>25s} {'폴리머':>25s} {'유형':>8s}")
print(f"  {'─'*20} {'─'*25} {'─'*25} {'─'*8}")

for name, info in monomer_library.items():
    mol = Chem.MolFromSmiles(info["smiles"])
    if mol:
        info["mol"] = mol
        info["valid"] = True
        status = "✅"
    else:
        info["mol"] = None
        info["valid"] = False
        status = "❌"
    print(f"  {status} {name:18s} {info['smiles']:>25s} {info['polymer']:>25s} {info['type']:>8s}")

valid_monomers = {k: v for k, v in monomer_library.items() if v.get("valid")}
print(f"\n  총 {len(valid_monomers)}/{len(monomer_library)}개 단량체 파싱 성공")

# ─────────────────────────────────────────────────────────────
# 2. 단량체 물성 분석
# ─────────────────────────────────────────────────────────────
print("\n\n📊 Step 2: 단량체 물성 분석")
print("-" * 50)

print(f"  {'단량체':18s} {'분자량':>8s} {'LogP':>7s} {'TPSA':>7s} {'HBD':>5s} {'HBA':>5s} {'이중결합':>8s}")
print(f"  {'─'*18} {'─'*8} {'─'*7} {'─'*7} {'─'*5} {'─'*5} {'─'*8}")

monomer_props = {}
for name, info in valid_monomers.items():
    mol = info["mol"]
    props = {
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "DoubleBonds": sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2.0),
    }
    monomer_props[name] = props
    print(f"  {name:18s} {props['MW']:>8.2f} {props['LogP']:>7.2f} {props['TPSA']:>7.2f} {props['HBD']:>5d} {props['HBA']:>5d} {props['DoubleBonds']:>8d}")

# ─────────────────────────────────────────────────────────────
# 3. 폴리머 물성 예측 (경험적 그룹 기여법)
# ─────────────────────────────────────────────────────────────
print("\n\n🔬 Step 3: 폴리머 물성 예측 (Group Contribution Method)")
print("-" * 50)

# 실험값 기반 데이터베이스
polymer_properties = {
    "PE":    {"Tg": -125, "Tm": 137, "density": 0.95,  "tensile_MPa": 30,  "thermal_cond": 0.46},
    "PP":    {"Tg": -10,  "Tm": 170, "density": 0.90,  "tensile_MPa": 35,  "thermal_cond": 0.22},
    "PS":    {"Tg": 100,  "Tm": 240, "density": 1.05,  "tensile_MPa": 45,  "thermal_cond": 0.14},
    "PVC":   {"Tg": 82,   "Tm": 210, "density": 1.40,  "tensile_MPa": 55,  "thermal_cond": 0.16},
    "PMMA":  {"Tg": 105,  "Tm": 160, "density": 1.18,  "tensile_MPa": 70,  "thermal_cond": 0.19},
    "PAN":   {"Tg": 95,   "Tm": 317, "density": 1.18,  "tensile_MPa": 60,  "thermal_cond": 0.26},
    "PTFE":  {"Tg": -97,  "Tm": 327, "density": 2.15,  "tensile_MPa": 25,  "thermal_cond": 0.25},
    "Nylon6": {"Tg": 47,  "Tm": 225, "density": 1.14,  "tensile_MPa": 80,  "thermal_cond": 0.25},
    "PEO":   {"Tg": -67,  "Tm": 65,  "density": 1.13,  "tensile_MPa": 15,  "thermal_cond": 0.21},
    "PLA":   {"Tg": 60,   "Tm": 175, "density": 1.25,  "tensile_MPa": 50,  "thermal_cond": 0.13},
    "PC":    {"Tg": 150,  "Tm": 267, "density": 1.20,  "tensile_MPa": 65,  "thermal_cond": 0.20},
    "PDMS":  {"Tg": -127, "Tm": -40, "density": 0.97,  "tensile_MPa": 6,   "thermal_cond": 0.15},
}

print(f"  {'폴리머':8s} {'Tg(°C)':>8s} {'Tm(°C)':>8s} {'밀도':>8s} {'인장강도':>10s} {'열전도율':>10s}")
print(f"  {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*10} {'─'*10}")

for poly, props in polymer_properties.items():
    print(f"  {poly:8s} {props['Tg']:>7d}° {props['Tm']:>7d}° {props['density']:>7.2f} {props['tensile_MPa']:>8d}MPa {props['thermal_cond']:>8.2f}W/mK")

# ─────────────────────────────────────────────────────────────
# 4. Tg 예측 모델 (Fox-Flory 방정식)
# ─────────────────────────────────────────────────────────────
print("\n\n📈 Step 4: Tg 분자량 의존성 (Fox-Flory 방정식)")
print("-" * 50)
print("  Tg(M) = Tg(∞) - K / Mn")
print("  (M→∞일 때 Tg가 포화됨)")

# PS (폴리스타이렌) 예시
Tg_inf = 373  # K (100°C)
K_ff = 1.2e5  # Fox-Flory 상수 (K·g/mol)

print(f"\n  폴리스타이렌 (PS):")
print(f"  {'분자량(Da)':>12s} {'Tg(K)':>8s} {'Tg(°C)':>8s} {'Tg/Tg∞':>8s}")
print(f"  {'─'*12} {'─'*8} {'─'*8} {'─'*8}")

mol_weights = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 500000, 1000000]
for Mn in mol_weights:
    Tg = Tg_inf - K_ff / Mn
    Tg_C = Tg - 273.15
    ratio = Tg / Tg_inf
    bar = "█" * int(ratio * 30)
    print(f"  {Mn:>10,d} {Tg:>8.1f} {Tg_C:>7.1f}° {ratio:>7.3f} {bar}")

# ─────────────────────────────────────────────────────────────
# 5. 공중합체 Tg 예측 (Fox 방정식)
# ─────────────────────────────────────────────────────────────
print("\n\n🔀 Step 5: 공중합체 Tg 예측 (Fox Equation)")
print("-" * 50)
print("  1/Tg = w1/Tg1 + w2/Tg2  (랜덤 공중합체)")

# PS-PMMA 공중합체
Tg1 = 373  # PS Tg (K)
Tg2 = 378  # PMMA Tg (K)

print(f"\n  PS/PMMA 랜덤 공중합체:")
print(f"  {'PS wt%':>8s} {'PMMA wt%':>9s} {'Tg(K)':>7s} {'Tg(°C)':>8s}")
print(f"  {'─'*8} {'─'*9} {'─'*7} {'─'*8}")

for ps_frac in np.arange(0, 1.05, 0.1):
    pmma_frac = 1.0 - ps_frac
    if ps_frac == 0:
        Tg = Tg2
    elif pmma_frac == 0:
        Tg = Tg1
    else:
        Tg = 1.0 / (ps_frac / Tg1 + pmma_frac / Tg2)
    print(f"  {ps_frac*100:>6.0f}% {pmma_frac*100:>7.0f}%  {Tg:>6.1f} {Tg-273.15:>7.1f}°")

# 더 흥미로운 조합: PS/PDMS (큰 Tg 차이)
Tg_ps = 373   # K
Tg_pdms = 146 # K (-127°C)

print(f"\n  PS/PDMS 공중합체 (큰 Tg 차이):")
print(f"  {'PS wt%':>8s} {'PDMS wt%':>10s} {'Tg(°C)':>8s} {'그래프':>20s}")
print(f"  {'─'*8} {'─'*10} {'─'*8} {'─'*20}")

for ps_frac in np.arange(0, 1.05, 0.1):
    pdms_frac = 1.0 - ps_frac
    if ps_frac == 0:
        Tg = Tg_pdms
    elif pdms_frac == 0:
        Tg = Tg_ps
    else:
        Tg = 1.0 / (ps_frac / Tg_ps + pdms_frac / Tg_pdms)
    Tg_C = Tg - 273.15
    bar_pos = int((Tg_C + 130) / 230 * 30)
    bar = " " * max(bar_pos, 0) + "█"
    print(f"  {ps_frac*100:>6.0f}% {pdms_frac*100:>8.0f}%  {Tg_C:>7.1f}° {bar}")

# ─────────────────────────────────────────────────────────────
# 6. 중합도(DP) 및 분자량 분포
# ─────────────────────────────────────────────────────────────
print("\n\n📊 Step 6: 중합도 및 분자량 분포 시뮬레이션")
print("-" * 50)

# 축합 중합 (Flory 분포) 시뮬레이션
print("\n  [6-1] 축합 중합 — Flory 분포")
print("  P(x) = p^(x-1) * (1-p),  p = 전환율")

for p in [0.90, 0.95, 0.99, 0.995, 0.999]:
    Xn = 1 / (1 - p)                  # 수평균 중합도
    Xw = (1 + p) / (1 - p)            # 중량평균 중합도
    PDI = Xw / Xn                     # 다분산지수
    print(f"  p={p:.3f}: Xn={Xn:>8.1f}, Xw={Xw:>8.1f}, PDI={PDI:.3f}")

# 자유 라디칼 중합 — Schulz-Flory 분포 시뮬레이션
print(f"\n  [6-2] 자유 라디칼 중합 — 분자량 분포 시뮬레이션")
np.random.seed(42)
n_chains = 10000
avg_dp = 500  # 평균 중합도

# 지수 분포로 체인 길이 생성 (종결: 불균등화)
chain_lengths = np.random.exponential(avg_dp, n_chains).astype(int)
chain_lengths = chain_lengths[chain_lengths > 0]

# 단량체 분자량 (스타이렌: 104 g/mol)
M0 = 104
mol_weights_dist = chain_lengths * M0

Mn = np.mean(mol_weights_dist)
Mw = np.sum(mol_weights_dist**2) / np.sum(mol_weights_dist)
PDI = Mw / Mn

print(f"  체인 수: {len(chain_lengths):,d}")
print(f"  Mn (수평균 분자량): {Mn:,.0f} g/mol")
print(f"  Mw (중량평균 분자량): {Mw:,.0f} g/mol")
print(f"  PDI (다분산지수): {PDI:.3f}")
print(f"  Xn (수평균 중합도): {Mn/M0:.0f}")

# 분자량 분포 히스토그램 (텍스트)
print(f"\n  [분자량 분포 히스토그램]")
hist, bin_edges = np.histogram(mol_weights_dist / 1000, bins=15)
max_count = max(hist)
for i in range(len(hist)):
    low = bin_edges[i]
    high = bin_edges[i+1]
    bar_len = int(hist[i] / max_count * 40)
    bar = "█" * bar_len
    print(f"  {low:>6.0f}-{high:>5.0f}kDa | {bar} ({hist[i]:>4d})")

# ─────────────────────────────────────────────────────────────
# 7. 용해도 파라미터 (Solubility Parameter)
# ─────────────────────────────────────────────────────────────
print("\n\n💧 Step 7: Hildebrand 용해도 파라미터")
print("-" * 50)
print("  |δ1 - δ2| < 2 (MPa^0.5) → 상용성 우수")

solubility_params = {
    "PE":     16.2,
    "PP":     16.6,
    "PS":     18.5,
    "PVC":    19.5,
    "PMMA":   18.6,
    "PC":     19.4,
    "Nylon6": 22.5,
    "PLA":    20.2,
    "PDMS":   15.1,
    "PEO":    20.2,
}

# 용해도 파라미터 시각화
print(f"  {'폴리머':8s} {'δ(MPa^0.5)':>12s} {'스케일':>30s}")
print(f"  {'─'*8} {'─'*12} {'─'*30}")

for poly, delta in sorted(solubility_params.items(), key=lambda x: x[1]):
    bar_pos = int((delta - 14) / 10 * 40)
    bar = " " * bar_pos + "█"
    print(f"  {poly:8s} {delta:>10.1f}   {bar}")

# 상용성 매트릭스 (상위 6개)
print(f"\n  [상용성 매트릭스 — |δ1 - δ2|]")
top_polymers = list(solubility_params.keys())[:6]
print(f"  {'':8s}", end="")
for p in top_polymers:
    print(f" {p:>6s}", end="")
print()

for p1 in top_polymers:
    print(f"  {p1:8s}", end="")
    for p2 in top_polymers:
        diff = abs(solubility_params[p1] - solubility_params[p2])
        marker = "  ●  " if diff < 2 else f"{diff:>5.1f}"
        print(f" {marker:>6s}", end="")
    print()

print(f"  (● = 상용 가능, |δ1-δ2| < 2)")

# ─────────────────────────────────────────────────────────────
# 8. 열분해 안정성 분석
# ─────────────────────────────────────────────────────────────
print("\n\n🌡️ Step 8: 열분해 안정성 분석")
print("-" * 50)

thermal_stability = {
    "PTFE":   {"Td_5pct": 508, "Td_50pct": 580, "char_yield": 0},
    "PI":     {"Td_5pct": 500, "Td_50pct": 600, "char_yield": 55},
    "PEEK":   {"Td_5pct": 490, "Td_50pct": 570, "char_yield": 48},
    "PPS":    {"Td_5pct": 470, "Td_50pct": 560, "char_yield": 45},
    "PC":     {"Td_5pct": 420, "Td_50pct": 510, "char_yield": 25},
    "Nylon6": {"Td_5pct": 350, "Td_50pct": 430, "char_yield": 5},
    "PMMA":   {"Td_5pct": 270, "Td_50pct": 360, "char_yield": 0},
    "PS":     {"Td_5pct": 300, "Td_50pct": 400, "char_yield": 0},
    "PE":     {"Td_5pct": 390, "Td_50pct": 470, "char_yield": 0},
    "PLA":    {"Td_5pct": 290, "Td_50pct": 350, "char_yield": 2},
}

print(f"  {'폴리머':8s} {'Td5%(°C)':>10s} {'Td50%(°C)':>10s} {'잔탄율(%)':>10s} {'내열등급':>10s}")
print(f"  {'─'*8} {'─'*10} {'─'*10} {'─'*10} {'─'*10}")

for poly, data in sorted(thermal_stability.items(), key=lambda x: x[1]["Td_5pct"], reverse=True):
    grade = "🟢 우수" if data["Td_5pct"] >= 450 else ("🟡 양호" if data["Td_5pct"] >= 350 else "🔴 보통")
    bar = "█" * int(data["Td_5pct"] / 25)
    print(f"  {poly:8s} {data['Td_5pct']:>8d}° {data['Td_50pct']:>8d}° {data['char_yield']:>8d}% {grade:>10s} {bar}")

# ─────────────────────────────────────────────────────────────
# 9. 폴리머 응용 분류 및 선택 가이드
# ─────────────────────────────────────────────────────────────
print("\n\n🎯 Step 9: 폴리머 응용 분류")
print("-" * 50)

applications = {
    "범용 플라스틱": {
        "polymers": ["PE", "PP", "PS", "PVC"],
        "특징": "대량 생산, 저비용",
        "용도": "포장재, 건축, 가전",
    },
    "엔지니어링 플라스틱": {
        "polymers": ["PC", "Nylon6", "POM", "PBT"],
        "특징": "고강도, 내열성",
        "용도": "자동차, 전자, 기계 부품",
    },
    "슈퍼 엔지니어링 플라스틱": {
        "polymers": ["PEEK", "PPS", "PI", "LCP"],
        "특징": "극한 내열·내화학",
        "용도": "항공우주, 반도체, 의료",
    },
    "생분해성 폴리머": {
        "polymers": ["PLA", "PHA", "PCL", "PBS"],
        "특징": "환경 친화적",
        "용도": "포장재, 의료, 농업",
    },
    "엘라스토머": {
        "polymers": ["PDMS", "NR", "SBR", "EPDM"],
        "특징": "높은 탄성",
        "용도": "타이어, 씰, 의료기기",
    },
}

for category, info in applications.items():
    print(f"\n  [{category}]")
    print(f"    폴리머: {', '.join(info['polymers'])}")
    print(f"    특징: {info['특징']}")
    print(f"    용도: {info['용도']}")

# ─────────────────────────────────────────────────────────────
# 10. 종합 폴리머 특성 비교 레이더
# ─────────────────────────────────────────────────────────────
print("\n\n📋 Step 10: 종합 폴리머 특성 비교")
print("=" * 65)

def score_polymer(name, props, thermal):
    """폴리머 종합 점수"""
    # 정규화
    tg_score = min(max((props.get("Tg", 0) + 130) / 280, 0), 1)
    strength_score = min(props.get("tensile_MPa", 0) / 80, 1)
    thermal_score = min(thermal.get("Td_5pct", 300) / 500, 1) if thermal else 0.5
    cost_score = 0.8 if props.get("density", 1) < 1.2 else 0.5
    
    total = tg_score * 0.25 + strength_score * 0.30 + thermal_score * 0.30 + cost_score * 0.15
    return total

print(f"  {'폴리머':8s} {'Tg(°C)':>8s} {'강도(MPa)':>10s} {'Td5%(°C)':>10s} {'밀도':>6s} {'종합점수':>10s}")
print(f"  {'─'*8} {'─'*8} {'─'*10} {'─'*10} {'─'*6} {'─'*10}")

rankings = []
for poly, props in polymer_properties.items():
    thermal = thermal_stability.get(poly, {})
    score = score_polymer(poly, props, thermal)
    rankings.append((poly, score, props, thermal))

rankings.sort(key=lambda x: x[1], reverse=True)
for rank, (poly, score, props, thermal) in enumerate(rankings, 1):
    td = thermal.get("Td_5pct", "N/A")
    td_str = f"{td}°" if isinstance(td, int) else td
    grade = "🥇" if rank == 1 else ("🥈" if rank == 2 else ("🥉" if rank == 3 else "  "))
    bar = "█" * int(score * 25)
    print(f"  {poly:8s} {props['Tg']:>7d}° {props['tensile_MPa']:>8d}MPa {td_str:>10s} {props['density']:>5.2f} {grade} {score:.3f} {bar}")

print("\n" + "=" * 65)
print(" 🧪 폴리머/고분자 데모 완료!")
print("=" * 65)
