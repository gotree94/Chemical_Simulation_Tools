#!/usr/bin/env python3
"""
=============================================================
 신약/약물 설계 데모 (Drug Discovery Demo)
=============================================================
 RDKit 기반 약물 후보 분자 분석, ADMET 물성 예측,
 분자 유사도 검색, Lipinski Rule of Five 필터링
=============================================================
"""

import json
from rdkit import Chem
from rdkit.Chem import (
    Descriptors, AllChem, Draw, rdMolDescriptors,
    Crippen, Lipinski, rdFingerprintGenerator
)
from rdkit import DataStructs
import numpy as np

print("=" * 65)
print(" 🧬 신약/약물 설계 데모 — Drug Discovery Pipeline")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# 1. 약물 분자 라이브러리 정의
# ─────────────────────────────────────────────────────────────
print("\n📦 Step 1: 약물 분자 라이브러리 구축")
print("-" * 50)

drug_library = {
    "Aspirin":       {"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",       "class": "NSAID (소염진통제)"},
    "Ibuprofen":     {"smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "class": "NSAID (소염진통제)"},
    "Caffeine":      {"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  "class": "CNS 자극제"},
    "Paracetamol":   {"smiles": "CC(=O)NC1=CC=C(O)C=C1",          "class": "해열진통제"},
    "Metformin":     {"smiles": "CN(C)C(=N)NC(=N)N",              "class": "당뇨병 치료제"},
    "Penicillin_G":  {"smiles": "CC1(C)SC2C(NC(=O)CC3=CC=CC=C3)C(=O)N2C1C(=O)O", "class": "항생제"},
    "Diazepam":      {"smiles": "CN1C(=O)CN=C(C2=CC=CC=C2)C1=CC=C(Cl)C=1",       "class": "벤조디아제핀"},
    "Omeprazole":    {"smiles": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=CC=C3N2",    "class": "양성자펌프억제제"},
    "Sildenafil":    {"smiles": "CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C", "class": "PDE5 억제제"},
    "Atorvastatin":  {"smiles": "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "class": "스타틴"},
}

for name, info in drug_library.items():
    mol = Chem.MolFromSmiles(info["smiles"])
    if mol:
        info["mol"] = mol
        info["valid"] = True
        print(f"  ✅ {name:15s} | {info['class']}")
    else:
        info["mol"] = None
        info["valid"] = False
        print(f"  ❌ {name:15s} | 파싱 실패")

valid_drugs = {k: v for k, v in drug_library.items() if v.get("valid")}
print(f"\n  총 {len(valid_drugs)}/{len(drug_library)}개 분자 파싱 성공")

# ─────────────────────────────────────────────────────────────
# 2. 분자 물성 계산 (Molecular Descriptors)
# ─────────────────────────────────────────────────────────────
print("\n\n📊 Step 2: 분자 물성 계산 (Molecular Descriptors)")
print("-" * 50)

print(f"  {'약물':15s} {'분자량':>8s} {'LogP':>7s} {'HBD':>5s} {'HBA':>5s} {'TPSA':>7s} {'RotBond':>8s}")
print(f"  {'─'*15} {'─'*8} {'─'*7} {'─'*5} {'─'*5} {'─'*7} {'─'*8}")

descriptor_results = {}
for name, info in valid_drugs.items():
    mol = info["mol"]
    props = {
        "MW":       round(Descriptors.MolWt(mol), 2),
        "LogP":     round(Crippen.MolLogP(mol), 2),
        "HBD":      Descriptors.NumHDonors(mol),
        "HBA":      Descriptors.NumHAcceptors(mol),
        "TPSA":     round(Descriptors.TPSA(mol), 2),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
    }
    descriptor_results[name] = props
    print(f"  {name:15s} {props['MW']:>8.2f} {props['LogP']:>7.2f} {props['HBD']:>5d} {props['HBA']:>5d} {props['TPSA']:>7.2f} {props['RotBonds']:>8d}")

# ─────────────────────────────────────────────────────────────
# 3. Lipinski Rule of Five 필터링
# ─────────────────────────────────────────────────────────────
print("\n\n💊 Step 3: Lipinski Rule of Five (경구 약물 적합성)")
print("-" * 50)
print("  규칙: MW≤500, LogP≤5, HBD≤5, HBA≤10")
print()

lipinski_results = {}
for name, props in descriptor_results.items():
    violations = 0
    details = []
    if props["MW"] > 500:
        violations += 1
        details.append(f"MW={props['MW']}")
    if props["LogP"] > 5:
        violations += 1
        details.append(f"LogP={props['LogP']}")
    if props["HBD"] > 5:
        violations += 1
        details.append(f"HBD={props['HBD']}")
    if props["HBA"] > 10:
        violations += 1
        details.append(f"HBA={props['HBA']}")
    
    passed = violations <= 1
    lipinski_results[name] = {"passed": passed, "violations": violations}
    
    status = "✅ PASS" if passed else "❌ FAIL"
    violation_str = f"위반 {violations}개: {', '.join(details)}" if details else "위반 없음"
    print(f"  {name:15s} {status}  ({violation_str})")

pass_count = sum(1 for v in lipinski_results.values() if v["passed"])
print(f"\n  결과: {pass_count}/{len(lipinski_results)}개 약물이 Lipinski 규칙 통과")

# ─────────────────────────────────────────────────────────────
# 4. Veber 규칙 (경구 생체이용률 추가 필터)
# ─────────────────────────────────────────────────────────────
print("\n\n🔬 Step 4: Veber Rules (경구 생체이용률 추가 평가)")
print("-" * 50)
print("  규칙: TPSA≤140, RotBonds≤10")
print()

for name, props in descriptor_results.items():
    tpsa_ok = props["TPSA"] <= 140
    rot_ok = props["RotBonds"] <= 10
    passed = tpsa_ok and rot_ok
    
    status = "✅ PASS" if passed else "❌ FAIL"
    issues = []
    if not tpsa_ok:
        issues.append(f"TPSA={props['TPSA']}")
    if not rot_ok:
        issues.append(f"RotBonds={props['RotBonds']}")
    issue_str = f"위반: {', '.join(issues)}" if issues else "위반 없음"
    print(f"  {name:15s} {status}  ({issue_str})")

# ─────────────────────────────────────────────────────────────
# 5. 분자 유사도 분석 (Tanimoto Similarity)
# ─────────────────────────────────────────────────────────────
print("\n\n🔍 Step 5: 분자 유사도 분석 (Tanimoto Similarity)")
print("-" * 50)

# Morgan Fingerprint 생성 (radius=2, 2048 bits)
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fingerprints = {}
for name, info in valid_drugs.items():
    fingerprints[name] = fpgen.GetFingerprint(info["mol"])

# 아스피린과 다른 약물들의 유사도 계산
reference = "Aspirin"
print(f"  기준 분자: {reference}")
print(f"  {'약물':15s} {'Tanimoto':>10s} {'유사도':>8s}")
print(f"  {'─'*15} {'─'*10} {'─'*8}")

similarities = []
for name in valid_drugs:
    if name == reference:
        continue
    sim = DataStructs.TanimotoSimilarity(fingerprints[reference], fingerprints[name])
    similarities.append((name, sim))

similarities.sort(key=lambda x: x[1], reverse=True)
for name, sim in similarities:
    bar = "█" * int(sim * 20)
    print(f"  {name:15s} {sim:>10.4f}   {bar}")

# ─────────────────────────────────────────────────────────────
# 6. 분자 유사도 매트릭스
# ─────────────────────────────────────────────────────────────
print("\n\n📈 Step 6: 분자 유사도 매트릭스 (상위 5개)")
print("-" * 50)

top_drugs = list(valid_drugs.keys())[:5]
print(f"  {'':15s}", end="")
for name in top_drugs:
    print(f" {name[:8]:>8s}", end="")
print()

for name1 in top_drugs:
    print(f"  {name1:15s}", end="")
    for name2 in top_drugs:
        sim = DataStructs.TanimotoSimilarity(fingerprints[name1], fingerprints[name2])
        print(f" {sim:>8.3f}", end="")
    print()

# ─────────────────────────────────────────────────────────────
# 7. 약물 유사성 점수 (QED — Quantitative Estimate of Drug-likeness)
# ─────────────────────────────────────────────────────────────
print("\n\n⭐ Step 7: QED 약물 유사성 점수")
print("-" * 50)
print("  QED: 0(약물 비유사) ~ 1(약물 유사)")
print()

from rdkit.Chem import QED

qed_results = []
for name, info in valid_drugs.items():
    qed_score = QED.qed(info["mol"])
    qed_results.append((name, qed_score))

qed_results.sort(key=lambda x: x[1], reverse=True)
for rank, (name, score) in enumerate(qed_results, 1):
    bar = "█" * int(score * 30)
    cls = valid_drugs[name]["class"]
    print(f"  {rank:2d}. {name:15s} QED={score:.4f} {bar}  ({cls})")

# ─────────────────────────────────────────────────────────────
# 8. 합성 접근성 점수 (SA Score)
# ─────────────────────────────────────────────────────────────
print("\n\n🏭 Step 8: 합성 접근성 점수 (Synthetic Accessibility)")
print("-" * 50)
print("  SA Score: 1(쉬움) ~ 10(어려움)")
print()

from rdkit.Chem import RDConfig
import os, sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
try:
    import sascorer
    sa_available = True
except ImportError:
    sa_available = False

if sa_available:
    sa_results = []
    for name, info in valid_drugs.items():
        sa_score = sascorer.calculateScore(info["mol"])
        sa_results.append((name, sa_score))
    
    sa_results.sort(key=lambda x: x[1])
    for rank, (name, score) in enumerate(sa_results, 1):
        difficulty = "쉬움" if score < 3 else ("보통" if score < 5 else "어려움")
        bar = "█" * int(score * 3)
        print(f"  {rank:2d}. {name:15s} SA={score:.2f} {bar}  ({difficulty})")
else:
    print("  ⚠️ SA Score 모듈을 사용할 수 없습니다.")
    print("  대체: 분자 복잡도 기반 평가")
    for name, info in valid_drugs.items():
        mol = info["mol"]
        complexity = (Descriptors.NumRotatableBonds(mol) + 
                     Descriptors.RingCount(mol) * 2 + 
                     Descriptors.NumHeteroatoms(mol))
        difficulty = "쉬움" if complexity < 8 else ("보통" if complexity < 15 else "어려움")
        print(f"  {name:15s} 복잡도={complexity:3d}  ({difficulty})")

# ─────────────────────────────────────────────────────────────
# 9. 약물-약물 상호작용 기반 구조 클러스터링
# ─────────────────────────────────────────────────────────────
print("\n\n🗂️ Step 9: 구조 기반 약물 클러스터링")
print("-" * 50)

from rdkit.ML.Cluster import Butina

# 거리 매트릭스 계산
drug_names = list(valid_drugs.keys())
n = len(drug_names)
dist_matrix = []
for i in range(1, n):
    for j in range(i):
        sim = DataStructs.TanimotoSimilarity(fingerprints[drug_names[i]], fingerprints[drug_names[j]])
        dist_matrix.append(1 - sim)  # 거리 = 1 - 유사도

# Butina 클러스터링 (cutoff=0.6)
clusters = Butina.ClusterData(dist_matrix, n, distThresh=0.6, isDistData=True)
print(f"  클러스터 수: {len(clusters)} (cutoff=0.6)")
print()
for i, cluster in enumerate(clusters):
    members = [drug_names[idx] for idx in cluster]
    print(f"  Cluster {i+1}: {', '.join(members)}")

# ─────────────────────────────────────────────────────────────
# 10. 종합 약물 후보 평가 리포트
# ─────────────────────────────────────────────────────────────
print("\n\n📋 Step 10: 종합 약물 후보 평가 리포트")
print("=" * 65)

print(f"  {'약물':13s} {'MW':>6s} {'LogP':>6s} {'QED':>6s} {'Lipinski':>9s} {'종합':>6s}")
print(f"  {'─'*13} {'─'*6} {'─'*6} {'─'*6} {'─'*9} {'─'*6}")

for name, info in valid_drugs.items():
    props = descriptor_results[name]
    qed_score = QED.qed(info["mol"])
    lip = lipinski_results[name]
    
    # 종합 점수 (QED * Lipinski 패스 여부)
    overall = qed_score * (1.0 if lip["passed"] else 0.5)
    lip_str = "PASS" if lip["passed"] else f"FAIL({lip['violations']})"
    
    grade = "🟢" if overall > 0.5 else ("🟡" if overall > 0.3 else "🔴")
    print(f"  {name:13s} {props['MW']:>6.1f} {props['LogP']:>6.2f} {qed_score:>6.3f} {lip_str:>9s} {grade} {overall:.3f}")

print("\n" + "=" * 65)
print(" 🧬 신약/약물 설계 데모 완료!")
print("=" * 65)
