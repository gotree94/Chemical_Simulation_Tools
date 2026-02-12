# 화학물질 합성 및 성능 시뮬레이션 도구 가이드

## 개요

AI 기반 분자 설계와 시뮬레이션 도구를 활용하여 화학물질의 합성 경로 예측, 물성 시뮬레이션, 신규 물질 탐색을 수행하기 위한 환경 가이드입니다. 단백질 구조 예측(AlphaFold) 환경에서 확장하여 4개 분야에 걸친 도구들을 정리합니다.

### 대상 분야

| 분야 | 핵심 목표 | 대표 도구 |
|------|-----------|-----------|
| 🧬 신약/약물 설계 | 약물 후보 탐색, ADMET 예측, 도킹 | AutoDock, DeepChem, AiZynthFinder |
| 🔋 배터리/에너지 소재 | 전해질 설계, 전극 물질 탐색 | VASP, Pymatgen, GNoME |
| 🧪 폴리머/고분자 | 물성 예측, 중합 시뮬레이션 | LAMMPS, PolymerGenome, RDKit |
| 💎 반도체/전자소재 | 밴드갭 예측, 결정 구조 탐색 | Quantum ESPRESSO, AFLOW, MEGNet |

### 시스템 요구사항

| 항목 | 최소 | 권장 (현재 환경) |
|------|------|------------------|
| GPU | 8GB VRAM | RTX 5090 24GB ✅ |
| RAM | 16GB | 64GB+ |
| Storage | 100GB | 1TB+ (데이터베이스 포함) |
| CUDA | 12.1+ | 12.4 ✅ |
| OS | Ubuntu 22.04 | Ubuntu 22.04 ✅ |

---

## 공통 기반 도구 (모든 분야 공통)

모든 분야에서 기초로 사용되는 핵심 도구들입니다. 먼저 이 도구들을 설치한 후 분야별 도구를 추가합니다.

### RDKit — 분자 정보학 플랫폼

분자 구조 생성, 조작, 시각화, 물성 계산의 기본 도구입니다.

```bash
conda create -n cheminformatics python=3.11 -y
conda activate cheminformatics
conda install -c conda-forge rdkit -y
```

**주요 기능:**
- SMILES/SMARTS 기반 분자 표현 및 변환
- 분자 지문(Fingerprint) 계산 및 유사도 검색
- 2D/3D 좌표 생성 및 시각화
- 분자 물성 계산 (분자량, LogP, TPSA 등)
- 서브구조 검색 및 반응 SMARTS

**기본 사용 예제:**
```python
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem

# SMILES로 분자 생성 (아스피린)
mol = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")

# 물성 계산
print(f"분자량: {Descriptors.MolWt(mol):.2f}")
print(f"LogP: {Descriptors.MolLogP(mol):.2f}")
print(f"수소결합 공여체: {Descriptors.NumHDonors(mol)}")
print(f"수소결합 수용체: {Descriptors.NumHAcceptors(mol)}")

# 3D 좌표 생성
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol)
```

### DeepChem — AI 기반 화학 플랫폼

머신러닝/딥러닝 기반 분자 물성 예측, 독성 예측, 생성 모델을 위한 통합 플랫폼입니다.

```bash
pip install deepchem
pip install torch torchvision  # PyTorch 백엔드
```

**주요 기능:**
- 분자 물성 예측 (GCN, MPNN, AttentiveFP)
- 독성 예측 (Tox21, ClinTox 데이터셋)
- 용해도, 약물 유사성 예측
- 분자 생성 모델 (VAE, GAN)
- 분자 지문 기반 분류/회귀

### ASE (Atomic Simulation Environment) — 원자 시뮬레이션 인터페이스

다양한 시뮬레이터와 연동되는 Python 기반 원자 시뮬레이션 프레임워크입니다.

```bash
pip install ase
```

**주요 기능:**
- 원자/분자/결정 구조 생성 및 조작
- 다양한 계산기(Calculator) 인터페이스 (VASP, GPAW, Gaussian 등)
- 구조 최적화, 분자 동역학(MD)
- 열역학적 물성 계산

### Open Babel — 분자 파일 형식 변환기

```bash
conda install -c conda-forge openbabel -y
```

**주요 기능:**
- 110+ 화학 파일 형식 간 변환 (SMILES, SDF, MOL2, PDB, CIF 등)
- 3D 좌표 생성, 에너지 최소화
- 명령줄 기반 배치 처리

---

## 1. 🧬 신약/약물 설계 (Drug Discovery)

약물 후보 물질의 탐색, 합성 경로 예측, 약물-단백질 상호작용 분석을 위한 도구입니다.

### 1.1 합성 경로 예측

#### AiZynthFinder — AI 역합성 분석

AstraZeneca에서 개발한 AI 기반 역합성 경로 예측 도구입니다. 목표 분자를 입력하면 합성 가능한 경로를 자동으로 탐색합니다.

```bash
conda create -n aizynthfinder python=3.10 -y
conda activate aizynthfinder
pip install aizynthfinder[all]

# 데이터 파일 다운로드
download_public_data /path/to/aizynthfinder/data
```

**주요 기능:**
- Monte Carlo Tree Search 기반 역합성 경로 탐색
- 반응 템플릿 기반 변환 예측
- 합성 가능성 점수 (SCScore) 제공
- 다단계 합성 경로 시각화

**사용 예제:**
```python
from aizynthfinder.aizynthfinder import AiZynthFinder

finder = AiZynthFinder(configfile="config.yml")
finder.target_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # 아스피린
finder.tree_search()
finder.build_routes()

stats = finder.extract_statistics()
print(f"찾은 경로 수: {stats['number_of_routes']}")
```

#### IBM RXN for Chemistry — 화학 반응 예측 API

Transformer 기반 화학 반응 예측 서비스입니다 (무료 API 제공).

```bash
pip install rxn4chemistry
```

**주요 기능:**
- 정방향 반응 예측 (반응물 → 생성물)
- 역합성 분석 (생성물 → 반응물)
- 반응 조건 예측
- 무료 웹 API 제공 (https://rxn.res.ibm.com)

### 1.2 분자 도킹 (약물-단백질 상호작용)

#### AutoDock Vina — 분자 도킹

약물 후보 분자가 단백질 활성 부위에 어떻게 결합하는지 예측합니다.

```bash
conda install -c conda-forge autodock-vina -y
pip install meeko  # PDBQT 파일 준비 도구
pip install vina   # Python 바인딩
```

**주요 기능:**
- 약물-단백질 결합 친화도 예측
- 결합 포즈 탐색 및 시각화
- 가상 스크리닝 (대량 후보 물질 평가)

**사용 예제:**
```python
from vina import Vina

v = Vina(sf_name='vina')
v.set_receptor("protein.pdbqt")
v.set_ligand_from_file("ligand.pdbqt")
v.compute_vina_maps(center=[0, 0, 0], box_size=[20, 20, 20])
v.dock(exhaustiveness=32, n_poses=10)
v.write_poses("output.pdbqt", n_poses=5)
```

#### DiffDock — AI 기반 분자 도킹

MIT에서 개발한 확산 모델(Diffusion Model) 기반 도킹 도구입니다.

```bash
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
conda env create -f environment.yml
conda activate diffdock
```

**주요 기능:**
- 결합 부위를 사전에 지정할 필요 없음 (Blind docking)
- 확산 모델 기반 포즈 생성
- AutoDock Vina보다 높은 성공률 보고

### 1.3 ADMET 예측

#### ADMETlab / pkCSM — 약물 동태 예측

약물의 흡수(Absorption), 분포(Distribution), 대사(Metabolism), 배설(Excretion), 독성(Toxicity) 예측 도구입니다.

```bash
# DeepChem 기반 ADMET 예측
pip install deepchem

# 또는 ADMETlab 웹 서비스 활용
# https://admetmesh.scbdd.com/
```

**주요 예측 항목:**
- 경구 생체이용률 (Oral Bioavailability)
- 혈뇌장벽 투과성 (BBB Permeability)
- CYP450 억제/대사 예측
- hERG 독성 (심장 독성)
- Lipinski Rule of Five 평가

### 1.4 분자 생성 모델

#### REINVENT — 강화학습 기반 분자 생성

AstraZeneca에서 개발한 강화학습 기반 신규 약물 분자 생성 도구입니다.

```bash
git clone https://github.com/MolecularAI/REINVENT4.git
cd REINVENT4
conda env create -f environment.yml
conda activate reinvent4
```

**주요 기능:**
- RNN/Transformer 기반 SMILES 생성
- 강화학습으로 원하는 물성 최적화
- 다목적 최적화 (활성 + 독성 + 합성 용이성)
- 스캐폴드 기반 분자 설계

#### MolGAN — GAN 기반 분자 생성

```bash
pip install deepchem  # MolGAN 구현 포함
```

### 1.5 신약 설계 워크플로우 요약

```
[타겟 단백질 선택]
     ↓
[AlphaFold → 단백질 구조 예측]
     ↓
[REINVENT/MolGAN → 약물 후보 분자 생성]
     ↓
[DeepChem → ADMET 물성 예측 & 필터링]
     ↓
[AutoDock Vina/DiffDock → 분자 도킹 시뮬레이션]
     ↓
[AiZynthFinder → 합성 경로 예측]
     ↓
[GROMACS/OpenMM → 분자 동역학 시뮬레이션 (결합 안정성 확인)]
```

---

## 2. 🔋 배터리/에너지 소재 (Battery & Energy Materials)

배터리 전해질, 전극 물질, 촉매 소재의 탐색과 성능 시뮬레이션을 위한 도구입니다.

### 2.1 제일원리 계산 (First-Principles Calculation)

#### VASP (Vienna Ab initio Simulation Package) — DFT 계산

고체 물질의 전자 구조를 밀도범함수이론(DFT)으로 계산하는 표준 도구입니다.

```bash
# 상용 소프트웨어 — 라이선스 필요
# 학술 라이선스: https://www.vasp.at/
# 설치 후 환경 변수 설정
export PATH=/opt/vasp/bin:$PATH
```

**주요 기능:**
- 전자 구조 계산 (밴드 구조, 상태 밀도)
- 구조 최적화 및 격자 상수 예측
- 이온 전도도, 확산 경로 계산
- 전압 프로파일 예측 (배터리 전극)

> ⚠️ **참고:** VASP는 상용 소프트웨어입니다. 무료 대안으로 Quantum ESPRESSO를 사용할 수 있습니다.

#### Quantum ESPRESSO — 오픈소스 DFT 계산

```bash
# Ubuntu 패키지 설치
sudo apt install quantum-espresso -y

# 또는 소스 빌드 (GPU 가속)
git clone https://github.com/QEF/q-e.git
cd q-e
./configure --with-cuda=yes
make all
```

**주요 기능:**
- 평면파 기반 DFT 계산
- 구조 최적화, 포논 계산
- 밴드 구조, 상태 밀도 (DOS)
- GPU 가속 지원

### 2.2 소재 데이터베이스 & 탐색

#### Pymatgen — Materials Project Python API

Materials Project 데이터베이스와 연동되는 소재 분석 Python 라이브러리입니다.

```bash
pip install pymatgen mp-api
```

**주요 기능:**
- Materials Project DB (15만+ 화합물) 접근
- 결정 구조 생성, 변환, 분석
- 상평형도(Phase Diagram) 계산
- 전기화학 분석 (배터리 전압, 용량)
- 밴드갭, 탄성 계수 등 물성 조회

**사용 예제:**
```python
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

# Materials Project API 키 필요 (무료 등록)
with MPRester("YOUR_API_KEY") as mpr:
    # LiFePO4 (리튬인산철 배터리 양극재) 검색
    docs = mpr.materials.summary.search(
        formula="LiFePO4",
        fields=["material_id", "formula_pretty", "band_gap", "energy_above_hull"]
    )
    for doc in docs:
        print(f"{doc.formula_pretty}: 밴드갭={doc.band_gap:.2f} eV")
```

#### GNoME (Google DeepMind) — AI 신소재 탐색

Graph Neural Network 기반 신소재 결정 구조 탐색 도구입니다. 2024년 Nature에 발표되어 220만개의 안정한 결정 구조를 예측했습니다.

```bash
git clone https://github.com/google-deepmind/materials_discovery.git
cd materials_discovery
pip install -e .
```

**주요 기능:**
- GNN 기반 결정 구조 안정성 예측
- 신규 결정 구조 생성
- 에너지/안정성 기반 후보 필터링
- Materials Project DB와 통합

### 2.3 분자 동역학 (배터리 전해질 시뮬레이션)

#### LAMMPS — 대규모 분자 동역학

```bash
# Ubuntu 패키지 설치
sudo apt install lammps -y

# 또는 GPU 가속 빌드
git clone https://github.com/lammps/lammps.git
cd lammps
mkdir build && cd build
cmake ../cmake -D PKG_GPU=on -D GPU_API=cuda
make -j$(nproc)
```

**주요 기능:**
- 대규모 원자/분자 시스템 MD 시뮬레이션
- 전해질 이온 전도도 계산
- 계면 특성 분석 (전극-전해질)
- ReaxFF 반응 시뮬레이션
- GPU(CUDA) 가속 지원

#### GROMACS — GPU 가속 MD 시뮬레이션

```bash
sudo apt install gromacs -y

# GPU 가속 빌드
sudo apt install gromacs-gpu -y
```

**주요 기능:**
- 고속 GPU 가속 MD 시뮬레이션
- 전해질 용매화 구조 분석
- 자유 에너지 계산
- 이온 확산 계수 계산

### 2.4 배터리 시뮬레이션 전용 도구

#### PyBaMM — 배터리 모델링

Python 기반 배터리 수학 모델링 프레임워크입니다.

```bash
pip install pybamm
```

**주요 기능:**
- SPM (Single Particle Model), DFN 모델
- 충방전 프로파일 시뮬레이션
- 열 거동 분석
- 열화 모델링

**사용 예제:**
```python
import pybamm

model = pybamm.lithium_ion.DFN()
sim = pybamm.Simulation(model)
sim.solve([0, 3600])  # 1시간 시뮬레이션
sim.plot()
```

### 2.5 배터리 소재 탐색 워크플로우 요약

```
[Materials Project DB → 후보 소재 탐색 (Pymatgen)]
     ↓
[GNoME → 신규 결정 구조 생성/안정성 예측]
     ↓
[Quantum ESPRESSO/VASP → DFT 전자 구조 계산]
     ↓
[Pymatgen → 전압 프로파일, 상평형도 분석]
     ↓
[LAMMPS/GROMACS → 이온 전도도, 확산 MD 시뮬레이션]
     ↓
[PyBaMM → 셀 레벨 성능 시뮬레이션]
```

---

## 3. 🧪 폴리머/고분자 (Polymer Science)

고분자 물질의 물성 예측, 중합 시뮬레이션, 신규 폴리머 설계를 위한 도구입니다.

### 3.1 폴리머 물성 예측

#### Polymer Genome — AI 폴리머 물성 예측

머신러닝 기반 폴리머 물성 예측 플랫폼입니다 (Georgia Tech 개발).

```bash
# 웹 서비스 제공
# https://www.polymergenome.org/

# Python API
pip install polygnn  # Polymer Graph Neural Network
```

**주요 예측 물성:**
- 유리전이온도 (Tg)
- 밴드갭
- 유전 상수
- 열전도율
- 용해도 파라미터

#### polyBERT — 트랜스포머 기반 폴리머 표현

```bash
pip install polyBERT
```

**주요 기능:**
- PSMILES (Polymer SMILES) 기반 폴리머 표현
- 사전학습된 트랜스포머로 폴리머 임베딩 생성
- 다운스트림 물성 예측 파인튜닝

### 3.2 분자 동역학 시뮬레이션

#### LAMMPS — 폴리머 MD 시뮬레이션

```bash
# 폴리머 전용 패키지 포함 빌드
cd lammps/build
cmake ../cmake -D PKG_MOLECULE=on -D PKG_KSPACE=on -D PKG_RIGID=on
make -j$(nproc)
```

**폴리머 관련 주요 기능:**
- 폴리머 체인 생성 및 평형화
- 기계적 물성 계산 (인장 강도, 탄성 계수)
- Tg 예측 (냉각 시뮬레이션)
- 결정화 시뮬레이션
- 폴리머 블렌드 상분리

#### Moltemplate — LAMMPS 입력 파일 생성기

```bash
pip install moltemplate
```

**주요 기능:**
- 복잡한 폴리머 체인 토폴로지 자동 생성
- 반복 단위 기반 폴리머 구조 구축
- LAMMPS 데이터 파일 자동 변환
- 다양한 역장(Force Field) 지원

### 3.3 양자화학 수준 계산

#### ORCA — 무료 양자화학 패키지

```bash
# 학술용 무료 다운로드
# https://orcaforum.kofo.mpg.de/
# 다운로드 후:
tar -xf orca_5_0_4_linux_x86-64_shared_openmpi411.tar.xz
export PATH=/opt/orca:$PATH
```

**주요 기능:**
- DFT, HF, MP2, CCSD(T) 계산
- 단량체/소량체 반응 에너지 계산
- UV-Vis 스펙트럼 예측
- NMR 화학적 이동 계산
- 병렬 계산 지원

#### Gaussian — 표준 양자화학 패키지

```bash
# 상용 소프트웨어 — 라이선스 필요
# https://gaussian.com/
```

### 3.4 폴리머 정보학

#### RDKit + BigSMILES — 폴리머 표현

```python
from rdkit import Chem

# 단량체 구조 (스타이렌)
styrene = Chem.MolFromSmiles("C=Cc1ccccc1")

# 폴리머 반복 단위 표현
# BigSMILES 표기법: {[]CC(c1ccccc1)[]}
# 폴리스타이렌의 반복 단위
repeat_unit = Chem.MolFromSmiles("CC(c1ccccc1)")
```

### 3.5 폴리머 설계 워크플로우 요약

```
[RDKit → 단량체/반복 단위 설계]
     ↓
[Polymer Genome / polyBERT → 물성 예측 (Tg, 밴드갭 등)]
     ↓
[ORCA → 단량체 반응성, 에너지 계산]
     ↓
[Moltemplate → 폴리머 체인 구조 생성]
     ↓
[LAMMPS → MD 시뮬레이션 (기계적 물성, Tg)]
     ↓
[DeepChem → ML 기반 스크리닝 & 최적화]
```

---

## 4. 💎 반도체/전자소재 (Semiconductor & Electronic Materials)

반도체 물질의 밴드갭 예측, 결정 구조 탐색, 전자 수송 특성 분석을 위한 도구입니다.

### 4.1 제일원리 전자 구조 계산

#### Quantum ESPRESSO — DFT 전자 구조

```bash
# 반도체 관련 주요 계산
# pw.x: Self-consistent field (SCF) 계산
# bands.x: 밴드 구조 계산
# dos.x: 상태 밀도 계산
# ph.x: 포논 계산
# pp.x: 후처리 (전하 밀도, 파동 함수 시각화)
```

**반도체 분석 주요 기능:**
- 밴드 구조 및 밴드갭 계산
- 유효 질량 추출
- 상태 밀도 (DOS, PDOS)
- 유전 함수, 광학 스펙트럼
- 포논 분산 관계

#### GPAW — GPU 가속 DFT

```bash
pip install gpaw
gpaw install-data  # 의사 포텐셜 데이터
```

**주요 기능:**
- 실공간/평면파/LCAO 기저 함수 지원
- GW 근사 (정확한 밴드갭)
- 시간의존 DFT (TD-DFT, 광학 특성)
- ASE와 완벽 통합

### 4.2 소재 데이터베이스 & ML 탐색

#### AFLOW — 자동화 소재 탐색 프레임워크

Duke University에서 개발한 고처리량(high-throughput) 소재 계산 프레임워크입니다.

```bash
pip install aflow
```

**주요 기능:**
- 360만+ 화합물 데이터베이스
- REST API 기반 데이터 접근
- 밴드갭, 탄성 계수, 열역학 물성 조회
- 자동화된 DFT 계산 워크플로우
- 머신러닝 물성 예측 (AFLOW-ML)

**사용 예제:**
```python
import aflow
from aflow import K  # Keyword 모듈

# 밴드갭 1~2 eV인 반도체 검색
results = aflow.search().filter(
    K.Egap > 1.0,
    K.Egap < 2.0,
    K.nspecies == 2  # 이원계 화합물
).select(
    K.compound, K.Egap, K.crystal_system
)

for material in results:
    print(f"{material.compound}: 밴드갭={material.Egap:.2f} eV")
```

#### MEGNet — Graph Neural Network 소재 물성 예측

```bash
pip install megnet
```

**주요 기능:**
- 결정 구조 → 그래프 변환 → 물성 예측
- 밴드갭, 형성 에너지, 탄성 계수 예측
- Materials Project 데이터로 사전학습
- 전이학습 가능

#### ALIGNN — Atomistic Line Graph Neural Network

NIST에서 개발한 최신 GNN 기반 소재 물성 예측 도구입니다.

```bash
pip install alignn
```

**주요 기능:**
- 결정 구조 기반 물성 예측 (SOTA 성능)
- 밴드갭, 형성 에너지, 열전도율 등
- JARVIS 데이터베이스 연동
- 50+ 물성 사전학습 모델 제공

### 4.3 전자 수송 특성

#### BoltzTraP2 — 볼츠만 수송 이론

DFT 밴드 구조 데이터에서 전자 수송 특성을 계산합니다.

```bash
pip install BoltzTraP2
```

**주요 기능:**
- 전기 전도도, 제벡 계수 계산
- 열전도율 (전자 기여분)
- 열전 성능지수 (ZT) 평가
- 온도 의존성 분석

#### Wannier90 — 맥시멀 국소화 와니어 함수

```bash
git clone https://github.com/wannier-developers/wannier90.git
cd wannier90
make
```

**주요 기능:**
- DFT 결과에서 타이트바인딩 해밀토니안 추출
- 정밀한 밴드 보간
- 베리 위상, 토폴로지 불변량 계산
- 표면 상태 계산

### 4.4 반도체 공정 시뮬레이션

#### DEVSIM — 반도체 소자 시뮬레이션

```bash
pip install devsim
```

**주요 기능:**
- TCAD 수준 소자 시뮬레이션
- 드리프트-확산 방정식 풀기
- PN 접합, MOSFET 시뮬레이션
- I-V 특성, C-V 특성

### 4.5 반도체 소재 탐색 워크플로우 요약

```
[AFLOW/Materials Project → 후보 소재 데이터베이스 탐색]
     ↓
[MEGNet/ALIGNN → ML 기반 물성 스크리닝 (밴드갭, 안정성)]
     ↓
[Quantum ESPRESSO/GPAW → DFT 전자 구조 정밀 계산]
     ↓
[Wannier90 → 타이트바인딩 모델 추출]
     ↓
[BoltzTraP2 → 전자 수송 특성 계산]
     ↓
[DEVSIM → 소자 레벨 시뮬레이션]
```

---

## 분야별 도구 비교 매트릭스

### 오픈소스 여부 & GPU 가속 지원

| 도구 | 분야 | 오픈소스 | GPU 가속 | Python API | 난이도 |
|------|------|----------|----------|------------|--------|
| RDKit | 공통 | ✅ | ❌ | ✅ | ⭐⭐ |
| DeepChem | 공통 | ✅ | ✅ | ✅ | ⭐⭐ |
| ASE | 공통 | ✅ | ❌ | ✅ | ⭐⭐ |
| AiZynthFinder | 신약 | ✅ | ✅ | ✅ | ⭐⭐⭐ |
| AutoDock Vina | 신약 | ✅ | ❌ | ✅ | ⭐⭐ |
| DiffDock | 신약 | ✅ | ✅ | ✅ | ⭐⭐⭐ |
| REINVENT | 신약 | ✅ | ✅ | ✅ | ⭐⭐⭐⭐ |
| VASP | 배터리 | ❌ (상용) | ✅ | ✅ | ⭐⭐⭐⭐ |
| Quantum ESPRESSO | 배터리/반도체 | ✅ | ✅ | ✅ | ⭐⭐⭐⭐ |
| Pymatgen | 배터리 | ✅ | ❌ | ✅ | ⭐⭐ |
| GNoME | 배터리 | ✅ | ✅ | ✅ | ⭐⭐⭐ |
| LAMMPS | 배터리/폴리머 | ✅ | ✅ | ✅ | ⭐⭐⭐ |
| GROMACS | 신약/배터리 | ✅ | ✅ | ✅ | ⭐⭐⭐ |
| PyBaMM | 배터리 | ✅ | ❌ | ✅ | ⭐⭐ |
| Polymer Genome | 폴리머 | ✅ (웹) | ❌ | ✅ | ⭐ |
| ORCA | 폴리머 | ✅ (학술) | ❌ | ❌ | ⭐⭐⭐⭐ |
| AFLOW | 반도체 | ✅ | ❌ | ✅ | ⭐⭐ |
| MEGNet | 반도체 | ✅ | ✅ | ✅ | ⭐⭐⭐ |
| ALIGNN | 반도체 | ✅ | ✅ | ✅ | ⭐⭐⭐ |
| BoltzTraP2 | 반도체 | ✅ | ❌ | ✅ | ⭐⭐⭐ |
| DEVSIM | 반도체 | ✅ | ❌ | ✅ | ⭐⭐⭐ |

### 추천 시작 순서

전 분야를 아우르는 학습 로드맵입니다:

```
Phase 1 — 공통 기반 (1~2주)
├── RDKit 설치 및 분자 표현 기초
├── DeepChem 설치 및 ML 물성 예측 기초
└── ASE 설치 및 원자 시뮬레이션 기초

Phase 2 — 분야별 탐색 (분야당 2~3주)
├── 신약: AiZynthFinder → AutoDock Vina → REINVENT
├── 배터리: Pymatgen → Quantum ESPRESSO → PyBaMM
├── 폴리머: Polymer Genome → LAMMPS → ORCA
└── 반도체: AFLOW → ALIGNN → Quantum ESPRESSO

Phase 3 — 고급 시뮬레이션 (분야당 4주+)
├── 신약: DiffDock + GROMACS MD
├── 배터리: VASP DFT + LAMMPS MD
├── 폴리머: LAMMPS 대규모 MD + ORCA 양자화학
└── 반도체: Wannier90 + BoltzTraP2 수송 계산
```

---

## 참고 자료

### 데이터베이스

| 데이터베이스 | URL | 설명 |
|-------------|-----|------|
| Materials Project | https://materialsproject.org | 15만+ 무기 화합물 |
| AFLOW | https://aflow.org | 360만+ 화합물 |
| JARVIS | https://jarvis.nist.gov | NIST 소재 데이터 |
| ChEMBL | https://www.ebi.ac.uk/chembl/ | 생물활성 약물 데이터 |
| ZINC | https://zinc.docking.org | 가상 스크리닝용 분자 DB |
| PubChem | https://pubchem.ncbi.nlm.nih.gov | 1.1억+ 화합물 |
| Polymer Genome | https://www.polymergenome.org | 폴리머 물성 DB |
| COD | https://www.crystallography.net | 결정 구조 DB |

### 온라인 학습 리소스

| 리소스 | URL | 설명 |
|--------|-----|------|
| DeepChem Tutorials | https://deepchem.io/tutorials | AI 화학 튜토리얼 |
| Materials Project Workshop | https://workshop.materialsproject.org | 소재 정보학 워크숍 |
| Quantum ESPRESSO School | https://www.quantum-espresso.org/resources/tutorials | DFT 계산 튜토리얼 |
| LAMMPS Tutorials | https://docs.lammps.org/Tutorials.html | MD 시뮬레이션 튜토리얼 |

---

## 환경 구성 요약

현재 시스템(RTX 5090 24GB, Ubuntu)에서 권장하는 conda 환경 구성:

```bash
# 1. 공통 기반 환경
conda create -n cheminformatics python=3.11 -y
conda activate cheminformatics
conda install -c conda-forge rdkit openbabel -y
pip install deepchem ase

# 2. 신약 설계 환경
conda create -n drugdesign python=3.10 -y
conda activate drugdesign
pip install aizynthfinder[all] vina meeko deepchem

# 3. 소재 과학 환경
conda create -n materials python=3.11 -y
conda activate materials
pip install pymatgen mp-api megnet alignn pybamm

# 4. 분자 동역학 환경
conda create -n molsim python=3.11 -y
conda activate molsim
pip install gromacs lammps-python ase
sudo apt install gromacs-gpu lammps -y

# 5. 양자화학 환경
conda create -n quantumchem python=3.11 -y
conda activate quantumchem
pip install ase gpaw BoltzTraP2
sudo apt install quantum-espresso -y
```

---

## 라이선스 정보

| 도구 | 라이선스 | 비용 |
|------|----------|------|
| RDKit | BSD | 무료 |
| DeepChem | MIT | 무료 |
| LAMMPS | GPL | 무료 |
| GROMACS | LGPL | 무료 |
| Quantum ESPRESSO | GPL | 무료 |
| ORCA | 학술 무료 | 학술 무료 / 상용 유료 |
| VASP | 상용 | 유료 (학술 할인) |
| Gaussian | 상용 | 유료 |

---

> **Note:** 이 문서는 AlphaFold/ColabFold 단백질 구조 예측 환경에서 화학물질 합성 및 시뮬레이션으로 확장하기 위한 가이드입니다. 각 도구의 상세 설치 및 사용법은 별도 문서로 정리할 예정입니다.
