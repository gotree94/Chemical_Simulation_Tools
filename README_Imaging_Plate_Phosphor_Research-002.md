# 화학물질 합성 및 성능 시뮬레이션 도구 가이드

   * 📊 실제 DQE 계산식 포함한 수식 버전 추가
   * 📈 Trap energy band diagram 그림 생성
   * 🧪 Monte Carlo Python 코드 샘플 생성
   * 📘 석사 논문 목차 구조 설계
   * 🏭 산업 적용 관점 정리

## 📊 1️⃣ DQE 계산식 포함 수식 정리
* docs/DQE_and_MTF_theory.md

```
# DQE, MTF 및 영상 물리 수식 정리

## 1. MTF (Modulation Transfer Function)

MTF는 공간 주파수 f에서의 대비 전달 능력을 나타낸다.

MTF 정의:

MTF(f) = |FFT(LSF(x))|

여기서,

LSF(x) = d(ESF)/dx

- ESF: Edge Spread Function
- LSF: Line Spread Function
- FFT: Fast Fourier Transform

---

## 2. NPS (Noise Power Spectrum)

NPS는 영상 노이즈의 공간 주파수 분포를 나타낸다.

NPS(f) = (ΔxΔy / MN) · |FFT(I(x,y) − mean)|²

---

## 3. DQE (Detective Quantum Efficiency)

DQE는 시스템의 양자 검출 효율을 의미한다.

DQE(f) = SNR_out²(f) / SNR_in²(f)

또는

DQE(f) = (MTF(f))² / [ q₀ · NPS(f) ]

여기서,

- q₀ = 단위 면적당 입사 X-ray photon 수
- MTF(f) = 공간 주파수 응답
- NPS(f) = 노이즈 파워 스펙트럼

---

## 4. PSL Intensity 모델

PSL intensity I는 트랩 밀도와 비례한다.

I ∝ N_trap · σ · Φ

- N_trap: 트랩 밀도
- σ: 자극 단면적
- Φ: 레이저 플럭스

---

## 5. Trap Depth 계산

Thermoluminescence(TL) 피크 온도 T_m 기준:

E ≈ 2kT_m

- E: trap depth (eV)
- k: Boltzmann constant
- T_m: TL peak temperature (Kelvin)
```

📈 2️⃣ Trap Energy Band Diagram (텍스트 기반 도식)

docs/trap_energy_band_diagram.md

```
# Trap Energy Band Diagram

Energy ↑

  Conduction Band (CB)
  ----------------------
        ↑   ↑
        |   |  (Laser stimulation)
        |   |
     Trap Level (F-center)
     ----------------------
            ↓
     Eu²⁺ excited state
     ----------------------
            ↓  (Blue emission ~390-420 nm)
     Eu²⁺ ground state
     ----------------------
  Valence Band (VB)

설명:

1. X-ray 조사 → 전자 CB로 여기
2. 전자 Trap level에 포획
3. Laser 자극 → 전자 CB로 재여기
4. Eu²⁺ 중심에서 발광


원하시면 실제 논문용 벡터 다이어그램도 제작 가능합니다.
```

🧪 3️⃣ Monte Carlo 기반 X-ray 흡수 시뮬레이션 (Python 샘플)


simulations/scripts/monte_carlo_xray.py

```
import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
num_photons = 100000
thickness = 0.03  # cm
mu = 20  # attenuation coefficient (cm^-1)

absorbed = 0
transmitted = 0

for _ in range(num_photons):
    penetration_depth = -np.log(np.random.rand()) / mu
    
    if penetration_depth < thickness:
        absorbed += 1
    else:
        transmitted += 1

absorption_efficiency = absorbed / num_photons

print("Absorption Efficiency:", absorption_efficiency)

plt.bar(["Absorbed", "Transmitted"], [absorbed, transmitted])
plt.title("Monte Carlo X-ray Interaction Simulation")
plt.show()
```

확장 가능:
   * 에너지 스펙트럼 고려
   * 다층 구조 모델링
   * 산란 이벤트 포함
   * CsBr needle 구조 반영

📘 4️⃣ 석사 논문 목차 구조 설계

```
# 석사 학위 논문 목차 예시

## 제1장 서론
1.1 연구 배경  
1.2 CR 시스템 개요  
1.3 연구 필요성  
1.4 연구 목적  

---

## 제2장 이론적 배경
2.1 Photostimulated Luminescence  
2.2 Storage Phosphor 메커니즘  
2.3 BaFBr:Eu²⁺ 특성  
2.4 CsBr:Eu²⁺ 특성  
2.5 DQE 및 MTF 이론  

---

## 제3장 실험 방법
3.1 시료 합성 방법  
3.2 구조 분석 (XRD, SEM)  
3.3 광 특성 측정  
3.4 영상 성능 평가  

---

## 제4장 결과 및 분석
4.1 PSL 특성  
4.2 Trap depth 분석  
4.3 MTF 결과  
4.4 DQE 계산 결과  
4.5 BaFBr vs CsBr 비교  

---

## 제5장 결론
5.1 연구 요약  
5.2 학술적 기여  
5.3 산업적 의의  
5.4 향후 연구 방향  

---

참고문헌
```

🏭 5️⃣ 산업 적용 관점 분석

* docs/industrial_application.md

```
# 산업 적용 관점 분석

## 1. 의료 영상 시장

- 디지털 전환 가속
- 고해상도 요구 증가
- 저선량 촬영 기술 중요

CsBr 구조는 고해상도 시장에 적합

---

## 2. 산업용 NDT

- 용접 검사
- 항공 부품 검사
- 배관 검사

BaFBr은 비용 효율적 대안

---

## 3. 경쟁 기술 비교

| 기술 | 장점 | 단점 |
|------|------|------|
| CR | 재사용 가능 | 해상도 한계 |
| DR | 실시간 영상 | 고가 |
| Flat Panel | 고해상도 | 제조 비용 높음 |

---

## 4. 상용화 핵심 요소

- 제조 공정 안정성
- 원가 절감
- 내구성
- 반복 사용 수명

---

## 5. 미래 연구 방향

- 저선량 고감도 phosphor
- AI 기반 영상 보정
- Hybrid phosphor 구조
```

🔥 현재 프로젝트 수준
   * 이제 이 레포는:
      * 📚 학술 연구 가능
      * 🎓 석사 논문화 가능
      * 🏭 산업 기술 분석 가능
      * 💻 시뮬레이션 확장 가능
      * 📊 정량적 성능 분석 가능
    
