# 화학물질 합성 및 성능 시뮬레이션 도구 가이드

   * 실제 Trap depth fitting 코드 작성
   * DQE 자동 계산 Python 코드
   * 논문 제출용 그래프 스타일 템플릿
   * SCI 논문 초록 샘플 작성
   * 연구 제안서(연구비 신청용) 작성


🧪 1️⃣ Trap Depth Fitting Python 코드 (Thermoluminescence 기반)

* simulations/scripts/trap_depth_fitting.py

```
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Boltzmann constant (eV/K)
k = 8.617e-5  

# TL glow peak model (first-order kinetics approximation)
def glow_curve(T, E, s, N0):
    beta = 1  # heating rate (K/s) - adjustable
    return N0 * np.exp(-E/(k*T)) * np.exp(-(s/beta) * np.exp(-E/(k*T)))

# Example temperature data (Kelvin)
T = np.linspace(300, 600, 200)

# Simulated TL signal (replace with real data)
true_E = 0.9  # eV
true_s = 1e12
true_N0 = 1e5
TL_signal = glow_curve(T, true_E, true_s, true_N0)
TL_signal += np.random.normal(0, 1000, size=len(T))

# Curve fitting
popt, _ = curve_fit(glow_curve, T, TL_signal, p0=[1.0, 1e12, 1e5])

E_fit, s_fit, N0_fit = popt

print("Fitted Trap Depth (E):", E_fit, "eV")
print("Frequency Factor (s):", s_fit)

# Plot
plt.figure()
plt.scatter(T, TL_signal, s=10, label="Experimental Data")
plt.plot(T, glow_curve(T, *popt), 'r', label="Fitted Curve")
plt.xlabel("Temperature (K)")
plt.ylabel("TL Intensity")
plt.legend()
plt.title("Trap Depth Fitting")
plt.show()
```

✔ 실제 TL 데이터로 교체하면 trap depth 자동 계산 가능

📊 2️⃣ DQE 자동 계산 Python 코드

* simulations/scripts/dqe_calculation.py

```
import numpy as np
import matplotlib.pyplot as plt

# Spatial frequency axis
f = np.linspace(0.1, 5, 100)

# Example MTF (replace with measured data)
MTF = np.exp(-0.5 * f)

# Example NPS (replace with measured data)
NPS = 0.01 + 0.02 * f

# Incident photon fluence
q0 = 10000  

DQE = (MTF**2) / (q0 * NPS)

plt.figure()
plt.plot(f, DQE)
plt.xlabel("Spatial Frequency (lp/mm)")
plt.ylabel("DQE")
plt.title("DQE vs Spatial Frequency")
plt.grid(True)
plt.show()

print("DQE at 0 frequency:", DQE[0])
```

✔ 실제 측정 MTF, NPS 데이터 입력 시 자동 계산

📈 3️⃣ 논문 제출용 그래프 스타일 템플릿 (Matplotlib SCI 스타일)

* simulations/scripts/scientific_plot_style.py
```
import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-whitegrid')

plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 12,
    "axes.linewidth": 1.2,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "lines.linewidth": 2,
})

# Example plot
x = [1,2,3,4,5]
y = [2,4,6,8,10]

plt.plot(x, y)
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Publication Quality Figure")
plt.tight_layout()
plt.savefig("publication_figure.tif", format="tif")
plt.show()
```

✔ SCI 저널 제출용 600 dpi TIFF 자동 저장

📘 4️⃣ SCI 논문 초록 샘플 (CR 형광체 연구용)

```
Abstract

Photostimulable phosphors play a crucial role in computed radiography (CR) systems by enabling efficient X-ray energy storage and subsequent luminescence readout. In this study, we investigated the trap characteristics and imaging performance of Eu²⁺-activated BaFBr and CsBr phosphors. Thermoluminescence analysis was performed to determine trap depth and frequency factors, while modulation transfer function (MTF) and detective quantum efficiency (DQE) were evaluated to assess imaging performance.

The results indicate that CsBr:Eu²⁺ exhibits superior spatial resolution due to its columnar crystal structure, leading to reduced light scattering and improved DQE at higher spatial frequencies. In contrast, BaFBr:Eu²⁺ demonstrated stable trap characteristics and cost-effective fabrication advantages. The correlation between trap density and PSL intensity was quantitatively analyzed, providing insights into phosphor optimization strategies.

These findings contribute to the development of high-resolution, low-dose CR imaging systems and offer practical guidelines for next-generation storage phosphor design.
```

🏆 5️⃣ 연구 제안서 (연구비 신청용 요약본)

```
연구 과제명:
고해상도 저선량 Computed Radiography용 차세대 저장 형광체 개발

1. 연구 배경
디지털 의료 영상 기술의 발전과 함께 저선량 고해상도 영상 기술의 수요가 증가하고 있다. 기존 BaFBr 기반 CR 시스템은 안정성이 높으나 해상도 및 DQE 측면에서 한계가 존재한다.

2. 연구 목표
- Eu²⁺ 활성 저장 형광체의 트랩 구조 정량 분석
- Trap depth 및 PSL 효율 최적화
- CsBr needle 구조의 광 산란 저감 효과 분석
- DQE 20% 이상 향상 목표

3. 연구 내용
(1) 형광체 합성 및 공도핑 연구
(2) Thermoluminescence 기반 트랩 분석
(3) MTF 및 DQE 정량 평가
(4) Monte Carlo 기반 X-ray 상호작용 모델링

4. 기대 효과
- 의료 영상 저선량화
- 산업용 NDT 고해상도 검사
- 국산 CR 핵심 소재 기술 확보

5. 연구 기간
3년

6. 활용 계획
- SCI 논문 3편 이상
- 특허 출원 2건 이상
- 산업체 기술 이전 추진
```



