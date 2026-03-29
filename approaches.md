# 논문 접근법 제안서

> 결정: **Approach 1 채택**
> 날짜: 2026-03-28

---

## Approach 1 (채택): Cell-Cell Communication Network 기반 ML Framework로 HCC ICB 반응 예측

### 핵심 아이디어
HCC 공개 scRNA-seq 데이터를 통합하여 T cell exhaustion states와 TAM subtypes를 고해상도로 규명한 뒤, 이 둘 사이의 ligand-receptor communication network를 GNN(Graph Neural Network)으로 모델링하여 "communication-aware" ICB 반응 예측 signature를 발굴

### 왜 Novel한가
- 기존 연구들은 cell type abundance 또는 gene expression 기반 signature → 세포 간 **상호작용 패턴** 자체를 signature로 쓰는 연구는 거의 없음
- T cell exhaustion은 연속적 스펙트럼인데, 각 state별로 어떤 TAM subtype과 소통하는지 체계적으로 매핑한 HCC 연구 부재
- GNN을 cell communication에 적용하는 방법론적 novelty + 생물학적 발견의 균형

### 분석 파이프라인
1. HCC scRNA-seq 통합 (5-8개 공개 데이터셋, ~20-30만 세포)
2. T cell exhaustion trajectory (Monocle3/scVelo) → pre-, early-, terminal-exhaustion states
3. TAM subtype annotation (SPP1+, TREM2+, FOLR2+, FCN1+ 등)
4. CellChat/NicheNet L-R interaction → cell communication graph 구축
5. GNN으로 communication embedding 추출 → "Communication Score"
6. TCGA-LIHC + ICB bulk cohorts에서 deconvolution 기반 검증

### 타겟 저널
Genome Medicine, npj Precision Oncology, Cancer Research

---

## Approach 2 (대안): T Cell Exhaustion Continuum + TAM Interaction Atlas in HCC

### 핵심 아이디어
HCC에서 T cell exhaustion의 multi-resolution continuum을 규명하고, 각 exhaustion state와 상호작용하는 TAM subtype을 매핑하여 state-specific prognostic signature 발굴

### 왜 Novel한가
- HCC T cell exhaustion 연구는 있지만, exhaustion을 binary가 아닌 continuous spectrum으로 보고 각 state별 TAM 파트너를 체계적으로 규명한 연구 부재
- Biology 비중이 높아 Cancer Research, Clinical Cancer Research에 적합

### 분석 파이프라인
1. HCC scRNA-seq 통합 → T cell subclustering + pseudotime
2. Exhaustion state 정의 (gene module scoring, NMF)
3. TAM subtype 규명 + spatial co-localization (공개 ST 데이터 활용)
4. State-specific L-R interaction mapping
5. ML (XGBoost/RSF)로 multi-state prognostic signature
6. TCGA-LIHC + GEO bulk cohorts 검증

### 타겟 저널
Cancer Research, Clinical Cancer Research

---

## Approach 3 (대안): 범용 scRNA-seq → Communication Signature Framework (HCC discovery, CRC validation)

### 핵심 아이디어
scRNA-seq에서 cell communication network를 자동 추출하고 prognostic/predictive signature로 변환하는 general-purpose framework를 개발. HCC에서 발굴, CRC에서 cross-cancer 검증

### 왜 Novel한가
- 프레임워크 자체가 contribution → 재현성과 범용성 강조
- 두 암종에서 검증하여 일반화 가능성 입증
- 가장 methodology novelty가 높음

### 단점
- 생물학적 깊이가 분산될 수 있음
- 두 암종 모두 데이터 수집·처리해야 하므로 1인 작업에 부담

### 타겟 저널
Genome Biology, Briefings in Bioinformatics

---

## 비교표

| | Approach 1 | Approach 2 | Approach 3 |
|---|---|---|---|
| **방법론 novelty** | ★★★★ | ★★ | ★★★★★ |
| **생물학적 깊이** | ★★★★ | ★★★★★ | ★★★ |
| **1인 실현 가능성** | ★★★★ | ★★★★★ | ★★★ |
| **임팩트 잠재력** | ★★★★★ | ★★★★ | ★★★★ |
| **차별화** | ★★★★★ | ★★★ | ★★★★ |
