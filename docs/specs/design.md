# Design Specification: HCC Communication-Aware ICB Response Prediction

> **상태:** 확정 (2026-03-29)
> **저자:** 1인
> **타임라인:** 6개월 (2026-04 ~ 2026-09)
> **타겟 저널:** Genome Medicine / npj Precision Oncology / Cancer Research

---

## 1. 논문 가제

**"Cross-attention modeling of T cell exhaustion–macrophage communication networks predicts immunotherapy response in hepatocellular carcinoma"**

---

## 2. 가설

### 핵심 가설
> HCC의 ICB 반응은 개별 세포 유형의 abundance나 단일 유전자 발현이 아니라, **T cell exhaustion state와 TAM subtype 간의 ligand-receptor communication 패턴**에 의해 결정된다.

### 세부 가설
1. T cell exhaustion은 discrete states가 아닌 continuous spectrum이며, 각 state는 특이적인 TAM subtype과 우선적으로 소통한다
2. Terminal exhaustion을 유도하는 특정 TAM–T cell communication axis (예: SPP1+ TAM → terminal exhaustion, TREM2+ TAM → pre-exhaustion 유지)가 존재한다
3. 이 communication pattern을 cross-attention으로 정량화한 "Communication Score"가 기존 biomarker (TMB, PD-L1, TIDE, GEP)보다 ICB 반응을 더 정확히 예측한다
4. ICB responder에서는 치료 후 communication pattern이 rewiring되며, non-responder에서는 유지된다 (GSE206325로 검증)

---

## 3. 데이터셋

### 3-1. Discovery: scRNA-seq 통합

| Tier | 데이터셋 | 세포 수 | 환자 수 | 플랫폼 | 다운로드 | 저장 위치 |
|------|---------|--------|---------|--------|---------|---------|
| **1** | GSE149614 (Lu 2022) | ~70K | 10 | 10x v2 | 1.4 GB | `/mnt/e/.../raw/GSE149614/` |
| **1** | GSE140228 (Zhang 2019) | ~66K droplet | 16+ | 10x v2 | 305 MB | `/mnt/e/.../raw/GSE140228/` |
| **1** | GSE156625 (Sharma 2020) | ~212K | 다수 | 10x v2/v3 | 1.7 GB | `/mnt/e/.../raw/GSE156625/` |
| **2** | GSE151530 (Ma 2021) | ~57K | 46 | 10x v2 | 292 MB | `/mnt/e/.../raw/GSE151530/` |

**통합 예상:** ~70+ 환자, ~350K–410K 세포, T cell ~50K+, TAM ~25K+

**주의사항:**
- GSE140228: CD45+ sorted → 면역세포 subset 분석에만 사용, cell proportion 분석 제외
- GSE151530: HCC + iCCA 혼합 → Info.txt로 HCC만 필터
- GSE156625: Human + Mouse 혼합 → HCC-only h5ad 사용

### 3-2. Validation: Prognostic

| 데이터셋 | 타입 | 규모 | 용도 | 저장 위치 |
|---------|------|------|------|---------|
| TCGA-LIHC | RNA-seq | 371명 | Prognostic (KM, Cox), TIDE, deconvolution | `/mnt/e/.../external/TCGA_LIHC/` |
| GSE14520 | Microarray | 488명 | Independent survival validation | `/mnt/e/.../external/GSE14520/` |

### 3-3. Validation: ICB Response

| 우선순위 | 데이터셋 | 타입 | 규모 | 약제 | 반응 data | 저장 위치 |
|---------|---------|------|------|------|---------|---------|
| **1** | GSE285963 (Kanzaki 2025) | Bulk RNA-seq | 73명/84샘플 | Atezo+Bev | 논문 suppl.에서 추출 | `/mnt/e/.../ICB_cohorts/GSE285963/` |
| **2** | GSE206325 (Magen 2023) | scRNA-seq | 18명 | Cemiplimab | R vs NR (path. necrosis) | `/mnt/e/.../ICB_cohorts/GSE206325/` |
| **3** | GSE140901 (Hsu 2021) | NanoString 770 genes | 24명 | Nivo/Pembro+Ipi | GEO metadata에 직접 포함 | `/mnt/e/.../ICB_cohorts/GSE140901/` |

### 3-4. Reference (통합 미포함)

- GSE98638 (Zheng 2017): Smart-seq2, T cell exhaustion marker reference
- MacParland et al. (2018): 정상 간 reference

---

## 4. 분석 파이프라인

### Phase 1: Data Integration (Month 1)

```
1-1. 데이터 다운로드 및 QC
     ├── GEO에서 4개 데이터셋 다운로드 → /mnt/e/no_exp_paper/data/raw/
     ├── 각 데이터셋별 Scanpy로 로드
     ├── QC: nCount_RNA, nFeature_RNA, percent_mito, doublet detection (Scrublet)
     ├── GSE151530에서 HCC 환자만 필터링
     └── GSE156625에서 HCC-only h5ad 사용

1-2. 배치 통합
     ├── 방법: scVI (recommended) 또는 Harmony
     │   └── scVI 선호 이유: 대규모 데이터에 강건, probabilistic framework
     ├── Batch variable: dataset (GSE번호)
     ├── 통합 품질 평가: kBET, LISI, silhouette score
     └── 결과: integrated.h5ad → /mnt/e/.../processed/

1-3. Major cell type annotation
     ├── Leiden clustering → rough annotation
     ├── Marker-based: T cell (CD3D/E), Myeloid (CD68, CD14, LYZ),
     │   B cell (CD79A, MS4A1), NK (NKG7, GNLY), Fibroblast (COL1A1, DCN),
     │   Endothelial (PECAM1, VWF), Epithelial/Tumor (EPCAM, ALB)
     └── CellTypist 또는 scType으로 automated annotation → manual curation
```

### Phase 2: Cell State Discovery (Month 1–2)

```
2-1. T cell subset analysis
     ├── CD3+ T cell 분리 → re-clustering
     ├── CD4 vs CD8 분리 (CD4, CD8A/B)
     ├── CD8+ T cell exhaustion trajectory:
     │   ├── Marker genes: PDCD1, HAVCR2, LAG3, TIGIT, TOX, ENTPD1, LAYN, CTLA4
     │   ├── Pseudotime: Monocle3 (principal graph) + scVelo (RNA velocity)
     │   ├── State 정의:
     │   │   ├── Progenitor exhausted (TCF7+, PDCD1low)
     │   │   ├── Transitional exhausted (PDCD1+, HAVCR2-, TOX+)
     │   │   └── Terminal exhausted (PDCD1+, HAVCR2+, LAG3+, ENTPD1+)
     │   └── Gene module scoring (Tirosh et al. approach) + NMF for data-driven states
     └── TF activity: pySCENIC → exhaustion state별 regulatory program

2-2. TAM subset analysis
     ├── CD68+/CD14+/LYZ+ myeloid cell 분리 → re-clustering
     ├── Known TAM subtypes:
     │   ├── SPP1+ TAM (pro-tumorigenic, angiogenic)
     │   ├── TREM2+ TAM (lipid-associated, immunosuppressive)
     │   ├── FOLR2+ TAM (tissue-resident, variable)
     │   ├── FCN1+ TAM (inflammatory monocyte-derived)
     │   ├── C1Q+ TAM (complement-activating)
     │   └── ISG+ TAM (interferon-stimulated, rare)
     ├── Data-driven subtyping: consensus NMF + Leiden → novel states 탐색
     ├── Functional characterization: pathway scoring (GSVA), phagocytosis/antigen presentation score
     └── pySCENIC → TAM subtype별 TF regulatory program
```

### Phase 3: Communication Network Construction (Month 2–3)

```
3-1. Ligand-receptor interaction quantification
     ├── CellChat v2:
     │   ├── T cell exhaustion state (3+) × TAM subtype (5-6+) 간 L-R interaction
     │   ├── Patient-level communication strength matrix 계산
     │   └── Signaling pathway 수준 aggregation (e.g., PD-L1, TIGIT, GALECTIN, SPP1, CXCL)
     ├── NicheNet:
     │   ├── Upstream ligand → downstream target gene regulatory potential
     │   └── TAM ligand → T cell exhaustion gene program 활성화 정도
     └── LIANA+ (optional): 다중 method consensus

3-2. Communication tensor 구축
     ├── 3D tensor: [T cell state × TAM subtype × L-R pair] per patient/sample
     ├── 차원: ~3-5 T states × 5-6 TAM subtypes × 50-200 significant L-R pairs
     ├── 값: communication probability (CellChat) 또는 interaction strength
     └── Missing value 처리: 세포 수 부족 환자의 특정 state는 0 또는 imputation

3-3. Differential communication analysis
     ├── Tumor vs Normal 간 communication 차이
     ├── Exhaustion state별 preferential TAM partner 식별
     └── 환자 간 communication heterogeneity 평가
```

### Phase 4: Deep Learning Model — CACE (Month 3–4)

```
4-1. Architecture: Cross-Attention Communication Encoder (CACE)

     ┌─────────────────────────────────────────────────────┐
     │                    CACE Architecture                 │
     │                                                     │
     │  Input: L-R interaction tensor per sample            │
     │  [T_states × TAM_subtypes × LR_pairs]              │
     │                                                     │
     │  ┌─────────────┐     ┌──────────────┐              │
     │  │ T cell state │     │ TAM subtype  │              │
     │  │  Encoder     │     │  Encoder     │              │
     │  │ (MLP + LN)   │     │ (MLP + LN)   │              │
     │  └──────┬───────┘     └──────┬───────┘              │
     │         │ Query               │ Key, Value          │
     │         └──────┬──────────────┘                     │
     │                ▼                                     │
     │  ┌──────────────────────┐                           │
     │  │  Multi-Head Cross     │                           │
     │  │  Attention (h=4-8)    │  ← 핵심: 어떤 TAM이      │
     │  │                       │    어떤 T state와         │
     │  └──────────┬────────────┘    소통하는지 학습        │
     │             ▼                                        │
     │  ┌──────────────────────┐                           │
     │  │  Communication        │                           │
     │  │  Embedding (d=64-128) │                           │
     │  └──────────┬────────────┘                           │
     │             ▼                                        │
     │  ┌────────┐  ┌────────────┐                         │
     │  │Prognost│  │ICB Response│                         │
     │  │Head    │  │Head        │                         │
     │  │(Cox PH)│  │(BCE)       │                         │
     │  └────────┘  └────────────┘                         │
     └─────────────────────────────────────────────────────┘

4-2. Training strategy
     ├── Prognostic head: TCGA-LIHC (n=371) + GSE14520 (n=488)
     │   ├── Input: CIBERSORTx로 bulk → cell state proportion 추정
     │   │          → scRNA-seq 기반 L-R reference로 communication score 재구성
     │   ├── Loss: Cox partial likelihood loss
     │   └── 5-fold CV
     ├── ICB response head: GSE285963 (n=73)
     │   ├── Input: bulk RNA-seq → deconvolution → communication score
     │   ├── Loss: Binary cross-entropy (responder vs non-responder)
     │   └── Leave-one-out CV (소규모 데이터)
     └── Multi-task learning: 두 head 동시 학습 (weighted loss)

4-3. scRNA-seq → Bulk transfer 전략
     ├── Step 1: scRNA-seq에서 cell state별 L-R signature 정의 (reference)
     ├── Step 2: Bulk RNA-seq에서 CIBERSORTx로 cell state fraction 추정
     ├── Step 3: Fraction × reference L-R strength = pseudo-communication matrix
     └── Step 4: CACE에 pseudo-communication matrix 입력

4-4. Interpretability
     ├── Attention weight analysis → 핵심 TAM–T cell communication axis 식별
     ├── SHAP / Integrated Gradients → L-R pair importance ranking
     ├── Communication Score: embedding의 scalar projection (PCA component 1 또는 learned)
     └── Biological validation: top L-R pairs의 known function과 일치 여부
```

### Phase 5: Validation (Month 4–5)

```
5-1. Prognostic validation
     ├── TCGA-LIHC (n=371):
     │   ├── Communication Score → KM survival analysis (high vs low, median split)
     │   ├── Multivariate Cox: Communication Score + stage + AFP + etiology
     │   ├── Time-dependent AUC (1yr, 3yr, 5yr OS)
     │   └── 기존 signature 비교: TIDE, ESTIMATE, GEP, TMB
     └── GSE14520 (n=488):
         └── Independent cohort KM + Cox (microarray → communication score)

5-2. ICB response validation
     ├── GSE285963 (n=73, RNA-seq):
     │   ├── Communication Score → ROC-AUC for response prediction
     │   ├── 기존 biomarker 대비 비교: TIDE, GEP, TMB
     │   └── Pre vs Post-treatment communication pattern 변화 (matched samples)
     ├── GSE206325 (n=18, scRNA-seq):
     │   ├── 직접 scRNA-seq 수준에서 communication pattern 계산
     │   ├── Responder vs Non-responder 간 communication 차이
     │   ├── ICB 전후 communication rewiring 관찰
     │   └── Attention weight 해석과 일치 여부 (CACE가 중요하다고 한 L-R axis가 실제로 변하는지)
     └── GSE140901 (n=24, NanoString):
         └── 770 immune genes 내에서 가능한 범위의 보조 검증

5-3. Spatial validation (if data available)
     ├── 공개 HCC Visium 데이터
     ├── Exhausted T cell ↔ TAM co-localization
     └── Communication axis의 spatial proximity 확인

5-4. Ablation study
     ├── Cross-attention 제거 → simple concatenation 대비 성능 차이
     ├── Communication feature vs cell abundance feature 단독 비교
     ├── T cell state 수 / TAM subtype 수 변화에 따른 robustness
     └── Random L-R pair permutation → communication specificity 확인
```

### Phase 6: Manuscript (Month 5–6)

```
6-1. Figure 구성 (8 Main + Supplementary)

Fig 1: Study overview
├── (A) Graphical abstract: 데이터 통합 → cell state → communication → CACE → prediction
├── (B) scRNA-seq 통합 UMAP (major cell types)
├── (C) 데이터셋 구성 summary (환자 수, 세포 수, 조직 유형)
└── (D) Validation cohort overview

Fig 2: T cell exhaustion continuum in HCC
├── (A) CD8+ T cell UMAP with exhaustion states
├── (B) Pseudotime trajectory (progenitor → transitional → terminal)
├── (C) Exhaustion state marker gene heatmap
├── (D) SCENIC: state-specific TF activity (TOX, TCF7, BATF, etc.)
└── (E) State proportion across patients (heterogeneity)

Fig 3: TAM subtype landscape
├── (A) Myeloid UMAP with TAM subtypes
├── (B) Subtype marker gene dot plot
├── (C) Functional pathway scoring (phagocytosis, angiogenesis, immunosuppression)
├── (D) pySCENIC: subtype-specific TF programs
└── (E) Tumor vs Normal 간 TAM composition shift

Fig 4: T cell–TAM communication network
├── (A) CellChat 결과: state-specific L-R interaction heatmap
├── (B) Sankey diagram: T cell state → L-R pair → TAM subtype
├── (C) 핵심 communication axis 개별 시각화 (e.g., SPP1–CD44, TIGIT–PVR)
├── (D) NicheNet: TAM ligand → T cell exhaustion gene regulatory potential
└── (E) Patient-level communication heterogeneity

Fig 5: CACE model architecture and performance
├── (A) Architecture schematic
├── (B) Training curves (loss, AUC)
├── (C) Ablation study results
└── (D) Comparison vs baseline models (XGBoost, simple MLP, cell abundance only)

Fig 6: Prognostic validation
├── (A) TCGA-LIHC: Communication Score KM plot
├── (B) Multivariate Cox forest plot
├── (C) Time-dependent AUC curves (vs TIDE, ESTIMATE, GEP)
├── (D) GSE14520: Independent KM validation
└── (E) Subgroup analysis (etiology: HBV vs HCV vs NASH)

Fig 7: ICB response prediction
├── (A) GSE285963: Communication Score → ROC curve
├── (B) Waterfall plot: Communication Score vs response
├── (C) 기존 biomarker 대비 AUC 비교
├── (D) Attention weight visualization: top communication axes
└── (E) SHAP: L-R pair importance for ICB prediction

Fig 8: Mechanistic validation
├── (A) GSE206325: Responder vs NR communication pattern 차이 (scRNA-seq)
├── (B) ICB 전후 communication rewiring (matched pre/post)
├── (C) Top L-R axis의 실제 발현 변화 (violin plots)
├── (D) Spatial co-localization (if available)
└── (E) Proposed mechanism model (graphical summary)

6-2. 논문 구조
├── Abstract (250 words)
├── Introduction
├── Results (Fig 1-8에 대응하는 8개 section)
├── Discussion
├── Methods (상세 computational methods)
├── Data & Code Availability
└── Supplementary (추가 figures, tables, 모델 hyperparameters)

6-3. 코드/데이터 공개
├── GitHub repository: 전체 분석 코드
├── Zenodo: processed data snapshot + model weights
└── Interactive portal (optional): Communication Score calculator
```

---

## 5. 기술 스택

| 카테고리 | 도구 |
|---------|------|
| **scRNA-seq 분석** | Scanpy, scVI, Scrublet, pySCENIC |
| **Trajectory** | Monocle3 (R) 또는 scVelo (Python) |
| **Communication** | CellChat v2 (R), NicheNet (R), LIANA+ (Python) |
| **Deconvolution** | CIBERSORTx, BayesPrism, 또는 SCDC |
| **DL Framework** | PyTorch |
| **ML 비교** | XGBoost, scikit-learn (Random Survival Forest) |
| **Interpretability** | SHAP, captum (Integrated Gradients) |
| **통계** | lifelines (survival), scipy, statsmodels |
| **시각화** | matplotlib, seaborn, scanpy.pl, CellChat plot |

---

## 6. 타임라인

| 기간 | Phase | 주요 산출물 |
|------|-------|-----------|
| **Month 1** (2026-04) | Data Integration | integrated.h5ad, annotated cell types |
| **Month 1–2** (2026-04~05) | Cell State Discovery | T cell exhaustion states, TAM subtypes, SCENIC results |
| **Month 2–3** (2026-05~06) | Communication Network | CellChat/NicheNet results, communication tensor |
| **Month 3–4** (2026-06~07) | CACE Model | Trained model, communication score |
| **Month 4–5** (2026-07~08) | Validation | All validation results, figures draft |
| **Month 5–6** (2026-08~09) | Manuscript | Submission-ready manuscript |

---

## 7. 리스크 및 대응

| 리스크 | 심각도 | 대응 |
|--------|--------|------|
| scVI integration 품질 부족 | 🟡 | Harmony, scANVI 대안; batch별 cell type 비율 비교로 품질 평가 |
| T cell exhaustion state 구분 불명확 | 🟡 | Gene module scoring + NMF 병행; binary (exhausted/non) fallback |
| Communication score의 bulk 전환 정확도 | 🟡 | CIBERSORTx 외 BayesPrism 비교; ablation에서 전환 방법 간 비교 |
| CACE 모델 overfitting (소규모 ICB 데이터) | 🟡 | LOOCV, 강한 regularization, pre-training on TCGA → fine-tune on ICB |
| GSE285963 response label 추출 실패 | 🟢 | 저자 contact; TIDE score로 대체 가능 |
| 가설 기각 (communication이 abundance보다 못 함) | 🟡 | Negative result도 가치 있음; 생물학적 발견(cell state atlas)으로 방향 전환 |
| 외장하드 접근 불안정 | 🟢 | 핵심 processed 파일은 로컬 백업; 스크립트로 재생성 가능하게 |

---

## 8. 성공 기준

| 기준 | 목표 |
|------|------|
| scRNA-seq 통합 세포 수 | >200K cells, >50 환자 |
| T cell exhaustion state 수 | 3개 이상 biologically distinct states |
| TAM subtype 수 | 5개 이상 |
| Communication Score prognostic value | TCGA-LIHC에서 KM p < 0.001, multivariate Cox에서 independent |
| ICB response prediction AUC | GSE285963에서 >0.70 (기존 TIDE 대비 우위) |
| Ablation: communication vs abundance | Communication feature가 단독 abundance보다 AUC 0.05+ 향상 |
| 논문 투고 | 2026-09까지 1차 투고 완료 |
