# Reference Papers: Computational-Only Biomarker Discovery via Public scRNA-seq Reanalysis

> 2025년 이후(또는 근접) 출판된 논문 중, wet-lab 실험 없이 공개 scRNA-seq 데이터 재분석으로 novel tumor/immune biomarker를 발굴한 peer-reviewed 논문 목록

---

## 1. TabulaTIME — Nature Cancer (2025.08)

**Han et al.** "Spatiotemporal analyses of the pan-cancer single-cell landscape reveal widespread profibrotic ecotypes associated with tumor immunity"

| 항목 | 내용 |
|------|------|
| **데이터** | 103개 공개 scRNA-seq 연구, 746명 환자, 36개 암종, ~450만 세포 |
| **방법** | MAESTRO workflow, MetaCell, CCA integration, spatial transcriptomics 통합 |
| **핵심 바이오마커** | CTHRC1+ CAF 및 SLPI+ macrophage profibrotic ecotype → 면역배제(immune exclusion) 및 불량 예후 |
| **특이점** | 순수 공개데이터 메타분석으로 Nature Cancer 게재. Pan-cancer 리소스(TabulaTIME) 구축 |

---

## 2. Pan-cancer TME Atlas — Cell Reports Medicine (2025.10)

**Lodi et al.** "Decoding tumor heterogeneity: A spatially informed pan-cancer analysis of the tumor microenvironment"

| 항목 | 내용 |
|------|------|
| **데이터** | 9개 암종, 230개 treatment-naïve 샘플 scRNA-seq + TCGA bulk RNA-seq (7,493 샘플, 14개 암종) |
| **핵심 바이오마커** | TLS-like hub 및 immune-reactive hub score → checkpoint immunotherapy 반응 예측 |
| **검증** | CD8+ exhausted T세포 시그니처와 높은 상관관계 (R=0.84–0.89) |
| **실험 여부** | 완전히 computational only |

---

## 3. Pan-cancer DC Atlas — Cancer Research (2025.10)

**Zemin Zhang 그룹** "Pan-Cancer Analyses Refine the Single-Cell Portrait of Tumor-Infiltrating Dendritic Cells"

| 항목 | 내용 |
|------|------|
| **데이터** | 158개 공개 10x Genomics 데이터셋, 1,507명 donor, 33개 암종, 2,500+ 샘플 |
| **방법** | CellTypist 기반 DC subtype classifier 개발, ensemble prediction |
| **핵심 바이오마커** | AXL+SIGLEC6+ DC, Langerhans cell-like DC → LAMP3+ DC의 추가적 기원, 암종별 예후 및 면역치료 반응 연관 |
| **특이점** | Cancer Research (AACR) 게재, 순수 재분석 |

---

## 4. Pan-cancer Tumor-Normal Ecosystem — Nature Communications (2025.03)

**저자 미상** "Systematic dissection of tumor-normal single-cell ecosystems across a thousand tumors of 30 cancer types"

| 항목 | 내용 |
|------|------|
| **데이터** | 1,070 종양 + 493 정상 샘플 (~490만 세포), 137 spatial transcriptomics, 8,887 TCGA, 1,261 ICB-treated bulk |
| **핵심 바이오마커** | AKR1C1+ vs WNT5A+ inflammatory fibroblast 구분, interferon-enriched community → ICB 반응 예측 |
| **특이점** | 역대 최대급 pan-cancer 통합 분석, 실험 없음 |

---

## 5. Myeloid Cell States as Prognostic Markers — Nature Communications (2024.07)

**저자 미상** "Single-cell resolution characterization of myeloid-derived cell states with implication in cancer outcome"

| 항목 | 내용 |
|------|------|
| **데이터** | 여러 암종 공개 scRNA-seq 통합 |
| **핵심 바이오마커** | TREM2+PD-1+ macrophage, FOLR2+PDL-2+ macrophage → 독립적 예후인자 |
| **검증** | 독립 코호트에서 FOLR2-발현 macrophage의 ovarian/TNBC 불량 예후 연관 확인 |
| **방법** | CIBERSORTx deconvolution 기반, 실험 없음 |

---

## 6. ML + scRNA-seq ICI Response Prediction — npj Precision Oncology (2025.04)

**저자 미상** "Uncovering gene and cellular signatures of immune checkpoint response via machine learning and single-cell RNA-seq"

| 항목 | 내용 |
|------|------|
| **데이터** | 공개 scRNA-seq (melanoma) + bulk cohort 검증 |
| **핵심 바이오마커** | 11-gene ICI response signature (pan-cancer 적용 가능) |
| **방법** | XGBoost (AUC 0.84→0.89), Boruta feature selection, SHAP value, 강화학습 기반 informative cell selection |
| **특이점** | 방법론적 차별점 — reinforcement learning으로 세포 선별 |

---

## 7. PDAC Single-Cell Atlas — Clinical Cancer Research (2025.02)

**Loveless, Steele et al.** "Human Pancreatic Cancer Single-Cell Atlas Reveals Association of CXCL10+ Fibroblasts and Basal Subtype Tumor Cells"

| 항목 | 내용 |
|------|------|
| **데이터** | 229명 환자 공개 scRNA-seq (~70만 세포) |
| **핵심 바이오마커** | CXCL10+ CAF – basal tumor cell niche → PDAC 예후 불량 |
| **검증** | Spatial transcriptomics (n=22) + bulk RNA-seq (n=744) + multiplex immunostaining |
| **저널** | Clinical Cancer Research (AACR, IF ~13) |

---

## 8. CRC Spatial + scRNA-seq Prognostic Signature — 2025

**저자 미상** "Spatial transcriptomics and scRNA-seq: decoding tumor complexity and constructing prognostic models in colorectal cancer"

| 항목 | 내용 |
|------|------|
| **데이터** | 공개 scRNA-seq + spatial transcriptomics |
| **핵심 바이오마커** | MLXIPL+ neoplasm subtype, 13-gene PS → CRC 예후 및 치료 반응 예측 |
| **방법** | Machine learning (StepCox backward), 순수 bioinformatics |
| **특이점** | Limitation에서 wet-lab 실험 없음을 명시 |

---

## 9. Ferroptosis-Related Biomarkers in TNBC — eLife (2025.06)

**Gong, Gu et al.** "Ferroptosis-related genes mediate tumor microenvironment and prognosis in triple-negative breast cancer via integrated RNA-seq analysis"

| 항목 | 내용 |
|------|------|
| **데이터** | GSE176078 (9 TNBC scRNA-seq), GSE25066 (178 TNBC bulk, training), GSE86166 (validation), 471 ferroptosis genes from FerrDb |
| **핵심 바이오마커** | 23-gene ferroptosis risk score (TMEM160, EWSR1, BCAT2 등 protective; CDC25B, HMGCS1, SPC25, PTTG1 등 risk) |
| **검증** | 외부 코호트 AUC 0.65-0.71, SCENIC TF regulatory network 분석 |
| **실험 여부** | 완전히 computational — 논문에서 명시 |

---

## 10. Prostate Cancer Epithelial Gene Signature — Molecular Oncology (2025.02)

**Mou & Harries** "Integration of single-cell and bulk RNA-sequencing data reveals the prognostic potential of epithelial gene markers for prostate cancer"

| 항목 | 내용 |
|------|------|
| **데이터** | TCGA-PRAD (n=547), GSE21034 (n=160), GSE70768 (n=199), E-MTAB-6128 (n=141), DKFZ (n=118), GSE193337 scRNA-seq |
| **핵심 바이오마커** | 11-gene epithelial signature (APEX1, CACNA1D, ERBB3, FASN, FBP1, FKBP4, LAMTOR2, PHB2, PSMG4, RAB25, TMED3) |
| **방법** | 97개 ML 알고리즘 비교, Random Survival Forest 최적, 68개 기존 signature 능가 (Decipher, Prolaris, GPS 포함) |
| **실험 여부** | 완전히 computational |

---

## 공통 패턴 분석

### 성공적인 computational-only 논문의 핵심 요소
1. **대규모 공개 데이터 통합** — 다수 연구의 scRNA-seq 데이터를 체계적으로 수집·통합
2. **Novel cell state/subtype 발견** — 기존에 알려지지 않은 세포 아집단이나 상태 규명
3. **다층적 검증** — TCGA bulk RNA-seq, spatial transcriptomics, ICB 코호트 등으로 cross-validation
4. **임상적 relevance** — 예후 예측, 치료 반응 예측 등 translational value 제시
5. **리소스 구축** — Atlas, classifier, signature 등 커뮤니티가 재사용할 수 있는 도구/데이터 제공
