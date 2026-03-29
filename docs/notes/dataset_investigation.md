# HCC scRNA-seq 공개 데이터셋 조사

> 작성일: 2026-03-28 (GEO 실제 확인 완료)
> 목적: Approach 1 (Communication-aware DL framework) 수행을 위한 데이터셋 적합성 평가
> 상태: GEO/CNGBdb 직접 확인 기반 정보 (2026-03-28 기준)

---

## 1. scRNA-seq 데이터셋 상세 (GEO 확인 완료)

### GSE149614 — Lu et al., Nature Communications (2022)

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | A Single-Cell Atlas of the Multicellular Ecosystem of Primary and Metastatic HCC |
| **PMID** | 35933472 |
| **공개일** | 2020-04-30 |
| **플랫폼** | GPL24676 — Illumina NovaSeq 6000 (10x Chromium 3' V2) |
| **샘플 수** | 21 GSM entries (GSM4505944–GSM4505964), 10명 환자 |
| **세포 수** | ~70,000 cells |
| **조직** | Primary tumor + non-tumor liver + PVTT + metastatic lymph node |
| **Supplementary files** | `HCC.scRNAseq.S71915.count.txt.gz` (158 MB), `normalized.txt.gz` (1.2 GB), `metadata.updated.txt.gz` (478 KB) |
| **다운로드 크기** | ~1.4 GB |
| **Raw FASTQ** | ⚠️ EGA (EGAS00001004468) — controlled access. GEO에는 processed matrix만 |
| **접근성** | ✅ Processed data 공개 |

---

### GSE140228 — Zhang et al., Cell (2019)

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | Landscape and Dynamics of Single Immune Cells in HCC |
| **PMID** | 31675496 |
| **공개일** | 2019-11-13 |
| **플랫폼** | GPL20301 — HiSeq 4000 (10x Droplet + Smart-seq2) |
| **샘플 수** | 41 GSM entries, CD45+ sorted |
| **세포 수** | ~66,187 droplet + ~7,099 Smart-seq2 = ~73,000 cells |
| **조직** | Tumor + adjacent liver + blood + hepatic lymph node + ascites (5개 조직 부위) |
| **Supplementary files** | Droplet: `UMI_counts_Droplet.mtx.gz` (250 MB) + barcodes/genes/cellinfo; Smart-seq2: `read_counts_Smartseq2.csv.gz` (53 MB) + cell_info/gene_info |
| **다운로드 크기** | ~305 MB |
| **Raw FASTQ** | China GSA (HRA000069) — not SRA |
| **접근성** | ✅ Processed data 공개 |
| **주의** | CD45+ sorted → 면역세포 편향. Droplet과 Smart-seq2 분리 제공 |

---

### GSE156625 — Sharma et al., Cell (2020) ⚠️ 기존 오류 수정

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | Onco-fetal reprogramming of endothelial cells drives immunosuppressive macrophages in HCC |
| **PMID** | 32976798 |
| **⚠️ 주의** | **기존에 Sun et al. 2021로 잘못 기재했으나, 실제로는 Sharma et al. 2020 논문** |
| **공개일** | 2020-09-24 |
| **플랫폼** | GPL16791 / GPL17021 — HiSeq 2500 (10x Chromium 3' v2/v3), Human + Mouse |
| **샘플 수** | 78 GSM entries (human HCC + fetal liver + normal liver + mouse) |
| **세포 수** | ~212,000 total cells |
| **조직** | HCC tumor + adjacent normal + fetal liver (human) + mouse liver |
| **Supplementary files** | HCC-only: `HCCmatrix.mtx.gz` (562 MB) + `HCCscanpyobj.h5ad.gz` (800 MB); Combined: `HCCFmatrix.mtx.gz` (687 MB) + `HCCFscanpyobj.h5ad.gz` (1.1 GB) + `HCCF_integrated.RData.gz` (1.6 GB) |
| **다운로드 크기** | ~1.7 GB (HCC only), ~6.5 GB (전체) |
| **접근성** | ✅ 완전 공개. **h5ad + RData 제공 — 즉시 사용 가능** |
| **핵심 가치** | Endothelial–macrophage crosstalk 연구 → TAM subtyping에 직접 관련 |

---

### GSE151530 — Ma et al., J Hepatology (2021)

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | Understanding tumor evolution and its clonality by single-cell transcriptomic analysis in liver cancer |
| **PMID** | 34216724 |
| **공개일** | 2021-02-07 |
| **플랫폼** | GPL18573/GPL20301/GPL24676 — NextSeq 500 / HiSeq 4000 / NovaSeq 6000 (모두 10x 3' V2) |
| **샘플 수** | 46 GSM entries (GSM4581240–GSM4581285) |
| **세포 수** | ~57,000 cells |
| **조직** | Tumor samples (HCC + iCCA 혼합) |
| **Supplementary files** | `matrix.mtx.gz` (291 MB) + `barcodes.tsv.gz` + `genes.tsv.gz` + `Info.txt.gz` |
| **다운로드 크기** | ~292 MB |
| **접근성** | ✅ 완전 공개. 표준 10x MTX 포맷 |
| **주의** | ⚠️ HCC + iCCA 혼합 — Info.txt로 HCC만 필터링 필요. Pre/post-treatment 포함 |

---

### GSE125449 — Ma et al., Cancer Cell (2019)

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | Tumor Cell Biodiversity Drives Microenvironmental Reprogramming in Liver Cancer |
| **PMID** | 31588021 |
| **공개일** | 2019-10-06 |
| **플랫폼** | GPL18573/GPL20301 — NextSeq 500 / HiSeq 4000 (10x 3' V2) |
| **샘플 수** | 19 GSM entries (9 HCC + 10 iCCA) |
| **세포 수** | ⚠️ 매우 적음 (일부 샘플 124 cells) |
| **Supplementary files** | Set1: matrix.mtx.gz (25 MB); Set2: matrix.mtx.gz (18 MB) |
| **다운로드 크기** | ~43 MB |
| **접근성** | ✅ 완전 공개 |
| **주의** | GSE151530의 초기 버전 (같은 PI). 세포 수 극히 적어 통합 가치 낮음 |

---

### CNP0000650 — Sun et al., Cell (2021) ⚠️ GEO 아님

| 항목 | 확인 결과 |
|------|---------|
| **제목** | Single-cell landscape of the ecosystem in early-relapse HCC |
| **PMID** | 33357445 |
| **저장소** | ⚠️ **CNGBdb (China National GeneBank)** — GEO에 없음 |
| **URL** | https://db.cngb.org/search/project/CNP0000650 |
| **환자 수** | 19명 (primary + relapse) |
| **세포 수** | ~20,000 cells |
| **접근성** | CNGBdb에서 다운로드 가능하나, GEO 대비 접근 불편 |
| **판정** | ⚠️ 접근성 이슈로 Tier 2 이하로 하향. 필수는 아님 |

---

## 2. 검증용 Bulk 데이터셋 (GEO 확인 완료)

### GSE14520 — Roessler/Wang LCI Cohort

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | Gene expression data of human HCC |
| **공개일** | 2010-12-01 |
| **플랫폼** | GPL571 / GPL3921 — Affymetrix U133A (⚠️ **Microarray, RNA-seq 아님**) |
| **샘플 수** | 488 GSM entries (paired tumor + non-tumor) |
| **임상 metadata** | ✅ 풍부 — survival, staging, etiology |
| **다운로드** | `RAW.tar` (1.0 GB) + clinical supplement |
| **접근성** | ✅ 완전 공개 |
| **가치** | HCC 분야에서 가장 많이 사용되는 bulk 검증 코호트. Microarray이지만 ~20K gene 커버 |

### TCGA-LIHC

| 항목 | 내용 |
|------|------|
| **규모** | 371명 HCC, RNA-seq + WES + clinical |
| **접근** | GDC/cBioPortal ✅ 공개 |
| **임상 metadata** | OS, DFS, stage, grade, etiology ✅ |
| **가치** | Prognostic validation의 gold standard |

---

## 3. ICB Response 데이터 (2026-03-28 확인 완료)

> 🟢 이전에 🔴 최대 리스크로 분류했으나, GSE285963 확인으로 리스크 대폭 완화

### 3-1. GSE285963 — Kanzaki et al., JHEP Reports (2025) ⭐ 핵심 ICB 검증 데이터

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | Transcriptome profile of HCC in patients with atezolizumab plus bevacizumab therapy |
| **PMID** | 41362708 |
| **공개일** | 2025-10-02 |
| **플랫폼** | GPL24676 — Illumina NovaSeq 6000 (**Bulk RNA-seq, whole-transcriptome**) |
| **Library prep** | SMART-Seq Stranded Kit (Takara Bio) |
| **환자 수** | 73명 (84 샘플 — pre/post/on-treatment 포함) |
| **치료** | Atezolizumab + Bevacizumab (1L advanced HCC) |
| **반응 annotation** | ⚠️ GEO metadata에는 없음 → **논문 supplementary에서 추출 필요** (disease control, clinical benefit, PFS) |
| **생존 데이터** | PFS (논문에서) |
| **Supplementary** | `GSE285963_CU24_n84_readcount.gct.gz` (~3 MB) |
| **접근성** | ✅ **완전 공개 (GEO)** |
| **BioProject** | PRJNA1207067 |

**평가:** ⭐ **최우선 ICB 검증 데이터** — Whole-transcriptome RNA-seq + atezo+bev + 73명. Communication score 검증에 이상적. Response label은 논문 Table/Supplement에서 매핑 필요.

**TODO:**
- [ ] 논문 (PMID 41362708) supplementary에서 환자별 response annotation 추출
- [ ] Sample ID (CU24_XXpre/post) ↔ 논문 patient ID 매핑 테이블 작성
- [ ] Pre-treatment 샘플만 분리하여 ICB 반응 예측 검증에 사용

---

### 3-2. GSE206325 — Magen et al., Nature Medicine (2023) ⭐ ICB scRNA-seq

| 항목 | GEO 확인 결과 |
|------|-------------|
| **제목** | Intratumoral mregDC and Tfh cooperate for CD8 T cell reinvigoration after PD-1 blockade |
| **PMID** | 37322116 |
| **플랫폼** | 10x Genomics 5' scRNA-seq + CITE-seq + TCR-seq |
| **환자 수** | 18명 (neoadjuvant anti-PD-1, cemiplimab) |
| **반응** | Responder (n=6) vs Non-responder (n=12), pathological necrosis 기준 |
| **샘플** | Pre-treatment biopsy + post-treatment resection (matched) |
| **세포 수** | ~1,000,000 immune cells |
| **접근성** | ✅ **완전 공개 (GEO)** |

**평가:** ⭐ **메커니즘 검증에 최적** — ICB 전후 matched scRNA-seq + response annotation. Communication pattern 변화를 직접 관찰 가능. 단, 환자 수 적음(n=18).

---

### 3-3. GSE140901 — Hsu et al., Liver Cancer (2021)

| 항목 | GEO 확인 결과 |
|------|-------------|
| **플랫폼** | NanoString nCounter PanCancer Immune Panel (770 genes) |
| **환자 수** | 24명 |
| **치료** | Nivolumab mono 또는 Pembrolizumab + Ipilimumab |
| **반응 annotation** | ✅ **GEO metadata에 직접 포함** — best_response (PD/SD/PR), clinical_benefit (Y/N) |
| **생존** | ✅ PFS time/event, OS time/event, age, gender, etiology (HBV/HCV/Other) |
| **접근성** | ✅ 완전 공개 |

**평가:** 보조 검증용. Metadata가 가장 풍부하나 770 gene 제한 + 소규모.

---

### 3-4. Controlled Access (EGA)

| 데이터셋 | 타입 | 규모 | 약제 | EGA |
|---------|------|------|------|-----|
| **IMbrave150** (Zhu 2022) | Bulk RNA-seq | ~311명 RNA-seq | Atezo+Bev vs Sorafenib | EGAS00001005503 |
| **Haber et al. 2023** | Microarray | 83명 | Anti-PD-1 | EGAS00001005477 |

---

### 3-5. ICB 검증 전략 (수정)

| 우선순위 | 데이터셋 | 용도 | 접근성 |
|---------|---------|------|--------|
| **1** | **GSE285963** (n=73, RNA-seq) | Communication score → ICB 반응 예측 (primary ICB validation) | ✅ 즉시 |
| **2** | **GSE206325** (n=18, scRNA-seq) | ICB 전후 communication pattern 변화 메커니즘 검증 | ✅ 즉시 |
| **3** | **GSE140901** (n=24, NanoString) | 보조 검증 (immune gene subset) | ✅ 즉시 |
| **4** | TCGA-LIHC (n=371) | TIDE/ESTIMATE 간접 ICB 반응 추정 | ✅ 즉시 |
| **5** | IMbrave150 (n=311, EGA) | 최대 규모 검증 (추가 optional) | ⚠️ 신청 필요 |

**결론:** GSE285963 (whole-transcriptome, n=73, open) + GSE206325 (scRNA-seq, n=18, open)로 **ICB 검증이 공개 데이터만으로 충분히 가능**. EGA 신청은 optional bonus.

---

## 4. Spatial Transcriptomics 데이터

| 검색 키워드 | 결과 |
|-----------|------|
| "hepatocellular carcinoma Visium" on GEO | GSE203612 등 소수 존재 — 상세 확인 필요 |
| Guilliams et al. 2022 (Cell) | 정상 간 spatial atlas — reference용 |
| 10x Genomics 데모 데이터 | HCC Visium 데모 존재 가능 |

> Spatial 데이터는 추가 GEO 검색 필요

---

## 5. 수정된 통합 적합성 종합표

| 데이터셋 | Tier | 세포 수 | 환자 수 | 10x 호환 | Unsorted | 다운로드 | 판정 |
|---------|------|--------|---------|---------|---------|---------|------|
| **GSE149614** | **1** | ~70K | 10 | ✅ v2 | ✅ | 1.4 GB | ✅ 필수 |
| **GSE140228** | **1** | ~73K (droplet 66K) | 16+ | ✅ v2 | ❌ CD45+ | 305 MB | ✅ 필수 (면역 subset) |
| **GSE156625** | **1** | ~212K | 다수 | ✅ v2/v3 | ✅ | 1.7 GB (HCC) | ✅ 필수 (h5ad 제공) |
| **GSE151530** | **2** | ~57K | 46 | ✅ v2 | ✅ | 292 MB | ✅ 권장 (HCC 필터 필요) |
| GSE125449 | 3 | 극소 | 19 | ✅ v2 | ✅ | 43 MB | ❌ 제외 (GSE151530에 포함) |
| CNP0000650 | 3 | ~20K | 19 | 확인필요 | ✅ | 확인필요 | ⚠️ 접근성 이슈 |

### 최종 Discovery Cohort 구성

**Tier 1 (필수, ~355K cells):**
1. GSE149614 (~70K) — 전체 TME, multi-site
2. GSE140228 droplet (~66K) — 면역세포 gold standard
3. GSE156625 HCC-only (~212K) — 최대 규모, h5ad 즉시 사용

**Tier 2 (권장, +57K):**
4. GSE151530 (~57K) — 환자 수 추가 (HCC만 필터)

→ **총 예상: ~70+ 환자, ~350,000–410,000 세포**

### Validation Cohort
- TCGA-LIHC (n=371) — prognostic + TIDE/ESTIMATE
- GSE14520 (n=488) — independent survival (microarray)
- ICB: GSE140901 (n=24, NanoString) + EGA 신청 병행

---

## 6. 데이터 저장 위치

모든 raw/large 데이터는 외장하드에 저장:
```
/mnt/e/no_exp_paper/data/
├── raw/
│   ├── GSE149614/      ← count.txt.gz (158 MB) + metadata
│   ├── GSE140228/      ← mtx + barcodes + cellinfo (305 MB)
│   ├── GSE156625/      ← HCCscanpyobj.h5ad.gz (800 MB) + mtx
│   └── GSE151530/      ← mtx (292 MB) + Info.txt
├── external/
│   ├── TCGA_LIHC/      ← RNA-seq + clinical
│   ├── GSE14520/       ← microarray CEL + clinical (1 GB)
│   └── ICB_cohorts/    ← GSE140901, GSE285963, EGA data
└── spatial/            ← Visium datasets
```
