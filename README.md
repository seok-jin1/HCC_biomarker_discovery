# HCC Communication-Aware Biomarker Discovery

> **논문 가제:** Communication-aware deep learning framework identifies T cell exhaustion–macrophage interaction signatures predictive of immunotherapy response in hepatocellular carcinoma

## 저장소 구조

이 프로젝트는 **두 개의 물리적 위치**를 사용합니다:
- **`/home/laugh/no_exp_paper/`** — 코드, 문서, 결과물 (WSL Ubuntu 로컬)
- **`/mnt/e/no_exp_paper/data/`** — 대용량 원본 데이터 (외장하드, 용량 절약)

```
/home/laugh/no_exp_paper/          ← 프로젝트 루트 (코드 + 문서)
├── README.md                      ← 이 파일
├── reference_papers.md            ← 참고 논문 목록
├── approaches.md                  ← 접근법 비교 (Approach 1 채택)
├── config.py                      ← 경로 설정 (DATA_ROOT 등)
│
├── code/
│   ├── preprocessing/             ← QC, filtering, normalization, batch correction
│   │   ├── 00_download.sh         ← 데이터 다운로드 스크립트 (재현용)
│   │   ├── 01_qc_filtering.py
│   │   ├── 02_integration.py
│   │   └── 03_annotation.py
│   ├── analysis/                  ← 생물학적 분석
│   │   ├── 01_tcell_exhaustion.py
│   │   ├── 02_tam_subtyping.py
│   │   ├── 03_cellchat.py
│   │   ├── 04_nichenet.py
│   │   └── 05_scenic.py
│   ├── models/                    ← DL 모델 (CACE architecture)
│   │   ├── cace_model.py
│   │   ├── train.py
│   │   ├── evaluate.py
│   │   └── interpret.py
│   ├── utils/                     ← 공통 유틸리티
│   │   ├── io.py
│   │   ├── plotting.py
│   │   └── metrics.py
│   └── notebooks/                 ← 탐색적 분석, 시각화
│       ├── EDA_*.ipynb
│       └── Figure_*.ipynb
│
├── results/
│   ├── figures/
│   │   ├── main/                  ← Fig1–Fig8 (PDF 벡터)
│   │   └── supplementary/
│   └── tables/
│
├── manuscript/
│   ├── main/
│   │   ├── manuscript.tex
│   │   └── bibliography.bib
│   ├── supplementary/
│   └── revision/
│
└── docs/
    ├── specs/                     ← 확정된 설계/계획
    └── notes/                     ← 조사 노트, 메모
        └── dataset_investigation.md


/mnt/e/no_exp_paper/data/          ← 외장하드 (대용량 데이터)
├── raw/                           ← GEO 원본 다운로드 (읽기 전용)
│   ├── GSE149614/                 ← Lu et al. ~1.4 GB
│   ├── GSE140228/                 ← Zhang et al. ~305 MB
│   ├── GSE156625/                 ← Sharma et al. ~1.7 GB (HCC only)
│   └── GSE151530/                 ← Ma et al. ~292 MB
├── processed/                     ← 전처리 완료 데이터
│   ├── integrated.h5ad
│   ├── tcell_subset.h5ad
│   ├── myeloid_subset.h5ad
│   └── communication/            ← L-R interaction matrices
├── external/                      ← 검증용 데이터
│   ├── TCGA_LIHC/
│   ├── GSE14520/                  ← ~1.0 GB microarray
│   ├── ICB_cohorts/               ← GSE140901, GSE285963, EGA
│   └── spatial/                   ← Visium datasets
└── model_checkpoints/             ← DL 모델 체크포인트
```

## 규칙

### 0. 경로 관리 (최우선)
- **대용량 데이터는 반드시 외장하드(`/mnt/e/no_exp_paper/data/`)에 저장**
- 코드에서 경로를 직접 하드코딩하지 않음. `config.py`에서 경로 상수를 정의하고 import
- `config.py` 예시:
  ```python
  from pathlib import Path

  # 외장하드 마운트 확인
  DATA_ROOT = Path("/mnt/e/no_exp_paper/data")
  RAW_DIR = DATA_ROOT / "raw"
  PROCESSED_DIR = DATA_ROOT / "processed"
  EXTERNAL_DIR = DATA_ROOT / "external"

  # 프로젝트 로컬
  PROJECT_ROOT = Path("/home/laugh/no_exp_paper")
  CODE_DIR = PROJECT_ROOT / "code"
  RESULTS_DIR = PROJECT_ROOT / "results"
  ```
- 분석 시작 전 외장하드 마운트 상태 확인: `ls /mnt/e/no_exp_paper/data/`

### 1. 데이터 관리
- `raw/`는 **읽기 전용**으로 취급. 원본 데이터를 절대 수정하지 않음
- 모든 전처리 결과는 `processed/`에 저장
- 데이터 다운로드 스크립트는 `code/preprocessing/00_download.sh`에 기록하여 재현 가능하게
- GEO 다운로드 시 `wget` 또는 `prefetch` 명령어를 스크립트에 기록

### 2. 코드 규칙
- 파일명은 **번호_설명.py** 형식 (실행 순서 명시)
- 각 스크립트는 **하나의 분석 단계**만 수행
- 재현성: random seed 고정, 패키지 버전은 `requirements.txt` 또는 `environment.yml`로 관리
- notebooks는 **탐색/시각화 전용** — 최종 분석 로직은 반드시 `.py` 스크립트로 정리

### 3. Figure 관리
- Main figure: `Fig{N}_{설명}.pdf` (벡터 포맷)
- Supplementary: `SFig{N}_{설명}.pdf`
- Figure 생성 코드는 `code/notebooks/Figure_*.ipynb` 또는 `code/analysis/` 내 해당 스크립트에 포함
- 하나의 figure를 생성하는 코드와 결과물이 추적 가능해야 함

### 4. Manuscript
- 원고 작성은 LaTeX 또는 Word 중 타겟 저널 요구사항에 맞춤
- Bibliography는 `.bib` 파일로 관리 (Zotero 등 연동 권장)
- revision 폴더는 리비전 시 point-by-point response 포함

### 5. 버전 관리
- 주요 milestone마다 git commit
- `.gitignore`에 포함: `*.h5ad`, `*.rds`, `*.h5`, `__pycache__/`, `.ipynb_checkpoints/`, `data/` (로컬 data 폴더도 제외)
- 외장하드 데이터는 git 관리 대상이 아님

### 6. 문서화
- 분석 결정사항 및 근거는 `docs/notes/`에 markdown으로 기록
- 확정된 계획/설계는 `docs/specs/`에 보관
- 데이터셋 변경, 분석 방법 변경 시 해당 문서 업데이트
