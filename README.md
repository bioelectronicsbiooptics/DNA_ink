# Printable DNA data storage inks for object-bound digital provenance

이 저장소는 위 논문의 **DNA 데이터 인코딩·디코딩 코드**를 담고 있다.

- 원시 시퀀싱 데이터: NCBI Sequence Read Archive, BioProject **PRJNA1522504** (13 runs)

---

## 구성

```
encoding/    디지털 파일 → 150-nt oligonucleotide
decoding/    NGS 리드 → 원본 파일 복원
```

---

## encoding/

| 파일 | 내용 |
|---|---|
| `DNA_encoding_annotated.ipynb` | 인코딩 원본 노트북 (주석본) |
| `DNA_encoding_full.py` | 10개 데이터셋 일괄 인코딩 |
| `DNA_encoding_with_seeds.py` | seed sweep 적용 버전 |
| `DNA_encoding_withoutSHA256.py` | **대조군** — SHA-256 마스킹 미적용 |
| `DNA_encoding_noSeed_withoutSHA256.py` | **대조군** — seed·마스킹 모두 미적용 |

### 인코딩 구조

```
입력 파일 → 이진 변환 → 124-bit 단위 분할
  → 14-bit index 부여 = 138-bit fragment
  → fragment 간 XOR 중복 생성
  → Reed–Solomon (8 parity symbols = 64 bit)
  → 18-bit seed의 SHA-256 해시에서 202-bit 마스크 → XOR
  → seed 부착 = 220-bit codeword
  → 2 bit씩 염기 변환 (G=00, C=10, A=01, T=11) = 110-nt payload
  → 양 끝 20-nt primer = 150-nt oligonucleotide
```

**설계 제약** (110-nt payload 기준): homopolymer ≤ 2, GC 40–60 %, ΔG > −30 kcal mol⁻¹

제약을 만족하는 서열을 얻기 위해 0–262,143 범위의 18-bit seed를 무작위 배열한 뒤 순차 적용하여 세 조건을 처음으로 모두 만족하는 seed를 채택한다(seed sweep).

대조군 두 종은 이 마스킹의 효과를 비교하기 위한 것이다. 세 조건(Hash / NoSeed / WithoutSHA256)의 서열 특성 비교 결과는 논문 Figure 2에 수록되어 있다.

---

## decoding/

| 파일 | 내용 |
|---|---|
| `DNA_decoding_annotated.ipynb` | 디코딩 원본 노트북 (주석본) |
| `dna_analysis_v17.py` | NGS 리드 처리 및 복원 평가 |

### 리드 처리 흐름 (`dna_analysis_v17.py`)

```
FASTQ (forward + reverse) 병합
  → 길이 140 nt 미만 · 모호 염기(N) 포함 리드 제외
  → 앞 20 nt를 primer 및 그 역상보와 대조 (mismatch ≤ 3)
  → 역방향 리드를 정방향으로 변환
  → 양 끝 primer 제거 = 110-nt payload
  → payload 선두 index로 분류
  → index당 최빈 서열을 대표 서열로 선택
  → Reed–Solomon 복호 → 원본 복원
```

최소 리드 수 기준은 적용하지 않는다. 오류 집계는 primer 구간(1–20, 131–150 nt)과 payload 구간(21–130 nt)을 나누어 수행한다.

지원 데이터셋 10종: `Sheet`, `mp3`, `girl`, `game`, `S-K`, `S-E`, `S-I`, `H-K`, `H-C`, `H-E`

---

## 실행 환경

| 패키지 | 버전 |
|---|---|
| Python | 3.13.2 |
| NumPy | 2.2.6 |
| galois | 0.4.7 |
| reedsolo | 1.7.0 |
| ViennaRNA | 2.7.0 |

인코딩 단계의 서열 선별에는 ViennaRNA 2.5.0의 `RNAfold` 실행 파일을 사용하였다.

---

## 라이선스

[미정 — LICENSE 파일 참조]
