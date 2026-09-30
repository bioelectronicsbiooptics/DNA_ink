# -*- coding: utf-8 -*-
"""
================================================================================
DNA Encoding with Pre-defined Seeds
================================================================================

Original 인코딩 결과와 동일한 DNA 서열을 생성하기 위해
미리 추출된 시드(seed) 정보를 사용하는 인코딩 시스템

================================================================================
핵심 변경사항 (DNA_encoding_full.py 대비):
================================================================================

1. check_restrict() 함수 수정:
   - 기존: random.seed(0)으로 셔플된 시드 배열에서 제약조건 만족하는 첫 번째 시드 선택
   - 변경: Seed_Information_Report.json에서 로드한 시드를 직접 사용

2. 시드 적용 방식:
   - 각 DNA 서열 인덱스에 대응하는 시드 값을 JSON에서 가져옴
   - 해당 시드로 SHA-256 해시 → XOR 변환 수행
   - 제약조건 검사 없이 바로 시드 적용 (이미 검증된 시드이므로)

================================================================================
DNA 서열 구조 (150 nt):
================================================================================

┌─────────────────────────────────────────────────────────────────────────────┐
│  Forward Primer (20nt)  │     Payload (110nt)     │  Reverse Primer (20nt)  │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                    ┌───────────────┴───────────────┐
                    │ Seed (9nt) │ Data+RS (101nt)  │
                    └───────────────────────────────┘

================================================================================
시드(Seed) → XOR 변환 과정:
================================================================================

1. seed 값 (0 ~ 262143) 입력
2. SHA-256 해시 계산: hashlib.sha256(bytes(seed)).hexdigest()
3. 해시의 상위 55자리(16진수) → 220 bits 중 (220 - seed_len) bits 추출
4. 원본 데이터와 XOR 연산
5. 결과 앞에 seed의 binary 표현 (18 bits) 추가

================================================================================
"""

import numpy as np
import math
import hashlib
import time
import os
import json
from datetime import datetime

# 필요한 라이브러리 임포트
try:
    import galois
    from reedsolo import RSCodec
except ImportError:
    print("필수 라이브러리가 설치되지 않았습니다.")
    print("실행: pip install galois reedsolo")
    exit()


# =============================================================================
# 설정
# =============================================================================
# 출력 폴더: Encoded_output_ 접두사 + 프로그램 구분 + 날짜시간
OUTPUT_DIR = f"Encoded_output_with_seeds_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
SEED_JSON_PATH = "Seed_Information_Report.json"  # 시드 정보 JSON 파일
DATA_DIR = "DNAPrint_Py/data"  # 입력 데이터 폴더

# 파일명 → Original 파일명 매핑
FILE_MAPPING = {
    "Sonnet18_eng.txt": "Sonnet18_eng.txt",
    "Sonnet18_kor.txt": "Sonnet18_kor.txt",
    "Sonnet18_ital.txt": "Sonnet18_Ital.txt",  # 대소문자 차이 주의
    "hayeoga_kor.txt": "hayeoga_kor.txt",
    "hayeoga_eng.txt": "hayeoga_eng.txt",
    "hayeoga_chn.txt": "hayeoga_chn.txt",
    "Mozart_musicsheet.jpg": "Mozart_musicsheet_jpg.txt",
    "Mozart_mp3.mp3": "Mozart_mp3.txt",
    "girlwithballoon.jpg": "girlwithballoon.txt",
    "gamechanger.jpg": "gamechanger.txt",
}


# =============================================================================
# 1. 시드 정보 로더
# =============================================================================
class SeedLoader:
    """
    Seed_Information_Report.json에서 시드 정보를 로드하는 클래스

    JSON 구조:
    {
        "파일명.txt": {
            "primer_num": 6,
            "primer_desc": "Sonnet18 영어",
            "total_sequences": 66,
            "seeds": [
                {"index": 0, "decimal": 162683, "dna": "CATCTATCT"},
                {"index": 1, "decimal": 255138, "dna": "TTCAGCCGC"},
                ...
            ]
        }
    }
    """

    def __init__(self, json_path):
        """
        Args:
            json_path: Seed_Information_Report.json 파일 경로
        """
        self.json_path = json_path
        self.seed_data = None
        self._load()

    def _load(self):
        """JSON 파일 로드"""
        if not os.path.exists(self.json_path):
            raise FileNotFoundError(f"시드 파일을 찾을 수 없습니다: {self.json_path}")

        with open(self.json_path, 'r', encoding='utf-8') as f:
            self.seed_data = json.load(f)

        print(f"[SeedLoader] 시드 정보 로드 완료: {len(self.seed_data)}개 파일")

    def get_seeds_for_file(self, original_filename):
        """
        특정 파일의 시드 목록 반환

        Args:
            original_filename: Original encoded_DNA 폴더의 파일명

        Returns:
            list: [seed_decimal_0, seed_decimal_1, ...] 형태의 시드 리스트
        """
        if original_filename not in self.seed_data:
            raise KeyError(f"시드 정보가 없습니다: {original_filename}")

        file_info = self.seed_data[original_filename]
        seeds = [s['decimal'] for s in file_info['seeds']]

        print(f"[SeedLoader] {original_filename}: {len(seeds)}개 시드 로드")
        return seeds

    def get_primer_num(self, original_filename):
        """파일의 프라이머 번호 반환"""
        if original_filename not in self.seed_data:
            raise KeyError(f"프라이머 정보가 없습니다: {original_filename}")
        return self.seed_data[original_filename]['primer_num']


# =============================================================================
# 2. 파일 처리 클래스
# =============================================================================
class File_tobinary:
    """파일을 binary로 읽어오는 클래스"""

    @staticmethod
    def file_open(path):
        """
        파일을 binary 배열로 변환

        과정:
        1. 파일을 bytes로 읽음
        2. np.uint8 배열로 변환 (0~255)
        3. np.unpackbits로 이진수 배열로 변환

        Args:
            path: 파일 경로

        Returns:
            np.array: 이진수 배열 [0,1,1,0,1,0,...]
        """
        with open(path, 'rb') as f:
            file = f.read()
        file_np = np.frombuffer(file, dtype=np.uint8)
        unpacked = np.unpackbits(file_np)
        return unpacked


# =============================================================================
# 3. 타입 변환 클래스
# =============================================================================
class Typetransform:
    """데이터 타입 변환 함수 모음"""

    @staticmethod
    def deci_to_bi_null(array):
        """
        십진수 배열을 이진수 2D 배열로 변환 (최대값 기준 패딩)

        예: [3, 5, 7] → [[011], [101], [111]]
        """
        array = array.flatten()
        arr_max = array.max()
        arr_len = len(format(arr_max, 'b'))
        binary_arr = []
        for i in range(len(array)):
            tmp = np.fromiter(f'{array[i]:0{arr_len}b}', dtype=int)
            binary_arr.append(tmp)
        return np.array(binary_arr)

    @staticmethod
    def deci2bi(n, length=None):
        """
        십진수를 지정 길이의 이진 배열로 변환

        예: deci2bi(5, 8) → [0,0,0,0,0,1,0,1]
        """
        binary = bin(n)[2:]
        if length is not None:
            binary = '0' * (length - len(binary)) + binary
        return np.fromiter(binary, dtype=int)

    @staticmethod
    def transform(arr):
        """
        2D binary 배열을 DNA 서열로 변환

        변환 규칙:
            00 → G
            10 → C
            01 → A
            11 → T
        """
        bases = 'GCAT'
        result_DNA = np.zeros((np.size(arr, axis=0), np.size(arr, axis=1) // 2), dtype='<U3')
        for i in range(np.size(result_DNA, axis=0)):
            for j in range(np.size(result_DNA, axis=1)):
                DNA_index = arr[i, j*2] + 2 * arr[i, j*2+1]
                result_DNA[i, j] = bases[DNA_index]
        return result_DNA

    @staticmethod
    def transform_linear(arr):
        """1D binary 배열을 DNA 서열로 변환"""
        bases = 'GCAT'
        result_DNA = np.zeros((np.size(arr) // 2), dtype='<U3')
        for j in range(np.size(result_DNA)):
            DNA_index = arr[j*2] + 2 * arr[j*2+1]
            result_DNA[j] = bases[DNA_index]
        return result_DNA

    @staticmethod
    def DNA2bi(arr):
        """DNA 서열을 binary 배열로 변환"""
        bi_arr = np.zeros((len(arr) * 2), dtype=np.int64)
        for j in range(len(arr)):
            if arr[j] == 'G':
                bi_arr[2*j], bi_arr[2*j+1] = 0, 0
            elif arr[j] == 'C':
                bi_arr[2*j], bi_arr[2*j+1] = 1, 0
            elif arr[j] == 'A':
                bi_arr[2*j], bi_arr[2*j+1] = 0, 1
            elif arr[j] == 'T':
                bi_arr[2*j], bi_arr[2*j+1] = 1, 1
            else:
                raise ValueError(f"잘못된 염기: {arr[j]}")
        return bi_arr


# =============================================================================
# 4. 시드 변환 클래스 (핵심!)
# =============================================================================
class SeedTransform:
    """
    시드를 사용한 XOR 변환 클래스

    이 클래스가 Original과 동일한 결과를 만드는 핵심입니다.
    """

    @staticmethod
    def apply_seed(arr, seed, seed_len):
        """
        시드를 사용하여 데이터에 XOR 변환 적용

        과정:
        1. seed를 bytes로 변환 후 SHA-256 해시 계산
        2. 해시의 상위 55자리(16진수)를 binary로 변환
        3. (220 - seed_len) bits만 사용
        4. 원본 데이터와 XOR

        Args:
            arr: 원본 데이터 배열 (길이: 220 - seed_len)
            seed: 시드 값 (정수)
            seed_len: 시드 비트 길이 (기본 18)

        Returns:
            np.array: XOR 변환된 배열
        """
        # SHA-256 해시 계산
        hash_hex = hashlib.sha256(bytes(seed)).hexdigest()

        # 상위 55자리를 binary로 변환 (최대 220 bits)
        bi_hash = bin(int(hash_hex[:55], 16))[2:]

        # 길이 조정
        target_len = 220 - seed_len
        if len(bi_hash) < target_len:
            bi_hash = '0' * (target_len - len(bi_hash)) + bi_hash
        elif len(bi_hash) > target_len:
            bi_hash = bi_hash[:target_len]

        # numpy 배열로 변환
        bi_hash = np.fromiter(bi_hash, dtype=int)

        # XOR 연산
        xor_arr = np.logical_xor(arr, bi_hash).astype(int)

        return xor_arr

    @staticmethod
    def create_payload_with_seed(data, seed, seed_len):
        """
        데이터에 시드를 적용하여 220 bits 페이로드 생성

        결과 구조: [seed_binary (seed_len bits)] + [XOR된 데이터 (220-seed_len bits)]

        Args:
            data: 원본 데이터 (220 - seed_len bits)
            seed: 시드 값
            seed_len: 시드 비트 길이

        Returns:
            np.array: 220 bits 페이로드
        """
        # XOR 변환 적용
        xor_data = SeedTransform.apply_seed(data, seed, seed_len)

        # 220 bits 페이로드 생성
        payload = np.zeros(220, dtype=np.int64)
        payload[:seed_len] = Typetransform.deci2bi(seed, seed_len)  # 시드 binary
        payload[seed_len:] = xor_data  # XOR된 데이터

        return payload


# =============================================================================
# 5. 공통 유틸리티 클래스
# =============================================================================
class Common:
    """인코딩 공통 함수"""

    @staticmethod
    def arr_pad(arr, message_len):
        """데이터를 message_len 단위로 패딩하여 reshape"""
        arr_flat = np.ravel(arr)
        pad_cnt = math.ceil(arr.size / message_len) * message_len - arr.size

        for i in range(pad_cnt):
            arr_flat = np.pad(arr_flat, (0, 1), constant_values=0)

        np_arr = arr_flat.reshape(int(arr_flat.size / message_len), message_len)
        return np_arr, pad_cnt

    @staticmethod
    def xor_merge(arr):
        """XOR 패리티 데이터 생성 및 병합"""
        bool_pad = 0
        if np.size(arr, axis=0) % 2 == 1:
            arr = np.pad(arr, ((0, 1), (0, 0)), constant_values=0)
            bool_pad = 1

        xor_frame1 = arr[:np.size(arr, axis=0) // 2]
        xor_frame2 = arr[np.size(arr, axis=0) // 2:]
        xor_data = np.logical_xor(xor_frame1, xor_frame2).astype(int)

        arr_xor = np.concatenate((arr, xor_data), axis=0)
        return arr_xor, bool_pad

    @staticmethod
    def Reed(arr, bps):
        """Reed-Solomon 심볼 부착"""
        parity = 8
        gf_en = galois.GF(2**3)
        gf_en = np.array(gf_en(arr))

        rsc = RSCodec(parity)
        code_list = []
        for i in range(np.size(gf_en, axis=0)):
            temp = rsc.encode(gf_en[i])
            code_list.append(temp)
        arr_index_pari = np.array(code_list, dtype=np.int64)

        gf_parity = arr_index_pari[:, -parity:]
        gf_payload = arr_index_pari[:, :-parity]
        biparity = Typetransform.deci_to_bi_null(gf_parity)
        bi_parity = biparity.reshape(np.size(gf_en, axis=0), parity * 8)
        result_bi = np.concatenate((gf_payload, bi_parity), axis=1)
        return result_bi

    @staticmethod
    def primerSelect(select):
        """프라이머 선택"""
        primers = {
            1: ("CTGTCCATAGCCTTGTTCGT", "GCGGAAACGTAGTGAAGGTA"),   # Music Sheet
            2: ("CAAGTAACCGGCAACAACTG", "ACATAACAACCACCGCGAAA"),   # MP3
            3: ("TGTATTTCCTTCGGTGCTCC", "TTTCGACAACGGTCTGGTTT"),   # Girl with Balloon
            4: ("TCCTCAGCCGATGAAATTCC", "TGTACCATCCGTTTGACTGG"),   # Gamechanger
            5: ("AAGGCAAGTTGTTACCAGCA", "TGCGACCGTAATCAAACCAA"),   # Sonnet18 한글
            6: ("GAAGAGTTTAGCCACCTGGT", "AAGGCCAATTCGCGGTTATT"),   # Sonnet18 영어
            7: ("ATCCTGCAAACGCATTTCCT", "ATGCCTTTCCGAAGTTTCCA"),   # Sonnet18 이탈리아어
            8: ("AATCATGGCCTTCAAACCGT", "AACGCTCCGAAAGTCTTGTT"),   # Hayeoga 한글
            9: ("AATGGACGTTCCGCAATCAT", "AGAGCCGTGGCAATGTAAAT"),   # Hayeoga 한문
            10: ("GTCCAGGCAAAGATCCAGTT", "ACCACCGTTAGGCTAAAGTG"),  # Hayeoga 영어
            11: ("TAGCCTCCAGAATGAAACGG", "TTCAAGCCAAACCGTGTGTA"),  # NFT 서명
        }
        if select not in primers:
            raise ValueError(f"잘못된 프라이머 번호: {select}")
        return primers[select]

    @staticmethod
    def save_txtframe(arr, lang_type, data_pad, all_frags, seed_len, xor_pad, file_exp):
        """TXT 파일 헤더 저장"""
        file_type = Common.file_types(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4 + 4] = Typetransform.deci2bi(lang_type, 4)
        arr[0, 8:8 + 8] = Typetransform.deci2bi(data_pad, 8)
        arr[0, 16:16 + 20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 36:36 + 6] = Typetransform.deci2bi(seed_len, 6)
        arr[0, 42:42 + 2] = Typetransform.deci2bi(xor_pad, 2)
        return arr

    @staticmethod
    def save_MSmp3frame(arr, data_pad, all_frags, seed_len, xor_pad, file_exp):
        """MP3/JPG 파일 헤더 저장"""
        file_type = Common.file_types(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4 + 8] = Typetransform.deci2bi(data_pad, 8)
        arr[0, 12:12 + 20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 32:32 + 6] = Typetransform.deci2bi(seed_len, 6)
        arr[0, 38:38 + 2] = Typetransform.deci2bi(xor_pad, 2)
        return arr

    @staticmethod
    def save_jpgNFTframe(arr, jpg_pad, TXN_pad, jpg_frags, TXN_frags, all_frags, seed_len, xor_pad, file_exp):
        """NFT 파일 헤더 저장"""
        file_type = Common.file_types(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4 + 8] = Typetransform.deci2bi(jpg_pad, 8)
        arr[0, 12:12 + 8] = Typetransform.deci2bi(TXN_pad, 8)
        arr[0, 20:20 + 20] = Typetransform.deci2bi(jpg_frags, 20)
        arr[0, 40:40 + 20] = Typetransform.deci2bi(TXN_frags, 20)
        arr[0, 60:60 + 20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 80:80 + 6] = Typetransform.deci2bi(seed_len, 6)
        arr[0, 86:86 + 2] = Typetransform.deci2bi(xor_pad, 2)
        return arr

    @staticmethod
    def file_types(exp):
        """파일 타입 코드"""
        types = {"txt": 0, "jpg_TXN": 1, "mp3": 2, "jpg": 3}
        if exp not in types:
            raise ValueError(f"지원하지 않는 파일 유형: {exp}")
        return types[exp]


# =============================================================================
# 6. 조각 생성 클래스
# =============================================================================
class To_frags:
    """데이터를 DNA 조각으로 구성"""

    @staticmethod
    def file_tofrag_xor(arr, pay_len, index_len):
        """파일 1개의 XOR 조각 생성"""
        mess_len = pay_len - index_len

        arr_np, arr_pad_count = Common.arr_pad(arr, mess_len)
        arr_np = np.pad(arr_np, ((1, 0), (0, 0)), constant_values=0)

        xor_merged, boolxor_pad = Common.xor_merge(arr_np)

        index = np.arange(2**index_len)
        index_bi = Typetransform.deci_to_bi_null(index)

        index_arr = np.concatenate((index_bi[:np.size(xor_merged, axis=0)], xor_merged), axis=1)
        frag_count = np.size(index_arr, axis=0)

        return index_arr, arr_pad_count, frag_count, pay_len, index_len, boolxor_pad

    @staticmethod
    def file2_xor_tofrag(arr, arr2, pay_len, index_len):
        """파일 2개의 XOR 조각 생성"""
        mess_len = pay_len - index_len

        arr_np, arr_pad_count = Common.arr_pad(arr, mess_len)
        arr2_np, arr2_pad_count = Common.arr_pad(arr2, mess_len)

        merged_2 = np.concatenate((arr_np, arr2_np), axis=0)
        merged_2 = np.pad(merged_2, ((1, 0), (0, 0)), constant_values=0)

        xor_merged, boolxor_pad = Common.xor_merge(merged_2)

        index = np.arange(2**index_len)
        index_bi = Typetransform.deci_to_bi_null(index)

        xor_merged = np.concatenate((index_bi[:np.size(xor_merged, axis=0)], xor_merged), axis=1)
        arr_frag = np.size(arr_np, axis=0)
        arr2_frag = np.size(arr2_np, axis=0)
        frag_count = np.size(xor_merged, axis=0)

        return xor_merged, arr_pad_count, arr2_pad_count, arr_frag, arr2_frag, frag_count, pay_len, index_len, boolxor_pad


# =============================================================================
# 7. DNA 인코딩 클래스 (시드 기반 - 핵심!)
# =============================================================================
class DNA_Encoder_WithSeeds:
    """
    미리 정의된 시드를 사용하는 DNA 인코더

    기존 check_restrict() 함수가 제약조건을 만족하는 시드를 찾는 반면,
    이 클래스는 이미 알려진 시드를 직접 적용합니다.
    """

    @staticmethod
    def apply_seeds_to_data(arr, seeds, seed_len, primer_type):
        """
        시드 리스트를 사용하여 데이터를 DNA 서열로 변환

        핵심 로직:
        1. Reed-Solomon 심볼 부착
        2. 각 행에 대해 미리 정의된 시드 적용
        3. Binary → DNA 변환
        4. 프라이머 부착

        Args:
            arr: 입력 데이터 배열 (행: 서열 수, 열: 데이터 bits)
            seeds: 시드 리스트 [seed_0, seed_1, ...]
            seed_len: 시드 비트 길이
            primer_type: 프라이머 번호

        Returns:
            list: DNA 서열 리스트
        """
        # 1. Reed-Solomon 심볼 부착
        result_bi = Common.Reed(arr, 156 - seed_len)

        # 2. 프라이머 가져오기
        F_primer, R_primer = Common.primerSelect(primer_type)

        # 3. 각 서열에 시드 적용
        num_seqs = np.size(result_bi, axis=0)
        wo_homo = np.zeros((num_seqs, 220), dtype=np.int64)

        print(f"[DNA_Encoder] {num_seqs}개 서열에 시드 적용 중...")

        for i in range(num_seqs):
            if i >= len(seeds):
                raise ValueError(f"시드 부족: 인덱스 {i}, 시드 수 {len(seeds)}")

            seed = seeds[i]

            # 시드 적용하여 220 bits 페이로드 생성
            payload = SeedTransform.create_payload_with_seed(
                result_bi[i], seed, seed_len
            )
            wo_homo[i] = payload

            if i % 100 == 0 and i > 0:
                print(f"  ... {i}/{num_seqs} 완료")

        # 4. Binary → DNA 변환
        DNA_seq = Typetransform.transform(wo_homo)

        # 5. 프라이머 부착
        DNA_library = []
        for i in range(len(DNA_seq)):
            tmp_DNA = "".join(DNA_seq[i, :].tolist())
            full_seq = F_primer + tmp_DNA + R_primer
            DNA_library.append(full_seq)

        print(f"[DNA_Encoder] 완료! {len(DNA_library)}개 DNA 서열 생성")
        return DNA_library


# =============================================================================
# 8. 언어 타입 변환
# =============================================================================
def lang_num(language):
    """언어 코드 반환"""
    lang_map = {
        'ENGLISH': 0,
        'KOREAN': 1,
        'ITALIANO': 2,
        'CHINESE': 3
    }
    lang_upper = language.upper()
    if lang_upper not in lang_map:
        raise ValueError(f"지원하지 않는 언어: {language}")
    return lang_map[lang_upper]


# =============================================================================
# 9. 메인 인코딩 클래스
# =============================================================================
class MainEncoder:
    """시드 기반 메인 인코딩 함수"""

    def __init__(self, seed_loader):
        """
        Args:
            seed_loader: SeedLoader 인스턴스
        """
        self.seed_loader = seed_loader

    def txt_encoding(self, path, primer_type, seed_len, lang_type, index_len, original_filename):
        """
        TXT 파일 인코딩

        Args:
            path: 입력 파일 경로
            primer_type: 프라이머 번호
            seed_len: 시드 비트 길이
            lang_type: 언어 코드
            index_len: 인덱스 비트 길이
            original_filename: Original 파일명 (시드 로드용)

        Returns:
            list: DNA 서열 리스트
        """
        print(f"\n{'='*60}")
        print(f"TXT 인코딩: {path}")
        print(f"Original: {original_filename}")
        print(f"{'='*60}")

        # 시드 로드
        seeds = self.seed_loader.get_seeds_for_file(original_filename)

        # 파일 읽기
        txt_np = File_tobinary.file_open(path)

        # 조각 생성
        xor_data, pad_count, all_frags, pay_len, index_bits_len, xor_pad = To_frags.file_tofrag_xor(
            txt_np, 156 - seed_len, index_len
        )

        arr_frags = all_frags // 3 * 2
        only_data = xor_data[:arr_frags, index_bits_len:]

        # 헤더 저장
        arr_data = Common.save_txtframe(only_data, lang_type, pad_count, all_frags, seed_len, xor_pad, 'txt')
        re_xor, rexor_pad = Common.xor_merge(arr_data)
        xor_data[:, index_bits_len:] = re_xor

        # DNA 인코딩 (시드 사용)
        DNA_seqs = DNA_Encoder_WithSeeds.apply_seeds_to_data(
            xor_data, seeds, seed_len, primer_type
        )

        # 저장
        name = os.path.splitext(os.path.basename(path))[0]
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_path = os.path.join(OUTPUT_DIR, f"{name}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in DNA_seqs:
                f.write(seq + '\n')

        print(f"저장 완료: {save_path}")
        return DNA_seqs

    def MS_sheet_onefile(self, path, primer_type, seed_len, index_len, exp, original_filename):
        """JPG/MP3 파일 인코딩"""
        print(f"\n{'='*60}")
        print(f"{exp.upper()} 인코딩: {path}")
        print(f"Original: {original_filename}")
        print(f"{'='*60}")

        # 시드 로드
        seeds = self.seed_loader.get_seeds_for_file(original_filename)

        # 파일 읽기
        file_np = File_tobinary.file_open(path)

        # 조각 생성
        file_xor_data, file_pad_count, file_xor_frags, file_pay_len, file_index_bits_len, file_xor_pad = To_frags.file_tofrag_xor(
            file_np, 156 - seed_len, index_len
        )

        arr_frags = file_xor_frags // 3 * 2
        file_only_data = file_xor_data[:arr_frags, file_index_bits_len:]

        # 헤더 저장
        file_arr_data = Common.save_MSmp3frame(file_only_data, file_pad_count, file_xor_frags, seed_len, file_xor_pad, exp)
        re_xor, rexor_pad = Common.xor_merge(file_arr_data)
        file_xor_data[:, file_index_bits_len:] = re_xor

        # DNA 인코딩 (시드 사용)
        DNA_seqs = DNA_Encoder_WithSeeds.apply_seeds_to_data(
            file_xor_data, seeds, seed_len, primer_type
        )

        # 저장
        name = os.path.splitext(os.path.basename(path))[0]
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_path = os.path.join(OUTPUT_DIR, f"{name}_{exp}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in DNA_seqs:
                f.write(seq + '\n')

        print(f"저장 완료: {save_path}")
        return DNA_seqs


# =============================================================================
# 10. 메인 실행
# =============================================================================
if __name__ == "__main__":
    print("=" * 80)
    print("DNA Encoding with Pre-defined Seeds")
    print("Original과 동일한 DNA 서열 생성")
    print("=" * 80)

    # 설정
    seed_len = 18
    index_len = 14
    data_dir = DATA_DIR

    # 시드 로더 초기화
    try:
        seed_loader = SeedLoader(SEED_JSON_PATH)
    except FileNotFoundError as e:
        print(f"오류: {e}")
        print("먼저 extract_seeds.py를 실행하여 시드 정보를 추출하세요.")
        exit(1)

    # 인코더 초기화
    encoder = MainEncoder(seed_loader)

    # 인코딩할 파일 목록
    txt_files = {
        'hayeoga_kor.txt': (8, 'korean', 'hayeoga_kor.txt'),
        'hayeoga_eng.txt': (10, 'english', 'hayeoga_eng.txt'),
        'hayeoga_chn.txt': (9, 'chinese', 'hayeoga_chn.txt'),
        'Sonnet18_kor.txt': (5, 'korean', 'Sonnet18_kor.txt'),
        'Sonnet18_eng.txt': (6, 'english', 'Sonnet18_eng.txt'),
        'Sonnet18_ital.txt': (7, 'italiano', 'Sonnet18_Ital.txt'),
    }

    # TXT 파일 인코딩
    for filename, (primer, lang, orig_name) in txt_files.items():
        filepath = os.path.join(data_dir, filename)
        if os.path.exists(filepath):
            lang_type = lang_num(lang)
            encoder.txt_encoding(filepath, primer, seed_len, lang_type, index_len, orig_name)
        else:
            print(f"파일 없음: {filepath}")

    print("\n" + "=" * 80)
    print("인코딩 완료!")
    print(f"출력 폴더: {OUTPUT_DIR}")
    print("=" * 80)
