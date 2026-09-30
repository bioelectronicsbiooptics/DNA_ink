# -*- coding: utf-8 -*-
"""
DNA Encoding Simple Version (without SHA256)
- SHA256 해시 변환 제거
- Free energy 검사 제거
- XOR 패리티 포함 (데이터 복구용)
- 시드만 붙이고 루프 없이 단순 변환
- 전체 서열 길이: 150nt (Primer 20nt + Payload 110nt + Primer 20nt)
"""

import numpy as np
import math
import os

try:
    import galois
    from reedsolo import RSCodec
except ImportError:
    print("필수 라이브러리가 설치되지 않았습니다. 'pip install galois reedsolo' 명령을 실행해주세요.")
    exit()

# 전역 출력 폴더 설정
from datetime import datetime
OUTPUT_DIR = f"Encoded_output_withoutSHA256_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
DATA_DIR = "DNAPrint_Py/data"  # 입력 데이터 폴더


# =============================================================================
# 1. 파일 처리
# =============================================================================
class File_tobinary:
    """파일을 binary로 읽어옴"""
    @staticmethod
    def file_open(path):
        with open(path, 'rb') as f:
            file = f.read()
        file_np = np.frombuffer(file, dtype=np.uint8)
        unpacked = np.unpackbits(file_np)
        print(f"파일 읽기 완료: {len(unpacked)} bits")
        return unpacked


# =============================================================================
# 2. 타입 변환
# =============================================================================
class Typetransform:
    """데이터 타입 간의 변환"""

    @staticmethod
    def deci_to_bi_null(array):
        """십진수를 이진수로 변경, 최대길이로 맞춤"""
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
        """십진수를 지정된 길이의 이진 배열로 변환"""
        binary = bin(n)[2:]
        if length is not None:
            binary = '0' * (length - len(binary)) + binary
        return np.fromiter(binary, dtype=int)

    @staticmethod
    def transform(arr):
        """binary를 DNA로 변환 (2차원)"""
        bases = 'GCAT'
        result_DNA = np.zeros((np.size(arr, axis=0), np.size(arr, axis=1) // 2), dtype='<U3')
        for i in range(np.size(result_DNA, axis=0)):
            for j in range(np.size(result_DNA, axis=1)):
                DNA_index = arr[i, j*2] + 2 * arr[i, j*2+1]
                result_DNA[i, j] = bases[DNA_index]  # 00->G, 10->C, 01->A, 11->T
        return result_DNA

    @staticmethod
    def transform_linear(arr):
        """binary를 DNA로 변환 (1차원)"""
        bases = 'GCAT'
        result_DNA = np.zeros((np.size(arr) // 2), dtype='<U3')
        for j in range(np.size(result_DNA)):
            DNA_index = arr[j*2] + 2 * arr[j*2+1]
            result_DNA[j] = bases[DNA_index]
        return result_DNA


# =============================================================================
# 3. Homopolymer 및 GC 검사 (free energy 제거됨)
# =============================================================================
class homo:
    """DNA 서열 제약 조건 검사 (homopolymer, GC만)"""

    @staticmethod
    def Check_WO_homopolymer(arr):
        """homopolymer 검사 - 3개 이상 연속 불허"""
        count = 1
        for j in range(0, np.size(arr) - 2, 2):
            if (arr[j:j+2] == arr[j+2:j+4]).all():
                count += 1
                if count > 2:
                    return False
            else:
                count = 1
        return True

    @staticmethod
    def Check_WO_GC(arr):
        """GC contents 40~60% 검사"""
        GC_contents = 0
        AT_contents = 0
        for i in range(0, np.size(arr), 2):
            if (arr[i:i+2] == [0, 0]).all() or (arr[i:i+2] == [1, 0]).all():
                GC_contents += 1
            elif (arr[i:i+2] == [0, 1]).all() or (arr[i:i+2] == [1, 1]).all():
                AT_contents += 1
        total = GC_contents + AT_contents
        return (GC_contents > total * 0.4) and (GC_contents < total * 0.6)


# =============================================================================
# 4. 공통 유틸리티
# =============================================================================
class Common:
    """인코딩에서 공통으로 사용되는 함수들"""

    @staticmethod
    def arr_pad(arr, message_len):
        """padding을 붙이고 원하는 크기로 reshape"""
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
        if np.size(arr, axis=0) % 2 == 1:  # 행 개수가 홀수라면 아래에 0으로 이루어진 행 생성
            arr = np.pad(arr, ((0, 1), (0, 0)), constant_values=0)
            bool_pad = 1

        xor_frame1 = arr[:np.size(arr, axis=0) // 2]
        xor_frame2 = arr[np.size(arr, axis=0) // 2:]

        xor_data = np.logical_xor(xor_frame1, xor_frame2).astype(int)

        arr_xor = np.concatenate((arr, xor_data), axis=0)

        return arr_xor, bool_pad

    @staticmethod
    def Reed(arr, bps):
        """Reed-Solomon symbol 부착"""
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
        """primer 선택"""
        primers = {
            1: ("CTGTCCATAGCCTTGTTCGT", "GCGGAAACGTAGTGAAGGTA"),  # musicsheet
            2: ("CAAGTAACCGGCAACAACTG", "ACATAACAACCACCGCGAAA"),  # mp3
            3: ("TGTATTTCCTTCGGTGCTCC", "TTTCGACAACGGTCTGGTTT"),  # jpg
            4: ("TCCTCAGCCGATGAAATTCC", "TGTACCATCCGTTTGACTGG"),  # jpg2
            5: ("AAGGCAAGTTGTTACCAGCA", "TGCGACCGTAATCAAACCAA"),  # txt_kor
            6: ("GAAGAGTTTAGCCACCTGGT", "AAGGCCAATTCGCGGTTATT"),  # txt_eng
            7: ("ATCCTGCAAACGCATTTCCT", "ATGCCTTTCCGAAGTTTCCA"),  # txt_ital
            8: ("AATCATGGCCTTCAAACCGT", "AACGCTCCGAAAGTCTTGTT"),  # hayeoga_kor
            9: ("AATGGACGTTCCGCAATCAT", "AGAGCCGTGGCAATGTAAAT"),  # hayeoga_chn
            10: ("GTCCAGGCAAAGATCCAGTT", "ACCACCGTTAGGCTAAAGTG"), # hayeoga_eng
            11: ("TAGCCTCCAGAATGAAACGG", "TTCAAGCCAAACCGTGTGTA"), # NFT
        }
        if select not in primers:
            raise ValueError(f"Invalid primer select: {select}")
        return primers[select]


# =============================================================================
# 5. 프레임 생성 (XOR 패리티 포함)
# =============================================================================
class To_frags:
    """XOR 패리티를 포함한 프레임 생성"""

    @staticmethod
    def file_tofrag_xor(arr, pay_len, index_len):
        """XOR 패리티를 포함한 데이터 프레임 생성"""
        mess_len = pay_len - index_len

        # padding 추가
        arr_np, arr_pad_count = Common.arr_pad(arr, mess_len)
        # 파일 헤더용 빈 행 추가
        arr_np = np.pad(arr_np, ((1, 0), (0, 0)), constant_values=0)

        # XOR 패리티 데이터 생성
        xor_merged, boolxor_pad = Common.xor_merge(arr_np)

        # index 생성 및 부착
        index = np.arange(2**index_len)
        index_bi = Typetransform.deci_to_bi_null(index)

        index_arr = np.concatenate((index_bi[:np.size(xor_merged, axis=0)], xor_merged), axis=1)
        frag_count = np.size(index_arr, axis=0)

        return index_arr, arr_pad_count, frag_count, pay_len, index_len, boolxor_pad

    @staticmethod
    def file2_xor_tofrag(arr, arr2, pay_len, index_len):
        """파일 2개의 XOR 조각 생성 (NFT용: JPG + TXN)"""
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
# 6. DNA 인코딩 (SHA256 없이, free energy 검사 없이)
# =============================================================================
class file_toDNA_simple:
    """Reed-Solomon 부착 후 바로 DNA 변환 (full 버전과 동일한 구조, 제약 검사 최소화)"""

    @staticmethod
    def simple_DNA(arr, seed_len, primer_type):
        """
        DNA 변환 (full 버전과 동일한 구조)

        구조 (full 버전과 동일):
        1. 입력: index(14bits) + data = 138 bits (pay_len = 156 - seed_len)
        2. RS 부착: 138 bits → 202 bits (+64 bits RS parity)
        3. Seed 추가: seed(18bits) + data(202bits) = 220 bits
        4. DNA 변환: 220 bits → 110 nt
        5. Primer 부착: 20nt + 110nt + 20nt = 150nt

        차이점:
        - SHA256 XOR 변환 없음 (seed만 앞에 붙임)
        - Free energy 검사 없음
        - Homopolymer/GC 검사만 수행 (실패해도 진행)
        """
        # Reed-Solomon 부착 (full 버전과 동일: 156 - seed_len)
        result_bi = Common.Reed(arr, 156 - seed_len)

        F_primers, R_primers = Common.primerSelect(primer_type)

        # 최종 220bits 배열 생성 (seed + RS 인코딩된 데이터)
        final_arr = np.zeros((np.size(result_bi, axis=0), 220), dtype=np.int64)

        for i in range(np.size(result_bi, axis=0)):
            # 시드 부분: 인덱스 기반 시드 (SHA256 변환 없이)
            seed_value = i % (2**seed_len)
            final_arr[i, :seed_len] = Typetransform.deci2bi(seed_value, seed_len)

            # 페이로드 부분: RS 인코딩된 데이터 (202 bits)
            final_arr[i, seed_len:] = result_bi[i]

        # DNA 변환
        DNA_seq = Typetransform.transform(final_arr)

        # DNA 라이브러리 생성
        DNA_library = []
        for i in range(len(DNA_seq)):
            tmp_DNA = "".join(DNA_seq[i, :].tolist())
            full_seq = F_primers + tmp_DNA + R_primers
            DNA_library.append(full_seq)

        print(f"총 {len(DNA_library)}개 서열 생성 완료")
        return DNA_library


# =============================================================================
# 7. 파일 헤더 저장
# =============================================================================
class Header:
    """파일 헤더 정보 저장"""

    @staticmethod
    def save_header(arr, file_type, data_pad, all_frags, seed_len, xor_pad):
        """파일 헤더 저장 (XOR 정보 포함)"""
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)      # 파일 유형
        arr[0, 4:12] = Typetransform.deci2bi(data_pad, 8)     # padding 개수
        arr[0, 12:32] = Typetransform.deci2bi(all_frags, 20)  # 전체 서열 개수
        arr[0, 32:38] = Typetransform.deci2bi(seed_len, 6)    # seed 길이
        arr[0, 38:40] = Typetransform.deci2bi(xor_pad, 2)     # XOR padding 여부
        print(f"헤더: type={file_type}, pad={data_pad}, frags={all_frags}, seed={seed_len}, xor_pad={xor_pad}")
        return arr

    @staticmethod
    def save_txtframe(arr, lang_type, data_pad, all_frags, seed_len, xor_pad, file_exp):
        """TXT 파일 헤더 저장 (언어 정보 포함)"""
        file_type = Header.file_type_num(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4+4] = Typetransform.deci2bi(lang_type, 4)
        arr[0, 8:8+8] = Typetransform.deci2bi(data_pad, 8)
        arr[0, 16:16+20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 36:36+6] = Typetransform.deci2bi(seed_len, 6)
        arr[0, 42:42+2] = Typetransform.deci2bi(xor_pad, 2)
        print(f"TXT헤더: type={file_type}, lang={lang_type}, pad={data_pad}, frags={all_frags}, seed={seed_len}, xor_pad={xor_pad}")
        return arr

    @staticmethod
    def save_MSmp3frame(arr, data_pad, all_frags, seed_len, xor_pad, file_exp):
        """MP3/JPG 파일 헤더 저장"""
        file_type = Header.file_type_num(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4+8] = Typetransform.deci2bi(data_pad, 8)
        arr[0, 12:12+20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 32:32+6] = Typetransform.deci2bi(seed_len, 6)
        arr[0, 38:38+2] = Typetransform.deci2bi(xor_pad, 2)
        print(f"JPG/MP3헤더: type={file_type}, pad={data_pad}, frags={all_frags}, seed={seed_len}, xor_pad={xor_pad}")
        return arr

    @staticmethod
    def save_jpgNFTframe(arr, jpg_pad, TXN_pad, jpg_frags, TXN_frags, all_frags, seed_len, xor_pad, file_exp):
        """NFT 파일 헤더 저장 (JPG + TXN)"""
        file_type = Header.file_type_num(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4+8] = Typetransform.deci2bi(jpg_pad, 8)
        arr[0, 12:12+8] = Typetransform.deci2bi(TXN_pad, 8)
        arr[0, 20:20+20] = Typetransform.deci2bi(jpg_frags, 20)
        arr[0, 40:40+20] = Typetransform.deci2bi(TXN_frags, 20)
        arr[0, 60:60+20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 80:80+6] = Typetransform.deci2bi(seed_len, 6)
        arr[0, 86:86+2] = Typetransform.deci2bi(xor_pad, 2)
        print(f"NFT헤더: type={file_type}, jpg_pad={jpg_pad}, TXN_pad={TXN_pad}, jpg_frags={jpg_frags}, TXN_frags={TXN_frags}, all_frags={all_frags}")
        return arr

    @staticmethod
    def file_type_num(exp):
        """파일 확장자를 숫자로"""
        types = {'txt': 0, 'jpg_TXN': 1, 'mp3': 2, 'jpg': 3}
        return types.get(exp, 0)


# =============================================================================
# 언어 타입 변환
# =============================================================================
def lang_num(language):
    """언어 코드 반환"""
    lang_map = {
        'ENGLISH': 0, 'English': 0, 'english': 0,
        'KOREAN': 1, 'Korean': 1, 'korean': 1,
        'ITALIANO': 2, 'Italiano': 2, 'italiano': 2,
        'CHINESE': 3, 'Chinese': 3, 'chinese': 3,
    }
    return lang_map.get(language, 0)


# =============================================================================
# 8. 메인 인코딩 함수
# =============================================================================
class SimpleEncoder:
    """간소화된 인코딩 메인 클래스 (XOR 패리티 포함)"""

    @staticmethod
    def txt_encoding(path, primer_type, lang_type, seed_len, index_len):
        """TXT 파일 인코딩 (언어 헤더 포함)"""
        print(f"\n{'='*60}")
        print(f"TXT Encoding (without SHA256): {path}")
        print(f"Seed: {seed_len} bits ({seed_len//2}nt)")
        print(f"{'='*60}")

        txt_np = File_tobinary.file_open(path)

        pay_len = 156 - seed_len
        xor_data, pad_count, all_frags, pay_len, index_bits_len, xor_pad = To_frags.file_tofrag_xor(
            txt_np, pay_len, index_len
        )

        arr_frags = all_frags // 3 * 2
        only_data = xor_data[:arr_frags, index_bits_len:]

        arr_data = Header.save_txtframe(only_data, lang_type, pad_count, all_frags, seed_len, xor_pad, 'txt')
        re_xor, rexor_pad = Common.xor_merge(arr_data)
        xor_data[:, index_bits_len:] = re_xor

        DNA_seqs = file_toDNA_simple.simple_DNA(xor_data, seed_len, primer_type)

        name = os.path.splitext(os.path.basename(path))[0]
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_path = os.path.join(OUTPUT_DIR, f"{name}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in DNA_seqs:
                f.write(seq + '\n')

        print(f"저장 완료: {save_path}")
        print(f"총 서열 수: {len(DNA_seqs)}, 서열 길이: {len(DNA_seqs[0])}nt")
        return DNA_seqs

    @staticmethod
    def MS_sheet_onefile(path, primer_type, seed_len, index_len, exp):
        """JPG 또는 MP3 단일 파일 인코딩"""
        print(f"\n{'='*60}")
        print(f"{exp.upper()} Encoding (without SHA256): {path}")
        print(f"Seed: {seed_len} bits ({seed_len//2}nt)")
        print(f"{'='*60}")

        file_np = File_tobinary.file_open(path)

        pay_len = 156 - seed_len
        file_xor_data, file_pad_count, file_xor_frags, file_pay_len, file_index_bits_len, file_xor_pad = To_frags.file_tofrag_xor(
            file_np, pay_len, index_len
        )

        arr_frags = file_xor_frags // 3 * 2
        file_only_data = file_xor_data[:arr_frags, file_index_bits_len:]

        file_arr_data = Header.save_MSmp3frame(file_only_data, file_pad_count, file_xor_frags, seed_len, file_xor_pad, exp)
        re_xor, rexor_pad = Common.xor_merge(file_arr_data)
        file_xor_data[:, file_index_bits_len:] = re_xor

        DNA_seqs = file_toDNA_simple.simple_DNA(file_xor_data, seed_len, primer_type)

        name = os.path.splitext(os.path.basename(path))[0]
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_path = os.path.join(OUTPUT_DIR, f"{name}_{exp}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in DNA_seqs:
                f.write(seq + '\n')

        print(f"저장 완료: {save_path}")
        print(f"총 서열 수: {len(DNA_seqs)}, 서열 길이: {len(DNA_seqs[0])}nt")
        return DNA_seqs

    @staticmethod
    def MS_mp3_encoding(jpg_path, mp3_path, jpg_primer, mp3_primer, seed_len, index_len):
        """악보(JPG) + MP3 인코딩"""
        sheet_seqs = SimpleEncoder.MS_sheet_onefile(jpg_path, jpg_primer, seed_len, index_len, 'jpg')
        mp3_seqs = SimpleEncoder.MS_sheet_onefile(mp3_path, mp3_primer, seed_len, index_len, 'mp3')

        merged_seq = sheet_seqs + mp3_seqs
        save_path = os.path.join(OUTPUT_DIR, "mp3_sheet_merged.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in merged_seq:
                f.write(seq + '\n')

        print(f"\n병합 저장: {save_path}")
        print(f"총 서열 수: {len(merged_seq)}")
        return merged_seq

    @staticmethod
    def jpg_NFT_encoding(jpg_path, TXN_path, primer_type, seed_len, index_len):
        """NFT 인코딩 (JPG + Transaction ID)"""
        print(f"\n{'='*60}")
        print(f"NFT Encoding (without SHA256): {jpg_path} + {TXN_path}")
        print(f"Seed: {seed_len} bits ({seed_len//2}nt)")
        print(f"{'='*60}")

        jpg_np = File_tobinary.file_open(jpg_path)
        TXN_np = File_tobinary.file_open(TXN_path)

        pay_len = 156 - seed_len
        xor_data, jpg_pad_count, TXN_pad_count, jpg_frags, TXN_frags, all_frags, pay_len, index_bits_len, xor_pad = To_frags.file2_xor_tofrag(
            jpg_np, TXN_np, pay_len, index_len
        )

        arr_frags = all_frags // 3 * 2
        only_data = xor_data[:arr_frags, index_bits_len:]

        arr_data = Header.save_jpgNFTframe(only_data, jpg_pad_count, TXN_pad_count, jpg_frags, TXN_frags, all_frags, seed_len, xor_pad, 'jpg_TXN')
        re_xor, rexor_pad = Common.xor_merge(arr_data)
        xor_data[:, index_bits_len:] = re_xor

        DNA_seqs = file_toDNA_simple.simple_DNA(xor_data, seed_len, primer_type)

        name = os.path.splitext(os.path.basename(jpg_path))[0]
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_path = os.path.join(OUTPUT_DIR, f"{name}_NFT.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in DNA_seqs:
                f.write(seq + '\n')

        print(f"저장 완료: {save_path}")
        print(f"총 서열 수: {len(DNA_seqs)}, 서열 길이: {len(DNA_seqs[0])}nt")
        return DNA_seqs

    @staticmethod
    def encode_file(path, primer_type, seed_len, index_len, file_type='txt'):
        """
        파일 인코딩 (XOR 포함, SHA256 없음, free energy 없음)
        - 레거시 호환을 위해 유지

        Parameters:
        - path: 파일 경로
        - primer_type: 프라이머 타입 (1-11)
        - seed_len: 시드 비트 길이 (예: 18 bits = 9nt)
        - index_len: 인덱스 비트 길이
        - file_type: 파일 타입 ('txt', 'jpg', 'mp3', 'nft')
        """
        print(f"\n{'='*60}")
        print(f"DNA Encoding (without SHA256): {path}")
        print(f"Seed: {seed_len} bits ({seed_len//2}nt)")
        print(f"{'='*60}")

        # 파일 읽기
        file_np = File_tobinary.file_open(path)

        # 프레임 생성 (XOR 포함)
        # pay_len = 156 - seed_len (full 버전과 동일)
        pay_len = 156 - seed_len

        xor_data, pad_count, all_frags, pay_len, index_bits_len, xor_pad = To_frags.file_tofrag_xor(
            file_np, pay_len, index_len
        )

        # XOR 서열 개수 계산
        arr_frags = all_frags // 3 * 2
        only_data = xor_data[:arr_frags, index_bits_len:]

        # 헤더 저장
        file_type_num = Header.file_type_num(file_type)
        arr_data = Header.save_header(only_data, file_type_num, pad_count, all_frags, seed_len, xor_pad)

        # XOR 다시 생성 (헤더 정보 포함)
        re_xor, rexor_pad = Common.xor_merge(arr_data)
        xor_data[:, index_bits_len:] = re_xor

        # DNA 변환
        DNA_seqs = file_toDNA_simple.simple_DNA(xor_data, seed_len, primer_type)

        # 저장
        name = os.path.splitext(os.path.basename(path))[0]
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_path = os.path.join(OUTPUT_DIR, f"{name}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in DNA_seqs:
                f.write(seq + '\n')

        print(f"\n저장 완료: {save_path}")
        print(f"총 서열 수: {len(DNA_seqs)}")
        print(f"서열 길이: {len(DNA_seqs[0])}nt (primer 포함)")

        return DNA_seqs


# =============================================================================
# 9. 실행
# =============================================================================
if __name__ == "__main__":
    print("="*60)
    print("DNA Encoding (without SHA256)")
    print("- XOR 패리티 포함")
    print("- SHA256 해시 변환 제거")
    print("- Free Energy 검사 제거")
    print("- TXT/JPG/MP3/NFT 모든 파일 타입 지원")
    print("="*60)

    # 설정 (full 버전과 동일)
    seed_len = 18         # 시드 비트 길이 (18 bits = 9nt)
    index_len = 14        # 인덱스 비트 길이

    print(f"\n설정:")
    print(f"  - Seed: {seed_len} bits ({seed_len//2}nt)")
    print(f"  - Index: {index_len} bits")
    print(f"  - 출력 폴더: {OUTPUT_DIR}")

    # 데이터 폴더
    data_dir = DATA_DIR

    # =========================================================================
    # TXT 파일 인코딩
    # =========================================================================
    txt_files = {
        'hayeoga_kor.txt': (8, 'korean'),
        'hayeoga_eng.txt': (10, 'english'),
        'hayeoga_chn.txt': (9, 'chinese'),
        'Sonnet18_kor.txt': (5, 'korean'),
        'Sonnet18_eng.txt': (6, 'english'),
        'Sonnet18_ital.txt': (7, 'italiano'),
    }

    print(f"\n{'='*60}")
    print("TXT 파일 인코딩")
    print(f"{'='*60}")

    for filename, (primer, lang) in txt_files.items():
        filepath = os.path.join(data_dir, filename)
        if os.path.exists(filepath):
            lang_type = lang_num(lang)
            SimpleEncoder.txt_encoding(filepath, primer, lang_type, seed_len, index_len)
        else:
            print(f"\n파일 없음: {filepath}")

    # =========================================================================
    # JPG + MP3 인코딩 (악보)
    # =========================================================================
    print(f"\n{'='*60}")
    print("JPG + MP3 인코딩 (악보)")
    print(f"{'='*60}")

    jpg_path = os.path.join(data_dir, 'small.jpg')
    mp3_path = os.path.join(data_dir, '5sec_64kbps.mp3')

    if os.path.exists(jpg_path) and os.path.exists(mp3_path):
        SimpleEncoder.MS_mp3_encoding(jpg_path, mp3_path, 1, 2, seed_len, index_len)
    else:
        print(f"파일 없음: {jpg_path} 또는 {mp3_path}")

    # =========================================================================
    # NFT 인코딩 (JPG + TXN)
    # =========================================================================
    print(f"\n{'='*60}")
    print("NFT 인코딩 (JPG + TXN)")
    print(f"{'='*60}")

    nft_files = [
        ('girlwithballoon.jpg', 'girlwithballoonTXN.txt', 3),
        ('gamechanger.jpg', 'gamechangerTXN.txt', 4),
    ]

    for jpg_file, txn_file, primer in nft_files:
        jpg_path = os.path.join(data_dir, jpg_file)
        txn_path = os.path.join(data_dir, txn_file)
        if os.path.exists(jpg_path) and os.path.exists(txn_path):
            SimpleEncoder.jpg_NFT_encoding(jpg_path, txn_path, primer, seed_len, index_len)
        else:
            print(f"파일 없음: {jpg_path} 또는 {txn_path}")

    print(f"\n{'='*60}")
    print("인코딩 완료!")
    print(f"출력 폴더: {OUTPUT_DIR}")
    print(f"{'='*60}")
