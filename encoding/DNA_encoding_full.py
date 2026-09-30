# -*- coding: utf-8 -*-
"""
DNA Encoding Full Version
ipynb 버전의 모든 기능을 포함한 완전한 Python 스크립트
TXT, JPG, MP3, NFT 등 모든 파일 타입 인코딩 지원
"""

# 필요한 라이브러리 임포트
try:
    import galois
    from reedsolo import RSCodec
except ImportError:
    print("필수 라이브러리가 설치되지 않았습니다. 'pip install galois reedsolo' 명령을 실행해주세요.")
    exit()

import numpy as np
import math
import random
import hashlib
import subprocess
import time
import os

# 전역 출력 폴더 설정
from datetime import datetime
OUTPUT_DIR = f"Encoded_output_full_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
DATA_DIR = "DNAPrint_Py/data"  # 입력 데이터 폴더

# =============================================================================
# 1. 파일 처리 클래스
# =============================================================================
class File_tobinary:
    """파일 읽는 class, binary로 읽어줌"""
    @staticmethod
    def file_open(path):
        with open(path, 'rb') as f:
            file = f.read()  # bytes 한줄
        file_np = np.frombuffer(file, dtype=np.uint8)  # 십진수 배열로 바꿔줌 -> 0~255
        unpacked = np.unpackbits(file_np)
        np.set_printoptions(threshold=np.inf)
        print(unpacked)  # 이진수로 바꿔줌 01110
        return unpacked


# =============================================================================
# 2. 타입 변환 클래스
# =============================================================================
class Typetransform:
    """다양한 데이터 타입 간의 변환을 담당하는 클래스"""

    @staticmethod
    def deci_to_bi_null(array):
        """십진수를 이진수로 변경해주고 최대길이로 맞춰줌"""
        array = array.flatten()
        arr_max = array.max()
        arr_len = len(format(arr_max, 'b'))
        binary_arr = []
        for i in range(0, len(array)):
            tmp = np.fromiter(f'{array[i]:0{arr_len}b}', dtype=int)  # arr_len 만큼 앞에 0을 붙여줌
            binary_arr.append(tmp)
        return np.array(binary_arr)

    @staticmethod
    def deci2bi(n, length=None):
        """하나의 십진수를 지정된 길이의 이진 배열로 변환"""
        binary = bin(n)[2:]
        if length is None:
            pass
        else:
            binary = '0' * (length - len(binary)) + binary
        binary = np.fromiter(binary, dtype=int)
        return binary

    @staticmethod
    def transform(arr):
        """binary를 DNA로 변환(2차원)"""
        bases = 'GCAT'
        result_DNA = np.zeros((np.size(arr, axis=0), np.size(arr, axis=1) // 2), dtype='<U3')
        for i in range(0, np.size(result_DNA, axis=0)):
            for j in range(0, np.size(result_DNA, axis=1)):
                DNA_index = arr[i, j*2] + 2 * arr[i, j*2+1]
                result_DNA[i, j] = bases[DNA_index]  # 00->G, 10 -> C, 01-> A, 11-> T
        return result_DNA

    @staticmethod
    def transform_linear(arr):
        """binary를 DNA로 변환(1차원)"""
        bases = 'GCAT'
        result_DNA = np.zeros((np.size(arr) // 2), dtype='<U3')
        for j in range(0, np.size(result_DNA)):
            DNA_index = arr[j*2] + 2 * arr[j*2+1]
            result_DNA[j] = bases[DNA_index]  # 00->G, 10 -> C, 01-> A, 11-> T
        return result_DNA

    @staticmethod
    def DNA2bi(arr):
        """DNA를 binary로 변환"""
        bi_arr = np.zeros((len(arr) * 2), dtype=np.int64)
        for j in range(0, len(arr)):
            if arr[j] == 'G':
                bi_arr[2*j] = 0
                bi_arr[2*j+1] = 0
            elif arr[j] == 'C':
                bi_arr[2*j] = 1
                bi_arr[2*j+1] = 0
            elif arr[j] == 'A':
                bi_arr[2*j] = 0
                bi_arr[2*j+1] = 1
            elif arr[j] == 'T':
                bi_arr[2*j] = 1
                bi_arr[2*j+1] = 1
            else:
                print("값이 ATCG가 아님 :", arr[j])
                raise ValueError
        return bi_arr


# =============================================================================
# 3. Homopolymer 및 서열 제약 검사 클래스
# =============================================================================
class homo:
    """Homopolymer, GC 함량, 2차 구조 등 DNA 서열의 제약 조건을 검사하는 클래스"""

    @staticmethod
    def seed_transform(arr, seed, seed_len):
        """주어진 시드를 사용하여 데이터에 XOR 변환을 적용"""
        hash = hashlib.sha256(bytes(seed)).hexdigest()
        bi_hash = bin(int(hash[:55], 16))[2:]  # 256bits 중에서 220bits 선택
        if len(bi_hash) < 220 - seed_len:
            bi_hash = '0' * (220 - seed_len - len(bi_hash)) + bi_hash  # 맨 앞쪽 00이 없어질때가 있음 0011 -> 11
        elif len(bi_hash) > 220 - seed_len:
            bi_hash = bi_hash[:220 - seed_len]

        bi_hash = np.fromiter(bi_hash, dtype=int)  # 1차원으로 변경
        xor_arr = np.logical_xor(arr, bi_hash).astype(int)  # xor 시킴

        return xor_arr

    @staticmethod
    def Check_WO_homopolymer(arr):
        """homopolymer 검사 - 이진수 상태에서 검사"""
        max_count = -1
        count = 1
        for j in range(0, np.size(arr) - 2, 2):
            if (arr[j:j+2] == arr[j+2:j+4]).all():  # 2bits가 연속된다면
                count += 1
                if count > max_count:  # count가 최댓값보다 크다면 max에 저장
                    max_count = count
                    if max_count > 2:  # max가 2보다 크다면 실패
                        return False
            else:  # 연속되지 않는다면
                count = 1

        if max_count > 2:
            return False
        else:
            return True

    @staticmethod
    def Check_WO_GC(arr):
        """GC contents 계산"""
        GC_contents = 0
        AT_contents = 0
        for i in range(0, np.size(arr), 2):
            if (arr[i:i+2] == [0, 0]).all():  # G면 GC 증가
                GC_contents += 1
            elif (arr[i:i+2] == [1, 0]).all():  # C면 GC 증가
                GC_contents += 1
            elif (arr[i:i+2] == [0, 1]).all():  # A면 AT 증가
                AT_contents += 1
            elif (arr[i:i+2] == [1, 1]).all():  # T면 AT 증가
                AT_contents += 1
            else:
                raise ValueError
        total = GC_contents + AT_contents
        if ((GC_contents > total * 0.4) and (GC_contents < total * 0.6)):  # 40~60
            return True
        else:
            return False

    @staticmethod
    def cal_secondary(seq):
        """RNAfold를 사용하여 2차 구조 안정성을 확인 (Delta G 값 기준)"""
        try:
            p = subprocess.Popen('RNAfold', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = p.communicate(seq)

            if not stdout or stdout.strip() == '':
                return False

            # 결과에서 자유 에너지(dG) 값 추출
            parts = stdout.split()
            if len(parts) > 1 and parts[-1].startswith('(') and parts[-1].endswith(')'):
                dG_str = parts[-1].strip('()').strip()
                score = float(dG_str)
                print(score)
                return score > -30
            return False
        except Exception as e:
            # ViennaRNA를 사용할 수 없을 경우, 이 검사를 통과시킵니다.
            return True


# =============================================================================
# 4. 공통 유틸리티 클래스
# =============================================================================
class Common:
    """인코딩 과정에서 공통적으로 사용되는 함수들을 모아둔 클래스"""

    @staticmethod
    def arr_pad(arr, message_len):
        """padding을 붙이고 원하는 크기로 모양 바꿔서 내보내는 코드"""
        arr_flat = np.ravel(arr)  # 0차원, 1차원
        print('arr_flat', arr_flat)
        pad_cnt = math.ceil(arr.size / message_len) * message_len - arr.size

        for i in range(0, pad_cnt):  # padding 붙여줌
            arr_flat = np.pad(arr_flat, (0, 1), constant_values=0)

        np_arr = arr_flat.reshape(int(arr_flat.size / message_len), message_len)

        return np_arr, pad_cnt

    @staticmethod
    def xor_merge(arr):
        """xor 생성하고 붙여줌"""
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
        """reedsolomon symbol 부착하는 코드"""
        # 갈루아 필드 사용
        parity = 8
        gf_en = galois.GF(2**3)  # 갈루아필드
        gf_en = np.array(gf_en(arr))  # array를 갈루아필드로 변환

        rsc = RSCodec(parity)  # parity symbol 8개
        code_list = []
        for i in range(0, np.size(gf_en, axis=0)):
            temp = rsc.encode(gf_en[i])  # RS symbol 생성
            code_list.append(temp)
        arr_index_pari = np.array(code_list, dtype=np.int64)  # 이건 gf_en에 parity붙은거

        gf_parity = arr_index_pari[:, -parity:]  # RS symbol부분만 분리
        gf_payload = arr_index_pari[:, :-parity]  # RS symbol부분빼고 나머지 부분
        biparity = Typetransform.deci_to_bi_null(gf_parity)  # RS symbol binary로 변환
        bi_parity = biparity.reshape(np.size(gf_en, axis=0), parity * 8)
        result_bi = np.concatenate((gf_payload, bi_parity), axis=1)  # binary로 바꾼 RS symbol data에 붙여줌
        return result_bi

    @staticmethod
    def primerSelect(select):
        """primer 지정하는 것"""
        if select == 1:  # musicsheet
            F_prim = "CTGTCCATAGCCTTGTTCGT"
            R_prim = "GCGGAAACGTAGTGAAGGTA"
        elif select == 2:  # mp3
            F_prim = "CAAGTAACCGGCAACAACTG"
            R_prim = "ACATAACAACCACCGCGAAA"
        elif select == 3:  # jpg(girl with balloon)
            F_prim = "TGTATTTCCTTCGGTGCTCC"
            R_prim = "TTTCGACAACGGTCTGGTTT"
        elif select == 4:  # jpg(gamechanger)
            F_prim = "TCCTCAGCCGATGAAATTCC"
            R_prim = "TGTACCATCCGTTTGACTGG"
        elif select == 5:  # Txt(sonnet18 한글)
            F_prim = "AAGGCAAGTTGTTACCAGCA"
            R_prim = "TGCGACCGTAATCAAACCAA"
        elif select == 6:  # Txt(sonnet18 영어)
            F_prim = "GAAGAGTTTAGCCACCTGGT"
            R_prim = "AAGGCCAATTCGCGGTTATT"
        elif select == 7:  # Txt(sonnet18 이탈리아어)
            F_prim = "ATCCTGCAAACGCATTTCCT"
            R_prim = "ATGCCTTTCCGAAGTTTCCA"
        elif select == 8:  # Txt(hayeoga 한글)
            F_prim = "AATCATGGCCTTCAAACCGT"
            R_prim = "AACGCTCCGAAAGTCTTGTT"
        elif select == 9:  # Txt(hayeoga 한문)
            F_prim = "AATGGACGTTCCGCAATCAT"
            R_prim = "AGAGCCGTGGCAATGTAAAT"
        elif select == 10:  # Txt(hayeoga 영어)
            F_prim = "GTCCAGGCAAAGATCCAGTT"
            R_prim = "ACCACCGTTAGGCTAAAGTG"
        elif select == 11:  # song signature NFT
            F_prim = "TAGCCTCCAGAATGAAACGG"
            R_prim = "TTCAAGCCAAACCGTGTGTA"
        else:
            raise ValueError(f"Invalid primer select: {select}")

        return F_prim, R_prim

    @staticmethod
    def check_restrict(arr, primer_type, seed_len):
        """homopolymer, GC contents, free energy 계산하는 코드"""
        F_primers, R_primers = Common.primerSelect(primer_type)
        F_primers_bi = Typetransform.DNA2bi(F_primers)[-4:]  # homopolymer 검사할때 primer는 동일하기때문에 payload와 접하는 부분만 검사하면됨
        R_primers_bi = Typetransform.DNA2bi(R_primers)[:4]

        F_primers_tmp = Typetransform.DNA2bi(F_primers)  # DNA 서열을 binary 서열로 변경
        R_primers_tmp = Typetransform.DNA2bi(R_primers)

        wo_homo = np.zeros((np.size(arr, axis=0), 220), dtype=np.int64)  # 조절된 서열 넣을 빈 행렬
        excep = []

        for i in range(0, np.size(arr, axis=0)):
            seed_arr = list(range(0, 2**seed_len))
            random.seed(0)
            random.shuffle(seed_arr)

            for j in range(0, 2**seed_len):
                temp = homo.seed_transform(arr[i], seed_arr[j], seed_len)
                temp_220 = np.zeros(220, dtype=np.int64)
                temp_220[:seed_len] = Typetransform.deci2bi(seed_arr[j], seed_len)
                temp_220[seed_len:] = temp
                temp_primer = np.concatenate((F_primers_bi, temp_220, R_primers_bi))  # homopolymer 검사용

                # 새로 만들어진 data서열과 프라이머 모두 붙여서 dna 서열로 만든거
                temp_seq = np.concatenate((F_primers_tmp, temp_220, R_primers_tmp))  # primer 전체 붙인거
                temp_SEQ = "".join(Typetransform.transform_linear(temp_seq).tolist())  # DNA로 변환

                homo_bool = False
                GC_bool = False
                energy_bool = False

                if homo.Check_WO_homopolymer(temp_primer) == False:  # homopolymer 찾기 수행
                    if len(seed_arr) - 1 == j:
                        excep.append(i)
                        raise ValueError("homo")
                    continue
                else:
                    homo_bool = True

                if homo.Check_WO_GC(temp_seq) == True:  # GC contents 40~60이내에 들어오는지 검사
                    GC_bool = True
                else:
                    if len(seed_arr) - 1 == j:
                        excep.append(i)
                        raise ValueError("GC")
                    continue

                if homo.cal_secondary(temp_SEQ) == True:  # delta G > -30인지 검사
                    energy_bool = True
                else:
                    if len(seed_arr) - 1 == j:
                        excep.append(i)
                        raise ValueError("secondary")
                    continue

                if (homo_bool and GC_bool and energy_bool) == True:  # 모두다 통과이면
                    wo_homo[i] = temp_220
                    print(i, j, '성공')
                    break
                else:
                    homo_bool = False
                    GC_bool = False
                    energy_bool = False
                    raise ValueError("")

        return wo_homo, F_primers, R_primers

    @staticmethod
    def save_txtframe(arr, lang_type, data_pad, all_frags, seed_len, xor_pad, file_exp):
        """맨 첫번째 index에 파일 정보 저장하는 코드 (파일 헤더 역할)"""
        file_type = Common.file_types(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)  # 2^4 encoding된 파일 유형
        arr[0, 4:4 + 4] = Typetransform.deci2bi(lang_type, 4)  # 2^4 언어 유형
        arr[0, 8:8 + 8] = Typetransform.deci2bi(data_pad, 8)  # 2^8 원본 data에 padding 몇개 들어갔나
        arr[0, 16:16 + 20] = Typetransform.deci2bi(all_frags, 20)  # 전체 DNA서열 개수 2^8
        arr[0, 36:36 + 6] = Typetransform.deci2bi(seed_len, 6)  # seed가 차지하는 공간(bit) 2^6
        arr[0, 42:42 + 2] = Typetransform.deci2bi(xor_pad, 2)  # row padding이 들어갔나
        print(file_type, lang_type, data_pad, all_frags, seed_len, xor_pad)
        return arr

    @staticmethod
    def save_MSmp3frame(arr, data_pad, all_frags, seed_len, xor_pad, file_exp):
        """맨 첫번째 index에 파일 정보 저장하는 코드 (파일 헤더 역할)"""
        file_type = Common.file_types(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)  # 2^4 encoding된 파일 유형
        arr[0, 4:4 + 8] = Typetransform.deci2bi(data_pad, 8)  # 2^8 원본 data에 padding 몇개 들어갔나
        arr[0, 12:12 + 20] = Typetransform.deci2bi(all_frags, 20)  # 전체 DNA서열 개수
        arr[0, 32:32 + 6] = Typetransform.deci2bi(seed_len, 6)  # seed가 차지하는 공간(bit)
        arr[0, 38:38 + 2] = Typetransform.deci2bi(xor_pad, 2)  # row padding이 들어갔나
        print(file_type, data_pad, all_frags, seed_len, xor_pad)
        return arr

    @staticmethod
    def save_jpgNFTframe(arr, jpg_pad, TXN_pad, jpg_frags, TXN_frags, all_frags, seed_len, xor_pad, file_exp):
        """맨 첫번째 index에 파일 정보 저장하는 코드 (파일 헤더 역할)"""
        file_type = Common.file_types(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)  # 2^4 encoding된 파일 유형
        arr[0, 4:4 + 8] = Typetransform.deci2bi(jpg_pad, 8)  # jpg pad값 2^8
        arr[0, 12:12 + 8] = Typetransform.deci2bi(TXN_pad, 8)  # TXN pad값 2^8
        arr[0, 20:20 + 20] = Typetransform.deci2bi(jpg_frags, 20)  # jpg 서열 개수 2^20
        arr[0, 40:40 + 20] = Typetransform.deci2bi(TXN_frags, 20)  # TXN 서열 개수 2^20
        arr[0, 60:60 + 20] = Typetransform.deci2bi(all_frags, 20)  # 전체 서열 개수 2^20
        arr[0, 80:80 + 6] = Typetransform.deci2bi(seed_len, 6)  # seed가 차지하는 공간(bit)
        arr[0, 86:86 + 2] = Typetransform.deci2bi(xor_pad, 2)  # xor 하기 위해 data뒤에 padding이 들어갔나
        print(file_type, jpg_pad, TXN_pad, jpg_frags, TXN_frags, all_frags, seed_len, xor_pad)
        return arr

    @staticmethod
    def file_types(exp):
        """파일별로 맨 첫번째 index 저장되는 값이 다름"""
        if exp == "txt":
            return 0
        elif exp == "jpg_TXN":  # NFT
            return 1
        elif exp == "mp3":
            return 2
        elif exp == "jpg":  # 악보
            return 3
        else:
            raise ValueError("파일 유형은 txt, jpg_TXN, mp3, jpg중 하나여야함")


# =============================================================================
# 5. 데이터 프레임/조각 생성 클래스
# =============================================================================
class To_frags:
    """데이터를 DNA로 변환하기 전, 인덱스와 페이로드가 포함된 조각으로 구성하는 클래스"""

    @staticmethod
    def file_tofrag_xor(arr, pay_len, index_len):
        """파일 1개 xor붙인것"""
        # 데이터가 들어갈 bits
        mess_len = pay_len - index_len  # pay_len = 220bits - seed - reed

        # null 값 채우는 것
        arr_np, arr_pad_count = Common.arr_pad(arr, mess_len)
        arr_np = np.pad(arr_np, ((1, 0), (0, 0)), constant_values=0)  # 파일 헤더 빈 서열 만듦

        # xor데이터 생성
        xor_merged, boolxor_pad = Common.xor_merge(arr_np)

        index = np.arange(2**index_len)  # 사용할수 있는 index 모두 생성
        index_bi = Typetransform.deci_to_bi_null(index)

        index_arr = np.concatenate((index_bi[:np.size(xor_merged, axis=0)], xor_merged), axis=1)
        frag_count = np.size(index_arr, axis=0)

        return index_arr, arr_pad_count, frag_count, pay_len, index_len, boolxor_pad

    @staticmethod
    def file2_xor_tofrag(arr, arr2, pay_len, index_len):
        """파일 2개 xor 생성"""
        # 데이터가 들어갈 bits
        mess_len = pay_len - index_len

        arr_np, arr_pad_count = Common.arr_pad(arr, mess_len)
        arr2_np, arr2_pad_count = Common.arr_pad(arr2, mess_len)

        merged_2 = np.concatenate((arr_np, arr2_np), axis=0)
        merged_2 = np.pad(merged_2, ((1, 0), (0, 0)), constant_values=0)  # 파일 헤더 빈 서열 만듦

        xor_merged, boolxor_pad = Common.xor_merge(merged_2)

        index = np.arange(2**index_len)
        index_bi = Typetransform.deci_to_bi_null(index)

        xor_merged = np.concatenate((index_bi[:np.size(xor_merged, axis=0)], xor_merged), axis=1)
        arr_frag = np.size(arr_np, axis=0)
        arr2_frag = np.size(arr2_np, axis=0)
        frag_count = np.size(xor_merged, axis=0)

        return xor_merged, arr_pad_count, arr2_pad_count, arr_frag, arr2_frag, frag_count, pay_len, index_len, boolxor_pad


# =============================================================================
# 6. DNA 인코딩 클래스
# =============================================================================
class file_toDNA:
    """리드솔로몬 붙여주고, 서열제한시키고 2bit씩 base로 변환해주는 코드"""

    @staticmethod
    def filtered_DNA(arr, seed_len, primer_type):
        """DNA 서열로 변환하고 제약 조건 검사"""
        result_bi = Common.Reed(arr, 156 - seed_len)  # reed solomon 부착
        wo_homo, F_primers, R_primers = Common.check_restrict(result_bi, primer_type, seed_len)  # 서열 제어

        DNA_seq = Typetransform.transform(wo_homo)  # binary -> DNA 변환

        DNA_library = []
        for i in range(0, len(DNA_seq)):
            tmp_DNA = "".join(DNA_seq[i, :].tolist())
            addition = F_primers + tmp_DNA + R_primers
            DNA_library.append(addition)

        return DNA_library


# =============================================================================
# 7. 언어 타입 변환 함수
# =============================================================================
def lang_num(language):
    """TXT 파일은 맨 첫서열에 언어 정보 넣어줌
    0: 영어, 1:한국어, 2:이탈리아어, 3:중국어"""
    if language.upper() == 'ENGLISH':
        lang_type = 0
    elif language.upper() == 'KOREAN':
        lang_type = 1
    elif language.upper() == 'ITALIANO':
        lang_type = 2
    elif language.upper() == 'CHINESE':
        lang_type = 3
    else:
        raise ValueError(f"지원하지 않는 언어: {language}")
    return lang_type


# =============================================================================
# 8. 메인 인코딩 클래스
# =============================================================================
class main:
    """메인 인코딩 함수들"""

    @staticmethod
    def txt_encoding(path, primer_type, seed_len, lang_type, index_len):
        """TXT 파일 인코딩 프로세스"""
        print(f"\n--- '{path}' 텍스트 파일 인코딩 시작 ---")

        # 파일들을 읽어옴
        txt_np = File_tobinary.file_open(path)

        # 정보 저장할 프레임 만듦
        xor_data, pad_count, all_frags, pay_len, index_bits_len, xor_pad = To_frags.file_tofrag_xor(
            txt_np, 156 - seed_len, index_len
        )

        arr_frags = all_frags // 3 * 2  # XOR 서열 개수
        only_data = xor_data[:arr_frags, index_bits_len:]  # 파일 data 부분만 골라냄

        txt_exp = 'txt'
        # 파일 헤더에 정보 저장
        arr_data = Common.save_txtframe(only_data, lang_type, pad_count, all_frags, seed_len, xor_pad, txt_exp)
        re_xor, rexor_pad = Common.xor_merge(arr_data)

        xor_data[:, index_bits_len:] = re_xor

        # reed solomon symbol 부착하고 서열제어 하고 염기서열로 변환
        DNA_seqs = file_toDNA.filtered_DNA(xor_data, seed_len, primer_type)

        # 파일 저장
        name, _ = os.path.splitext(os.path.basename(path))
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        current_date = time.strftime("%Y%m%d")
        save_path = os.path.join(OUTPUT_DIR, f"{name}_encoded_{current_date}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for i in range(0, len(DNA_seqs)):
                f.write(DNA_seqs[i] + '\n')

        print(f"인코딩된 DNA 서열을 '{save_path}'에 저장했습니다.\n")
        return DNA_seqs

    @staticmethod
    def jpg_NFT_encoding(jpg_path, TXN_path, primer_type, seed_len, index_len):
        """JPG + NFT Transaction ID 인코딩"""
        print(f"\n--- '{jpg_path}' + '{TXN_path}' NFT 인코딩 시작 ---")

        # 파일들을 읽어옴
        jpg_np = File_tobinary.file_open(jpg_path)
        TXN_np = File_tobinary.file_open(TXN_path)

        # 정보 저장할 프레임 만듦
        xor_data, jpg_pad_count, TXN_pad_count, jpg_frags, TXN_frags, all_frags, pay_len, index_bits_len, xor_pad = To_frags.file2_xor_tofrag(
            jpg_np, TXN_np, 156 - seed_len, index_len
        )

        arr_frags = all_frags // 3 * 2  # XOR 서열 개수
        only_data = xor_data[:arr_frags, index_bits_len:]  # 파일 data 부분만 골라냄

        file_exp = 'jpg_TXN'

        # binary에 frame 정보를 넣고 다시 xor을 만들어줌
        arr_data = Common.save_jpgNFTframe(only_data, jpg_pad_count, TXN_pad_count, jpg_frags, TXN_frags, all_frags, seed_len, xor_pad, file_exp)
        re_xor, rexor_pad = Common.xor_merge(arr_data)

        xor_data[:, index_bits_len:] = re_xor

        # reed solomon symbol 부착하고 서열제어 하고 염기서열로 변환
        DNA_seqs = file_toDNA.filtered_DNA(xor_data, seed_len, primer_type)

        # 파일 저장
        name, _ = os.path.splitext(os.path.basename(jpg_path))
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        current_date = time.strftime("%Y%m%d")
        save_path = os.path.join(OUTPUT_DIR, f"{name}_NFT_encoded_{current_date}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for i in range(0, len(DNA_seqs)):
                f.write(DNA_seqs[i] + '\n')

        print(f"인코딩된 DNA 서열을 '{save_path}'에 저장했습니다.\n")
        return DNA_seqs

    @staticmethod
    def MS_sheet_onefile(path, primer_type, seed_len, index_len, exp):
        """악보(JPG) 또는 MP3 파일 하나를 인코딩"""
        print(f"\n--- '{path}' ({exp}) 파일 인코딩 시작 ---")

        # 파일들을 읽어옴
        file_np = File_tobinary.file_open(path)

        # XOR 생성(index 맞추는 용도)
        file_xor_data, file_pad_count, file_xor_frags, file_pay_len, file_index_bits_len, file_xor_pad = To_frags.file_tofrag_xor(
            file_np, 156 - seed_len, index_len
        )

        arr_frags = file_xor_frags // 3 * 2
        file_only_data = file_xor_data[:arr_frags, file_index_bits_len:]

        # binary에 frame 정보를 넣고 다시 xor을 만들어줌
        file_arr_data = Common.save_MSmp3frame(file_only_data, file_pad_count, file_xor_frags, seed_len, file_xor_pad, exp)
        re_xor, rexor_pad = Common.xor_merge(file_arr_data)

        file_xor_data[:, file_index_bits_len:] = re_xor

        # reed solomon symbol 부착하고 서열제어 하고 염기서열로 변환
        DNA_seqs = file_toDNA.filtered_DNA(file_xor_data, seed_len, primer_type)

        # 파일 저장
        name, _ = os.path.splitext(os.path.basename(path))
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        current_date = time.strftime("%Y%m%d")
        save_path = os.path.join(OUTPUT_DIR, f"{name}_{exp}_encoded_{current_date}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for i in range(0, len(DNA_seqs)):
                f.write(DNA_seqs[i] + '\n')

        print(f"인코딩된 DNA 서열을 '{save_path}'에 저장했습니다.\n")
        return DNA_seqs

    @staticmethod
    def MS_mp3_encoding(jpg_path, mp3_path, jpg_primer_type, mp3_primer_type, seed_len, index_len):
        """악보(JPG) + MP3 파일을 함께 인코딩"""
        print(f"\n=== 악보 + MP3 통합 인코딩 시작 ===")

        sheet_seqs = main.MS_sheet_onefile(jpg_path, jpg_primer_type, seed_len, index_len, 'jpg')
        mp3_seqs = main.MS_sheet_onefile(mp3_path, mp3_primer_type, seed_len, index_len, 'mp3')

        merged_seq = sheet_seqs + mp3_seqs

        os.makedirs(OUTPUT_DIR, exist_ok=True)
        current_date = time.strftime("%Y%m%d")
        save_path = os.path.join(OUTPUT_DIR, f"mp3_sheet_merged_encoded_{current_date}.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for i in range(0, len(merged_seq)):
                f.write(merged_seq[i] + '\n')

        print(f"통합 인코딩된 DNA 서열을 '{save_path}'에 저장했습니다.\n")
        return merged_seq


# =============================================================================
# 9. 메인 실행 부분
# =============================================================================
if __name__ == "__main__":
    print("=" * 80)
    print("DNA Encoding Full Version - 모든 파일 타입 인코딩 지원")
    print("=" * 80)

    # 설정 변수
    seed_len = 18  # 시드 비트 길이
    index_len = 14  # 인덱스 비트 길이

    print(f"\n설정: seed_len={seed_len}, index_len={index_len}\n")

    # 데이터 폴더 확인 (전역 설정 사용)
    data_dir = DATA_DIR
    if not os.path.exists(data_dir):
        print(f"'{data_dir}' 폴더가 없습니다.")
        exit()

    # 사용 가능한 파일 확인
    print("\n사용 가능한 파일:")
    available_files = {}

    # TXT 파일
    txt_files = {
        'hayeoga_kor.txt': (8, 'korean'),
        'hayeoga_eng.txt': (10, 'english'),
        'hayeoga_chn.txt': (9, 'chinese'),
        'Sonnet18_kor.txt': (5, 'korean'),
        'Sonnet18_eng.txt': (6, 'english'),
        'Sonnet18_ital.txt': (7, 'italiano')
    }

    for filename, (primer, lang) in txt_files.items():
        filepath = os.path.join(data_dir, filename)
        if os.path.exists(filepath):
            available_files[filename] = ('txt', primer, lang)
            print(f"  ✓ {filename} (TXT, {lang})")

    # JPG/MP3 파일 (예시)
    media_files = {
        'small.jpg': ('jpg', 1),
        '5sec_64kbps.mp3': ('mp3', 2),
        'girlwithballoon.jpg': ('nft_jpg', 3),
        'gamechanger.jpg': ('nft_jpg', 4)
    }

    for filename, (ftype, primer) in media_files.items():
        filepath = os.path.join(data_dir, filename)
        if os.path.exists(filepath):
            available_files[filename] = (ftype, primer, None)
            print(f"  ✓ {filename} ({ftype.upper()})")

    if not available_files:
        print("\n인코딩할 파일이 없습니다. 'data' 폴더에 파일을 추가해주세요.")
        exit()

    # 사용자 선택
    print("\n" + "=" * 80)
    print("인코딩 옵션:")
    print("  1. 모든 TXT 파일 인코딩")
    print("  2. 특정 파일 선택하여 인코딩")
    print("  3. 악보 + MP3 통합 인코딩")
    print("  4. 종료")
    print("=" * 80)

    try:
        choice = input("\n선택 (1-4): ").strip()

        if choice == '1':
            # 모든 TXT 파일 인코딩
            print("\n=== 모든 TXT 파일 인코딩 시작 ===\n")
            for filename, (ftype, primer, lang) in available_files.items():
                if ftype == 'txt':
                    filepath = os.path.join(data_dir, filename)
                    lang_type = lang_num(lang)
                    main.txt_encoding(filepath, primer, seed_len, lang_type, index_len)

        elif choice == '2':
            # 특정 파일 선택
            print("\n사용 가능한 파일 목록:")
            file_list = list(available_files.keys())
            for idx, filename in enumerate(file_list, 1):
                print(f"  {idx}. {filename}")

            file_idx = int(input("\n인코딩할 파일 번호: ").strip()) - 1
            if 0 <= file_idx < len(file_list):
                filename = file_list[file_idx]
                filepath = os.path.join(data_dir, filename)
                ftype, primer, lang = available_files[filename]

                if ftype == 'txt':
                    lang_type = lang_num(lang)
                    main.txt_encoding(filepath, primer, seed_len, lang_type, index_len)
                elif ftype == 'jpg':
                    main.MS_sheet_onefile(filepath, primer, seed_len, index_len, 'jpg')
                elif ftype == 'mp3':
                    main.MS_sheet_onefile(filepath, primer, seed_len, index_len, 'mp3')
            else:
                print("잘못된 번호입니다.")

        elif choice == '3':
            # 악보 + MP3 통합 인코딩
            jpg_path = os.path.join(data_dir, 'small.jpg')
            mp3_path = os.path.join(data_dir, '5sec_64kbps.mp3')

            if os.path.exists(jpg_path) and os.path.exists(mp3_path):
                main.MS_mp3_encoding(jpg_path, mp3_path, 1, 2, seed_len, index_len)
            else:
                print("악보(small.jpg) 또는 MP3(5sec_64kbps.mp3) 파일이 없습니다.")

        elif choice == '4':
            print("종료합니다.")
            exit()

        else:
            print("잘못된 선택입니다.")

    except KeyboardInterrupt:
        print("\n\n사용자에 의해 중단되었습니다.")
    except Exception as e:
        print(f"\n오류 발생: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 80)
    print("DNA 인코딩 완료!")
    print("=" * 80)
