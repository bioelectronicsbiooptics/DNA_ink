# -*- coding: utf-8 -*-
"""
================================================================================
DNA Encoding without Seed (withoutSHA256)
================================================================================

150nt 서열을 꽉 채워서 데이터 저장 효율을 극대화한 인코딩 시스템

핵심 특징:
- 시드 없음 (9nt/18bits 제거) → 150nt 전체를 데이터 저장에 활용
- SHA256 XOR 변환 없음
- Free Energy 검사 없음
- XOR 패리티 포함 (데이터 복구용)

================================================================================
150nt 꽉 채우기 - 시드 제거로 데이터 저장량 증가:
================================================================================

기존 (시드 있음):
  - 시퀀스당 데이터: 124 bits (index 14 + data 110)
  - 시드 공간: 18 bits (9nt) → 낭비

본 버전 (시드 없음):
  - 시퀀스당 데이터: 142 bits (index 14 + data 128)
  - 시드 공간: 0 bits → 전부 데이터로 활용
  - 데이터 저장량 14.5% 증가!

================================================================================
DNA 서열 구조 (150 nt) - 꽉 채운 구조:
================================================================================

┌─────────────────────────────────────────────────────────────────────────────┐
│  Forward Primer (20nt)  │     Payload (110nt)     │  Reverse Primer (20nt)  │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                   ┌───────────────┴───────────────┐
                   │      Data + RS (110nt)        │ ← 시드 없이 전체 데이터!
                   │      (220 bits 전부 활용)      │
                   └───────────────────────────────┘

비교:
- 시드 있는 버전: Seed(9nt) + Data+RS(101nt) = 110nt  → 일부 낭비
- 시드 없는 버전: Data+RS(110nt) = 110nt             → 꽉 채움!

================================================================================
RS 인코딩 구조 (시드 없음 - 150nt 꽉 채움):
================================================================================

입력: index(14bits) + data(142bits) = 156 bits    ← 기존 138bits보다 18bits 증가
RS 부착: 156 bits → 220 bits (+64 bits RS parity)
DNA 변환: 220 bits → 110 nt
Primer 부착: 20nt + 110nt + 20nt = 150nt          ← 꽉 채운 150nt!

================================================================================
"""

import numpy as np
import math
import os
from datetime import datetime

try:
    import galois
    from reedsolo import RSCodec
except ImportError:
    print("필수 라이브러리가 설치되지 않았습니다. 'pip install galois reedsolo' 명령을 실행해주세요.")
    exit()

# =============================================================================
# 설정
# =============================================================================
OUTPUT_DIR = f"Encoded_output_noSeed_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
DATA_DIR = "DNAPrint_Py/data"


# =============================================================================
# 1. 파일 처리
# =============================================================================
class File_tobinary:
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
    @staticmethod
    def deci_to_bi_null(array):
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
        binary = bin(n)[2:]
        if length is not None:
            binary = '0' * (length - len(binary)) + binary
        return np.fromiter(binary, dtype=int)

    @staticmethod
    def transform(arr):
        bases = 'GCAT'
        result_DNA = np.zeros((np.size(arr, axis=0), np.size(arr, axis=1) // 2), dtype='<U3')
        for i in range(np.size(result_DNA, axis=0)):
            for j in range(np.size(result_DNA, axis=1)):
                DNA_index = arr[i, j*2] + 2 * arr[i, j*2+1]
                result_DNA[i, j] = bases[DNA_index]
        return result_DNA


# =============================================================================
# 3. 공통 유틸리티
# =============================================================================
class Common:
    @staticmethod
    def arr_pad(arr, message_len):
        arr_flat = np.ravel(arr)
        pad_cnt = math.ceil(arr.size / message_len) * message_len - arr.size
        for i in range(pad_cnt):
            arr_flat = np.pad(arr_flat, (0, 1), constant_values=0)
        np_arr = arr_flat.reshape(int(arr_flat.size / message_len), message_len)
        return np_arr, pad_cnt

    @staticmethod
    def xor_merge(arr):
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
        """
        Reed-Solomon 심볼 부착

        시드 없는 버전:
        - 입력: 156 bits (index 14 + data 142)
        - 출력: 220 bits (+64 RS parity)
        """
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
        primers = {
            1: ("CTGTCCATAGCCTTGTTCGT", "GCGGAAACGTAGTGAAGGTA"),
            2: ("CAAGTAACCGGCAACAACTG", "ACATAACAACCACCGCGAAA"),
            3: ("TGTATTTCCTTCGGTGCTCC", "TTTCGACAACGGTCTGGTTT"),
            4: ("TCCTCAGCCGATGAAATTCC", "TGTACCATCCGTTTGACTGG"),
            5: ("AAGGCAAGTTGTTACCAGCA", "TGCGACCGTAATCAAACCAA"),
            6: ("GAAGAGTTTAGCCACCTGGT", "AAGGCCAATTCGCGGTTATT"),
            7: ("ATCCTGCAAACGCATTTCCT", "ATGCCTTTCCGAAGTTTCCA"),
            8: ("AATCATGGCCTTCAAACCGT", "AACGCTCCGAAAGTCTTGTT"),
            9: ("AATGGACGTTCCGCAATCAT", "AGAGCCGTGGCAATGTAAAT"),
            10: ("GTCCAGGCAAAGATCCAGTT", "ACCACCGTTAGGCTAAAGTG"),
            11: ("TAGCCTCCAGAATGAAACGG", "TTCAAGCCAAACCGTGTGTA"),
        }
        if select not in primers:
            raise ValueError(f"Invalid primer select: {select}")
        return primers[select]


# =============================================================================
# 4. 프레임 생성 (시드 없음)
# =============================================================================
class To_frags_noSeed:
    """시드 없이 156 bits 전체를 데이터로 사용"""

    @staticmethod
    def file_tofrag_xor(arr, pay_len, index_len):
        """
        시드 없는 프레임 생성 (파일 1개)

        pay_len = 156 (시드 공간 없음)
        mess_len = 156 - 14 = 142 bits 데이터
        """
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
# 5. DNA 인코딩 (시드 없음)
# =============================================================================
class file_toDNA_noSeed:
    """시드 없이 220 비트 전체를 RS 데이터로 사용"""

    @staticmethod
    def encode_noSeed(arr, primer_type):
        """
        시드 없는 DNA 변환

        구조:
        1. 입력: index(14bits) + data(142bits) = 156 bits
        2. RS 부착: 156 bits → 220 bits (+64 bits RS parity)
        3. DNA 변환: 220 bits → 110 nt
        4. Primer 부착: 20nt + 110nt + 20nt = 150nt
        """
        # Reed-Solomon 부착 (156 bits → 220 bits)
        result_bi = Common.Reed(arr, 156)

        F_primers, R_primers = Common.primerSelect(primer_type)

        # DNA 변환 (220 bits 전체가 RS 인코딩된 데이터)
        DNA_seq = Typetransform.transform(result_bi)

        DNA_library = []
        for i in range(len(DNA_seq)):
            tmp_DNA = "".join(DNA_seq[i, :].tolist())
            full_seq = F_primers + tmp_DNA + R_primers
            DNA_library.append(full_seq)

        print(f"총 {len(DNA_library)}개 서열 생성 완료 (시드 없음)")
        return DNA_library


# =============================================================================
# 6. 헤더 저장 (시드 없음)
# =============================================================================
class Header:
    @staticmethod
    def save_txtframe(arr, lang_type, data_pad, all_frags, xor_pad, file_exp):
        """TXT 파일 헤더 저장 (시드 없음)"""
        file_type = Header.file_type_num(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4+4] = Typetransform.deci2bi(lang_type, 4)
        arr[0, 8:8+8] = Typetransform.deci2bi(data_pad, 8)
        arr[0, 16:16+20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 36:36+6] = Typetransform.deci2bi(0, 6)  # seed_len = 0
        arr[0, 42:42+2] = Typetransform.deci2bi(xor_pad, 2)
        print(f"TXT헤더: type={file_type}, lang={lang_type}, pad={data_pad}, frags={all_frags}, seed=0, xor_pad={xor_pad}")
        return arr

    @staticmethod
    def save_MSmp3frame(arr, data_pad, all_frags, xor_pad, file_exp):
        """MP3/JPG 파일 헤더 저장 (시드 없음)"""
        file_type = Header.file_type_num(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4+8] = Typetransform.deci2bi(data_pad, 8)
        arr[0, 12:12+20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 32:32+6] = Typetransform.deci2bi(0, 6)  # seed_len = 0
        arr[0, 38:38+2] = Typetransform.deci2bi(xor_pad, 2)
        print(f"JPG/MP3헤더: type={file_type}, pad={data_pad}, frags={all_frags}, seed=0, xor_pad={xor_pad}")
        return arr

    @staticmethod
    def save_jpgNFTframe(arr, jpg_pad, TXN_pad, jpg_frags, TXN_frags, all_frags, xor_pad, file_exp):
        """NFT 파일 헤더 저장 (JPG + TXN, 시드 없음)"""
        file_type = Header.file_type_num(file_exp)
        arr[0, :4] = Typetransform.deci2bi(file_type, 4)
        arr[0, 4:4+8] = Typetransform.deci2bi(jpg_pad, 8)
        arr[0, 12:12+8] = Typetransform.deci2bi(TXN_pad, 8)
        arr[0, 20:20+20] = Typetransform.deci2bi(jpg_frags, 20)
        arr[0, 40:40+20] = Typetransform.deci2bi(TXN_frags, 20)
        arr[0, 60:60+20] = Typetransform.deci2bi(all_frags, 20)
        arr[0, 80:80+6] = Typetransform.deci2bi(0, 6)  # seed_len = 0
        arr[0, 86:86+2] = Typetransform.deci2bi(xor_pad, 2)
        print(f"NFT헤더: type={file_type}, jpg_pad={jpg_pad}, TXN_pad={TXN_pad}, jpg_frags={jpg_frags}, TXN_frags={TXN_frags}, all_frags={all_frags}")
        return arr

    @staticmethod
    def file_type_num(exp):
        """파일 타입 코드"""
        types = {'txt': 0, 'jpg_TXN': 1, 'mp3': 2, 'jpg': 3}
        return types.get(exp, 0)


# =============================================================================
# 7. 언어 타입 변환
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
# 8. 메인 인코더
# =============================================================================
class NoSeedEncoder:
    """시드 없는 인코딩 클래스 - 150nt 꽉 채움"""

    @staticmethod
    def txt_encoding(path, primer_type, lang_type, index_len):
        """TXT 파일 인코딩 (시드 없음)"""
        print(f"\n{'='*60}")
        print(f"TXT Encoding (no Seed - 150nt 꽉 채움): {path}")
        print(f"{'='*60}")

        txt_np = File_tobinary.file_open(path)

        pay_len = 156  # 시드 공간 없음
        xor_data, pad_count, all_frags, pay_len, index_bits_len, xor_pad = To_frags_noSeed.file_tofrag_xor(
            txt_np, pay_len, index_len
        )

        arr_frags = all_frags // 3 * 2
        only_data = xor_data[:arr_frags, index_bits_len:]

        arr_data = Header.save_txtframe(only_data, lang_type, pad_count, all_frags, xor_pad, 'txt')
        re_xor, rexor_pad = Common.xor_merge(arr_data)
        xor_data[:, index_bits_len:] = re_xor

        DNA_seqs = file_toDNA_noSeed.encode_noSeed(xor_data, primer_type)

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
    def MS_sheet_onefile(path, primer_type, index_len, exp):
        """JPG 또는 MP3 단일 파일 인코딩 (시드 없음)"""
        print(f"\n{'='*60}")
        print(f"{exp.upper()} Encoding (no Seed - 150nt 꽉 채움): {path}")
        print(f"{'='*60}")

        file_np = File_tobinary.file_open(path)

        pay_len = 156  # 시드 공간 없음
        file_xor_data, file_pad_count, file_xor_frags, file_pay_len, file_index_bits_len, file_xor_pad = To_frags_noSeed.file_tofrag_xor(
            file_np, pay_len, index_len
        )

        arr_frags = file_xor_frags // 3 * 2
        file_only_data = file_xor_data[:arr_frags, file_index_bits_len:]

        file_arr_data = Header.save_MSmp3frame(file_only_data, file_pad_count, file_xor_frags, file_xor_pad, exp)
        re_xor, rexor_pad = Common.xor_merge(file_arr_data)
        file_xor_data[:, file_index_bits_len:] = re_xor

        DNA_seqs = file_toDNA_noSeed.encode_noSeed(file_xor_data, primer_type)

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
    def MS_mp3_encoding(jpg_path, mp3_path, jpg_primer, mp3_primer, index_len):
        """악보(JPG) + MP3 인코딩 (시드 없음)"""
        sheet_seqs = NoSeedEncoder.MS_sheet_onefile(jpg_path, jpg_primer, index_len, 'jpg')
        mp3_seqs = NoSeedEncoder.MS_sheet_onefile(mp3_path, mp3_primer, index_len, 'mp3')

        merged_seq = sheet_seqs + mp3_seqs
        save_path = os.path.join(OUTPUT_DIR, "mp3_sheet_merged.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in merged_seq:
                f.write(seq + '\n')

        print(f"\n병합 저장: {save_path}")
        print(f"총 서열 수: {len(merged_seq)}")
        return merged_seq

    @staticmethod
    def jpg_NFT_encoding(jpg_path, TXN_path, primer_type, index_len):
        """NFT 인코딩 (JPG + Transaction ID, 시드 없음)"""
        print(f"\n{'='*60}")
        print(f"NFT Encoding (no Seed - 150nt 꽉 채움): {jpg_path} + {TXN_path}")
        print(f"{'='*60}")

        jpg_np = File_tobinary.file_open(jpg_path)
        TXN_np = File_tobinary.file_open(TXN_path)

        pay_len = 156  # 시드 공간 없음
        xor_data, jpg_pad_count, TXN_pad_count, jpg_frags, TXN_frags, all_frags, pay_len, index_bits_len, xor_pad = To_frags_noSeed.file2_xor_tofrag(
            jpg_np, TXN_np, pay_len, index_len
        )

        arr_frags = all_frags // 3 * 2
        only_data = xor_data[:arr_frags, index_bits_len:]

        arr_data = Header.save_jpgNFTframe(only_data, jpg_pad_count, TXN_pad_count, jpg_frags, TXN_frags, all_frags, xor_pad, 'jpg_TXN')
        re_xor, rexor_pad = Common.xor_merge(arr_data)
        xor_data[:, index_bits_len:] = re_xor

        DNA_seqs = file_toDNA_noSeed.encode_noSeed(xor_data, primer_type)

        name = os.path.splitext(os.path.basename(jpg_path))[0]
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        save_path = os.path.join(OUTPUT_DIR, f"{name}_NFT.txt")

        with open(save_path, 'w', encoding='utf-8') as f:
            for seq in DNA_seqs:
                f.write(seq + '\n')

        print(f"저장 완료: {save_path}")
        print(f"총 서열 수: {len(DNA_seqs)}, 서열 길이: {len(DNA_seqs[0])}nt")
        return DNA_seqs


# =============================================================================
# 9. 실행
# =============================================================================
if __name__ == "__main__":
    print("="*60)
    print("DNA Encoding without Seed (withoutSHA256)")
    print("- 시드 없음 (9nt/18bits 제거) → 150nt 꽉 채움!")
    print("- 110nt 전체를 RS 인코딩 데이터로 사용")
    print("- XOR 패리티 포함")
    print("="*60)

    index_len = 14
    data_dir = DATA_DIR

    print(f"\n설정:")
    print(f"  - Seed: 없음 (0 bits)")
    print(f"  - Index: {index_len} bits")
    print(f"  - Data per sequence: 142 bits (vs 124 bits with seed)")
    print(f"  - 출력 폴더: {OUTPUT_DIR}")

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
            NoSeedEncoder.txt_encoding(filepath, primer, lang_type, index_len)
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
        NoSeedEncoder.MS_mp3_encoding(jpg_path, mp3_path, 1, 2, index_len)
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
            NoSeedEncoder.jpg_NFT_encoding(jpg_path, txn_path, primer, index_len)
        else:
            print(f"파일 없음: {jpg_path} 또는 {txn_path}")

    print(f"\n{'='*60}")
    print("인코딩 완료!")
    print(f"출력 폴더: {OUTPUT_DIR}")
    print(f"{'='*60}")
