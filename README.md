# Printable DNA data storage inks for object-bound digital provenance

Encoding and decoding code accompanying the paper.

Raw sequencing data are deposited in the NCBI Sequence Read Archive under BioProject [**PRJNA1522504**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1522504) (13 runs).

---

## Contents

```
encoding/    digital file  ->  150-nt oligonucleotides
decoding/    NGS reads     ->  reconstructed file
```

---

## encoding/

| File | Description |
|---|---|
| `DNA_encoding_annotated.ipynb` | Annotated source notebook for the encoder |
| `DNA_encoding_full.py` | Batch encoding of all ten datasets |
| `DNA_encoding_with_seeds.py` | Variant with the seed sweep enabled |
| `DNA_encoding_withoutSHA256.py` | **Control** — no SHA-256 masking |
| `DNA_encoding_noSeed_withoutSHA256.py` | **Control** — no seed and no masking |

### Encoding scheme

```
input file  ->  binary  ->  split into 124-bit blocks
  ->  prepend a 14-bit index            = 138-bit fragment
  ->  generate XOR redundancy fragments
  ->  Reed-Solomon encoding (8 parity symbols = 64 bit)
  ->  XOR with a 202-bit mask taken from the SHA-256 hash of an 18-bit seed
  ->  prepend the seed                  = 220-bit codeword
  ->  map 2 bits per base (G=00, C=10, A=01, T=11) = 110-nt payload
  ->  append 20-nt primers at both ends = 150-nt oligonucleotide
```

**Design constraints**, applied to the 110-nt payload only: maximum homopolymer length of 2, GC content of 40–60 %, and a minimum free energy ΔG > −30 kcal mol⁻¹.

To obtain compliant sequences, the 18-bit seeds from 0 to 262,143 are shuffled and applied in turn, and the first seed satisfying all three constraints is adopted for that fragment (seed sweep).

The two control encoders isolate the contribution of the hash masking. A comparison of the three conditions (Hash / NoSeed / WithoutSHA256) is reported in Figure 2 of the paper.

---

## decoding/

| File | Description |
|---|---|
| `DNA_decoding_annotated.ipynb` | Annotated source notebook for the decoder |
| `dna_analysis_v17.py` | Read processing and recovery evaluation |

### Read processing (`dna_analysis_v17.py`)

```
merge forward and reverse FASTQ files
  ->  discard reads shorter than 140 nt or containing ambiguous bases (N)
  ->  match the first 20 nt against the primer and its reverse complement
      (three or fewer mismatches)
  ->  reorient reverse-complement reads
  ->  trim the flanking primers          = 110-nt payload
  ->  group reads by the index at the start of the payload
  ->  take the most frequent sequence per index as the representative
  ->  Reed-Solomon decoding              -> reconstructed data
```

No minimum read count is imposed. Substitutions, insertions and deletions are tallied separately for the primer regions (positions 1–20 and 131–150) and the payload region (positions 21–130).

Supported datasets: `Sheet`, `mp3`, `girl`, `game`, `S-K`, `S-E`, `S-I`, `H-K`, `H-C`, `H-E`.

---

## Environment

| Package | Version |
|---|---|
| Python | 3.13.2 |
| NumPy | 2.2.6 |
| galois | 0.4.7 |
| reedsolo | 1.7.0 |
| ViennaRNA | 2.7.0 |

Sequence selection during encoding used the `RNAfold` executable from ViennaRNA 2.5.0.

---

## License

See `LICENSE`.
