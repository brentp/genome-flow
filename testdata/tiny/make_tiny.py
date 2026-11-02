#!/usr/bin/env python3
"""Synthetic no-alt FASTA + PE FASTQs for the align-pipe smoke test."""

from __future__ import annotations

import gzip
import random
from pathlib import Path

SEED = 7
READ_LEN = 150
INSERT = 400
N_UNIQUE = 4000
N_PCR_DUP = 200
N_OPTICAL_DUP = 80
CONTIG_LEN = 80_000

BASES = "ACGT"


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def qual(n: int) -> str:
    return "I" * n


def illumina_name(flowcell: str, lane: int, tile: int, x: int, y: int, end: int) -> str:
    return f"A00421:12:{flowcell}:{lane}:{tile}:{x}:{y} {end}:N:0:CGATGT"


def sra_name(i: int, end: int) -> str:
    return f"SRR000001.{i} {end}"


def write_fa(path: Path, contigs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def write_fq(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{qual(len(seq))}\n")


def sample_fragment(rng: random.Random, contigs: list[tuple[str, str]]) -> tuple[str, str]:
    name, seq = rng.choice(contigs)
    start = rng.randrange(0, len(seq) - INSERT)
    frag = seq[start : start + INSERT]
    return frag[:READ_LEN], rc(frag[-READ_LEN:])


def main() -> None:
    out = Path(__file__).resolve().parent
    rng = random.Random(SEED)

    contigs = {
        "chrTiny1": "".join(rng.choice(BASES) for _ in range(CONTIG_LEN)),
        "chrTiny2": "".join(rng.choice(BASES) for _ in range(CONTIG_LEN)),
    }
    write_fa(out / "ref.fa", contigs)
    contig_list = list(contigs.items())

    r1_ill: list[tuple[str, str]] = []
    r2_ill: list[tuple[str, str]] = []
    r1_sra: list[tuple[str, str]] = []
    r2_sra: list[tuple[str, str]] = []

    n = 0
    for _ in range(N_UNIQUE):
        n += 1
        r1, r2 = sample_fragment(rng, contig_list)
        tile, x, y = rng.randint(1101, 2216), rng.randint(1000, 30000), rng.randint(1000, 30000)
        r1_ill.append((illumina_name("HXXXXXXXX", 1, tile, x, y, 1), r1))
        r2_ill.append((illumina_name("HXXXXXXXX", 1, tile, x, y, 2), r2))
        r1_sra.append((sra_name(n, 1), r1))
        r2_sra.append((sra_name(n, 2), r2))

    for i in range(N_PCR_DUP):
        n += 1
        src = rng.randrange(N_UNIQUE)
        r1, r2 = r1_ill[src][1], r2_ill[src][1]
        tile, x, y = rng.randint(1101, 2216), rng.randint(1000, 30000), rng.randint(1000, 30000)
        r1_ill.append((illumina_name("HXXXXXXXX", 1, tile, x, y, 1), r1))
        r2_ill.append((illumina_name("HXXXXXXXX", 1, tile, x, y, 2), r2))
        r1_sra.append((sra_name(n, 1), r1))
        r2_sra.append((sra_name(n, 2), r2))

    for i in range(N_OPTICAL_DUP):
        n += 1
        src = rng.randrange(N_UNIQUE)
        r1, r2 = r1_ill[src][1], r2_ill[src][1]
        # parse original coords and nudge x/y within an optical radius
        parts = r1_ill[src][0].split()[0].split(":")
        tile, x, y = int(parts[4]), int(parts[5]), int(parts[6])
        r1_ill.append((illumina_name("HXXXXXXXX", 1, tile, x + 2, y + 1, 1), r1))
        r2_ill.append((illumina_name("HXXXXXXXX", 1, tile, x + 2, y + 1, 2), r2))
        r1_sra.append((sra_name(n, 1), r1))
        r2_sra.append((sra_name(n, 2), r2))

    write_fq(out / "illumina.R1.fq.gz", r1_ill)
    write_fq(out / "illumina.R2.fq.gz", r2_ill)
    write_fq(out / "simulated.R1.fq.gz", r1_sra)
    write_fq(out / "simulated.R2.fq.gz", r2_sra)
    print(
        f"wrote {out}/ref.fa and {n} PE pairs "
        f"({N_UNIQUE} unique + {N_PCR_DUP} PCR dups + {N_OPTICAL_DUP} optical dups)"
    )


if __name__ == "__main__":
    main()
