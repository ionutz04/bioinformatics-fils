'''
Take an arbitrary DNA sequence from NCBI (bacterial genome).  Use 5 restriction enzymes (name, recognized site and cleaveage site).

        5'GAATTC        5'---G      AATTC---3'
EcoR1   3'CTTAAG        3'---CTTAA  G---5'

        5'GGATCC        5'---G      GATCC---3'
BamHl   3'CCTAGG        3'---CCTAG  G---5'

        5'AAGCTT        5'---A      AGCTT---3'
Himdlll 3'TTCGAA        3'---TTCGA  A---5'

        5'TCGA          5'---T      CGA---3'
Taql    3'AGCT          3'---AGC    T---5'

        5'GGCC          5'---GG     CC---3'
Haelll* 3'CCGG          3'---CC     GG---5'

Input of the implementation:
1. The recognized sequence for each restriction enzyme.
2. The DNA sequence to be digested.

Output of the implementation:
1. Number of cleaveages ofr each restriction enzyme.
2. Cleaveage positions and length of fragments.
A simulation of the electrophoresis gen based on the number of restriction enzymes used.
'''
import math

import matplotlib.pyplot as plt

# name -> (recognized site 5'->3', cut offset inside the site on the top strand)
ENZYMES = {
    "EcoRI":   ("GAATTC", 1),  # G^AATTC
    "BamHI":   ("GGATCC", 1),  # G^GATCC
    "HindIII": ("AAGCTT", 1),  # A^AGCTT
    "TaqI":    ("TCGA",   1),  # T^CGA
    "HaeIII":  ("GGCC",   2),  # GG^CC
}


def read_fasta(path):
    with open(path) as f:
        lines = f.read().split("\n")
    return lines[0][1:], "".join(line.strip() for line in lines[1:]).upper()


def find_cuts(dna, site, offset):
    """Return the cut positions (0-based, cut happens before this index)."""
    cuts = []
    i = dna.find(site)
    while i != -1:
        cuts.append(i + offset)
        i = dna.find(site, i + 1)  # +1 so overlapping sites are also found
    return cuts


def fragments(length, cuts):
    """Lengths of the fragments of a linear DNA molecule cut at the given positions."""
    bounds = [0] + sorted(set(cuts)) + [length]
    return [b - a for a, b in zip(bounds, bounds[1:]) if b > a]


def simulate_gel(lanes, out_file="gel.png"):
    """Draw an agarose gel: migration distance is proportional to -log10(fragment length)."""
    ladder = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 300, 200, 100]
    lanes = {"Ladder": ladder, **lanes}

    lo, hi = math.log10(50), math.log10(15000)

    def y(bp):  # 0 = well (top), 1 = bottom of the gel
        return (hi - math.log10(bp)) / (hi - lo)

    fig, ax = plt.subplots(figsize=(1.4 * len(lanes) + 1, 7))
    ax.set_facecolor("#111")
    for x, (name, frags) in enumerate(lanes.items()):
        ax.add_patch(plt.Rectangle((x - 0.35, -0.04), 0.7, 0.025, color="#444"))  # well
        biggest = max(frags) if frags else 1
        for bp in frags:
            if bp < 50:  # very small fragments run off the gel
                continue
            # bigger fragments carry more DNA -> brighter band
            alpha = 0.35 + 0.65 * (bp / biggest) if name != "Ladder" else 0.8
            ax.plot([x - 0.3, x + 0.3], [y(bp)] * 2, color="#7CFC00", lw=3, alpha=alpha)
        if name == "Ladder":
            for bp in ladder:
                ax.text(x - 0.45, y(bp), f"{bp}", color="white", fontsize=7, ha="right", va="center")

    ax.set_xticks(range(len(lanes)), lanes.keys(), rotation=30)
    ax.set_xlim(-1, len(lanes) - 0.5)
    ax.set_ylim(1.05, -0.08)
    ax.set_yticks([])
    ax.set_title("Simulated agarose gel electrophoresis")
    fig.tight_layout()
    fig.savefig(out_file, dpi=150)
    print(f"\nGel image saved to {out_file}")
    plt.show()


def main():
    header, dna = read_fasta("sequence.fasta")
    print(f"Sequence: {header}")
    print(f"Length:   {len(dna)} bp\n")

    lanes = {}
    all_cuts = []
    for name, (site, offset) in ENZYMES.items():
        cuts = find_cuts(dna, site, offset)
        frags = fragments(len(dna), cuts)
        lanes[name] = frags
        all_cuts += cuts

        print(f"{name} ({site[:offset]}^{site[offset:]})")
        print(f"  Number of cleavages: {len(cuts)}")
        print(f"  Cleavage positions:  {cuts}")
        print(f"  Fragment lengths:    {frags}\n")

    lanes["All 5"] = fragments(len(dna), all_cuts)
    simulate_gel(lanes)


if __name__ == "__main__":
    main()