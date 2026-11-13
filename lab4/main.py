#Make an app which samples the genome in the fasta file such that 20000 samples are taken above the genome.
#  Each sample must have a length between 100 and 150 letters.

# design an app that is able to use the 20k samples in order to reconstruct the original genome. The case is that you consider that you don't have access 
# to the original genome.
import random

import numpy as np

AMINO_ACID_TO_CODONS = {
    "Phe":  ["UUU", "UUC"],
    "Leu":  ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "Ile":  ["AUU", "AUC", "AUA"],
    "Met":  ["AUG"],                     
    "Val":  ["GUU", "GUC", "GUA", "GUG"],
    "Ser":  ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "Pro":  ["CCU", "CCC", "CCA", "CCG"],
    "Thr":  ["ACU", "ACC", "ACA", "ACG"],
    "Ala":  ["GCU", "GCC", "GCA", "GCG"],
    "Tyr":  ["UAU", "UAC"],
    "His":  ["CAU", "CAC"],
    "Gln":  ["CAA", "CAG"],
    "Asn":  ["AAU", "AAC"],
    "Lys":  ["AAA", "AAG"],
    "Asp":  ["GAU", "GAC"],
    "Glu":  ["GAA", "GAG"],
    "Cys":  ["UGU", "UGC"],
    "Trp":  ["UGG"],
    "Arg":  ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Gly":  ["GGU", "GGC", "GGA", "GGG"],
    "Stop": ["UAA", "UAG", "UGA"],
}

CODON_TO_AMINO_ACID = {
    codon: aa for aa, codons in AMINO_ACID_TO_CODONS.items() for codon in codons
}

sequence = []

with open("ecoli_cft073.fasta") as f:
    lines = f.read().split('\n')
    genome = ''.join(line.strip() for line in lines[1:]).upper()  # Skip the first line (header) and join the rest
start = 0
length = len(genome) // 20000  # Calculate the length of each sample based on the genome length and number of samples
print(length)
for i in range(20000):
    # if start + length > len(genome):
    #     length = len(genome) - start
    sample = genome[start:start + length]
    start += length
    sequence.append(sample)

print(len(sequence))
print(sequence[:2])
seq2 = sequence.copy()
#shuffle the sequence to randomize the order of samples
random.shuffle(seq2)

random_genome = ''.join(seq2)


def compute_overlap(seq_a, seq_b, min_overlap=20):
    # Slide seq_b along seq_a and, for every offset, compare the overlapping
    # characters one by one (not seq_a == seq_b as whole substrings).
    # offset > 0: seq_b starts inside seq_a (overlap in the last part of seq_a)
    # offset < 0: seq_b starts before seq_a (overlap in the first part of seq_a)
    a = np.frombuffer(seq_a.encode(), dtype=np.uint8)
    b = np.frombuffer(seq_b.encode(), dtype=np.uint8)
    offsets = []
    rates = []
    best = None
    for offset in range(-(len(b) - 1), len(a)):
        a_start = max(0, offset)
        b_start = max(0, -offset)
        overlap_len = min(len(a) - a_start, len(b) - b_start)
        matches = np.count_nonzero(a[a_start:a_start + overlap_len] == b[b_start:b_start + overlap_len])
        offsets.append(offset)
        rates.append(int(matches) / overlap_len)
        # Ignore offsets where only a few chars overlap, their rate is just noise
        if overlap_len >= min_overlap and (best is None or rates[-1] > rates[best]):
            best = len(rates) - 1

    return offsets, rates, best


def plot_overlap(seq_a, seq_b, min_overlap=20):
    offsets, rates, best = compute_overlap(seq_a, seq_b, min_overlap)
    best_offset = offsets[best]

    # Per-character match along seq_a at the best offset (1 = same char, 0 = different)
    a_start = max(0, best_offset)
    b_start = max(0, -best_offset)
    overlap_len = min(len(seq_a) - a_start, len(seq_b) - b_start)
    positions = list(range(a_start, a_start + overlap_len))
    char_match = [int(seq_a[p] == seq_b[b_start + k]) for k, p in enumerate(positions)]

    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7))

    ax1.plot(offsets, rates)
    ax1.axhline(0.25, color="gray", linestyle="--", label="random (0.25)")
    ax1.axvline(best_offset, color="red", linestyle=":", label=f"best offset = {best_offset}")
    ax1.set_xlabel("offset of sequence 2 relative to sequence 1")
    ax1.set_ylabel("overlap rate (matching chars / overlap length)")
    ax1.set_title("Character overlap rate of the first 2 sequences")
    ax1.legend()

    ax2.bar(positions, char_match, width=1.0)
    ax2.set_xlabel("position in sequence 1")
    ax2.set_ylabel("match")
    ax2.set_yticks([0, 1])
    ax2.set_title(f"Per-character matches at offset {best_offset} "
                  f"(rate = {rates[best]:.2%}, overlap = {overlap_len} chars)")

    plt.tight_layout()
    # plt.savefig("output.png")
    # plt.show()


    return best_offset, rates[best]


print(plot_overlap(sequence[0], sequence[1]))



    
def plot_overlap_for_entire_sequance_array(sequences, min_overlap=20, min_rate=0.9):
    # Sort the sequences using the character overlap between every pair:
    # overlap in the first part of sequence i  -> sequence j goes before i
    # overlap in the last part of sequence i   -> sequence j goes after i
    edges = []  # (rate, overlap_len, before, after)
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            offsets, rates, best = compute_overlap(sequences[i], sequences[j], min_overlap)
            if best is None or rates[best] < min_rate:
                continue
            offset = offsets[best]
            overlap_len = min(len(sequences[i]) - max(0, offset), len(sequences[j]) - max(0, -offset))
            if offset < 0:
                edges.append((rates[best], overlap_len, j, i))
            else:
                edges.append((rates[best], overlap_len, i, j))

    # Link the best overlaps first; each sequence gets at most one neighbour
    # before it and one after it, and no link may close a cycle
    next_seq = [None] * len(sequences)
    prev_seq = [None] * len(sequences)
    for rate, overlap_len, before, after in sorted(edges, reverse=True):
        if next_seq[before] is not None or prev_seq[after] is not None:
            continue
        k = after
        while k is not None and k != before:
            k = next_seq[k]
        if k == before:
            continue
        next_seq[before] = after
        prev_seq[after] = before

    # Walk every chain from its first sequence (the one with nothing before it)
    order = []
    for start in range(len(sequences)):
        if prev_seq[start] is None:
            k = start
            while k is not None:
                order.append(k)
                k = next_seq[k]

    return [sequences[k] for k in order]


def merge_sequences(sorted_seqs, min_overlap=20, min_rate=0.9):
    # Glue the sorted sequences into one string, writing the overlapping part only once
    merged = sorted_seqs[0]
    for prev, cur in zip(sorted_seqs, sorted_seqs[1:]):
        offsets, rates, best = compute_overlap(prev, cur, min_overlap)
        if best is not None and rates[best] >= min_rate and offsets[best] >= 0:
            overlap_len = len(prev) - offsets[best]
            merged += cur[overlap_len:]
        else:
            merged += cur
    return merged


def compare_with_original(sorted_seqs, original_seqs, min_overlap=20, min_rate=0.9):
    # original_seqs = the same samples, in the order they appear in the genome
    position = {s: k for k, s in enumerate(original_seqs)}
    order = [position[s] for s in sorted_seqs]
    same_place = sum(k == p for k, p in enumerate(order))
    correct_neighbours = sum(b == a + 1 for a, b in zip(order, order[1:]))
    print(f"Samples in the same place as in the original: {same_place}/{len(order)}")
    print(f"Neighbours kept from the original: {correct_neighbours}/{len(order) - 1}")

    # Character by character comparison of the rebuilt genome with the original one
    rebuilt = merge_sequences(sorted_seqs, min_overlap, min_rate)
    original = merge_sequences(original_seqs, min_overlap, min_rate)
    compare_len = min(len(rebuilt), len(original))
    char_match = [int(rebuilt[k] == original[k]) for k in range(compare_len)]
    print(f"Rebuilt length: {len(rebuilt)}, original length: {len(original)}")
    print(f"Matching chars: {sum(char_match)}/{compare_len} ({sum(char_match) / compare_len:.2%})")

    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7))

    ax1.plot(range(len(order)), order, marker=".", linestyle="")
    ax1.plot([0, len(order) - 1], [0, len(order) - 1], color="gray", linestyle="--", label="perfect order")
    ax1.set_xlabel("position in the sorted sequences")
    ax1.set_ylabel("position in the original sequences")
    ax1.set_title("Sorted order vs original order")
    ax1.legend()

    ax2.plot(range(compare_len), char_match, linewidth=0.5)
    ax2.set_xlabel("position in the genome")
    ax2.set_ylabel("match")
    ax2.set_yticks([0, 1])
    ax2.set_title(f"Rebuilt vs original genome, char by char "
                  f"({sum(char_match) / compare_len:.2%} matching)")

    plt.tight_layout()
    # plt.savefig("comparison.png")
    # plt.show()

    return same_place, correct_neighbours, sum(char_match) / compare_len


# O(n^2) pairs, so only run it on a subset of the samples: take the first 200
# (a continuous part of the genome) and shuffle them, so there is an original to compare with
original_subset = sequence[:200]
shuffled_subset = original_subset.copy()
random.shuffle(shuffled_subset)

sorted_sequences = plot_overlap_for_entire_sequance_array(shuffled_subset)
print(len(sorted_sequences))

compare_with_original(sorted_sequences, original_subset)
