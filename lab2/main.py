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

s = input("Give me a list of genes, in upper format:\n")

s_new = [s[i:i+3] for i in range(0, len(s), 3)]
output = []
for gene in s_new:
    temp = CODON_TO_AMINO_ACID[gene]
    if temp == 'Stop':
        output.append(temp)
        break
    output.append(temp)

print(output)
