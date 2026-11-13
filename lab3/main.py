# design a software app that is able to process a fasta file of a genome either from a virus or a bacteria.
# find the frequencies of the codons over the genome and plot and show into the console the top 10 most frequent codons.
# use the AI in order to formulate a prompt in which the frequencies are shown to the AI and in which the AI is asked 
# about a list of foods that contain the least amount of the aminoacids associated with the codons found in the top 10 list detected by it.
# note: use a sliding window with 3 positions in order to measure the frequency of the codons.

#search on google NCBI genome download FTP

# hello sir, i am google support, you need to give me $100 to unlock your account

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



with open('GCF_023122875.1_ASM2312287v1_genomic.fna', 'r') as f:
    # skip FASTA header lines (starting with '>') and join the sequence lines
    s = ''.join(line.strip() for line in f if not line.startswith('>')).upper()

# DNA -> RNA, so the codons match CODON_TO_AMINO_ACID
s_new = s.replace('T', 'U')

frequency = {}

# len(s_new) - 2 so the last window is always a full codon
for i in range(0, len(s_new) - 2, 3):
    codon = s_new[i:i+3]
    frequency[codon] = frequency.get(codon, 0) + 1

print("Codon frequencies:")
sorted_frequency = sorted(frequency.items(), key=lambda x: x[1], reverse=True)
top10 = sorted_frequency[:10]
for codon, freq in top10:
    print(f"{codon} ({CODON_TO_AMINO_ACID.get(codon, '?')}): {freq}")

import matplotlib.pyplot as plt

plt.bar([codon for codon, freq in top10], [freq for codon, freq in top10])
plt.xlabel('Codons')
plt.ylabel('Frequency')
plt.title('Top 10 Most Frequent Codons')
plt.savefig('codon_frequencies.png')


# ---------------------------------------------------------------------------
# Foods with the least amount of the amino acids from the top 10 codons,
# using the USDA FoodData Central REST API (SR Legacy dataset).
# ---------------------------------------------------------------------------
import json
import os
import urllib.parse
import urllib.request

FDC_URL = "https://api.nal.usda.gov/fdc/v1/foods/list"
FDC_CACHE = "fdc_sr_legacy.json"

# amino acid -> nutrient name used by FoodData Central.
# Gln and Asn are not measured separately: they are reported inside Glu and Asp.
AMINO_ACID_TO_FDC_NUTRIENT = {
    "Phe": "Phenylalanine",
    "Leu": "Leucine",
    "Ile": "Isoleucine",
    "Met": "Methionine",
    "Val": "Valine",
    "Ser": "Serine",
    "Pro": "Proline",
    "Thr": "Threonine",
    "Ala": "Alanine",
    "Tyr": "Tyrosine",
    "His": "Histidine",
    "Gln": "Glutamic acid",
    "Asn": "Aspartic acid",
    "Lys": "Lysine",
    "Asp": "Aspartic acid",
    "Glu": "Glutamic acid",
    "Cys": "Cystine",
    "Trp": "Tryptophan",
    "Arg": "Arginine",
    "Gly": "Glycine",
}


def load_api_key(path=".env"):
    with open(path) as f:
        for line in f:
            key, _, value = line.strip().partition("=")
            if key == "API_KEY":
                return value.strip().strip('"').strip("'")
    raise RuntimeError("API_KEY not found in .env")


def fetch_all_foods(api_key):
    """Download every SR Legacy food (with nutrients), cached in FDC_CACHE."""
    if os.path.exists(FDC_CACHE):
        with open(FDC_CACHE) as f:
            return json.load(f)

    foods = []
    page = 1
    while True:
        params = urllib.parse.urlencode({
            "api_key": api_key,
            "dataType": "SR Legacy",
            "pageSize": 200,
            "pageNumber": page,
        })
        with urllib.request.urlopen(f"{FDC_URL}?{params}") as resp:
            batch = json.load(resp)
        if not batch:
            break
        foods.extend(batch)
        print(f"Downloaded page {page} ({len(foods)} foods)")
        page += 1

    with open(FDC_CACHE, "w") as f:
        json.dump(foods, f)
    return foods


# how many of the top 10 codons code for each amino acid (stop codons dropped)
aa_weight = {}
for codon, freq in top10:
    aa = CODON_TO_AMINO_ACID.get(codon)
    if aa and aa != "Stop":
        aa_weight[aa] = aa_weight.get(aa, 0) + freq

# amino acids sharing a nutrient (Gln/Glu, Asn/Asp) are added together
nutrient_weight = {}
for aa, weight in aa_weight.items():
    nutrient = AMINO_ACID_TO_FDC_NUTRIENT[aa]
    nutrient_weight[nutrient] = nutrient_weight.get(nutrient, 0) + weight

total_weight = sum(nutrient_weight.values())
print("\nAmino acids from the top 10 codons:", ", ".join(aa_weight))

foods = fetch_all_foods(load_api_key())

# score = grams of the target amino acids per 100 g, weighted by how often
# their codons appear in the top 10
ranked = []
for food in foods:
    amounts = {n["name"]: n.get("amount", 0) for n in food["foodNutrients"]}
    if all(nutrient in amounts for nutrient in nutrient_weight):
        score = sum(amounts[n] * w for n, w in nutrient_weight.items()) / total_weight
    elif amounts.get("Protein") == 0:
        # no amino acid rows, but no protein at all (oils, sugar, ...)
        score = 0.0
    else:
        # amino acids not reported: we can't rank this food
        continue
    ranked.append((score, food["description"], amounts))

ranked.sort(key=lambda x: x[0])

print(f"\nFoods with the least of these amino acids ({len(ranked)} foods ranked):")
header = "".join(f"{n[:6]:>8}" for n in nutrient_weight)
print(f"{'Food':<55}{'score':>8}{header}")
for score, name, amounts in ranked[:20]:
    cols = "".join(f"{amounts.get(n, 0):>8.3f}" for n in nutrient_weight)
    print(f"{name[:54]:<55}{score:>8.4f}{cols}")
print("(amounts in g per 100 g of food)")

plt.show()
