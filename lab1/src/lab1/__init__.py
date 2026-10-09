from collections import Counter
S_init = "cccccccccccccccc"
S = S_init.lower()
list_of_chars = ['a','b','c','g']
   
def count_combinations(S):
    pairs = Counter(S[i:i + 2] for i in range(len(S) - 1))
    triples = Counter(S[i:i + 3] for i in range(len(S) - 2))
    print(pairs)
    print(triples)
    return (pairs, triples)

def print_the_counts_to_histograms(pairs, triples):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    for ax, counter, title in zip(axes, (pairs, triples), ("Pairs", "Triples")):
        items = counter.most_common()
        ax.bar([k for k, _ in items], [v for _, v in items])
        ax.set_title(title)
        ax.set_ylabel("count")
    fig.tight_layout()
    plt.show()

def main() -> None:
    res = {
        key: val for key, val in Counter(S).items() if key in list_of_chars
    }
    print("ex1:")
    print(res)
    print("ex2:")
    pairs, triples = count_combinations(S)
    print_the_counts_to_histograms(pairs, triples)
    
    
