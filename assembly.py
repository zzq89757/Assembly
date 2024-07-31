

# cut sequences to kmer by common region length (len(common_region)/2)

# whether two sequence have same kmer in terminal




from collections import defaultdict
from itertools import product

def generate_kmers(sequence, k):
    """Generate k-mers for a given sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]

def build_de_bruijn_graph(sequences, k):
    """Build De Bruijn graph from k-mers."""
    edges = defaultdict(list)
    nodes = set()
    
    for sequence in sequences:
        for kmer in generate_kmers(sequence, k):
            prefix = kmer[:-1]
            suffix = kmer[1:]
            edges[prefix].append(suffix)
            nodes.add(prefix)
            nodes.add(suffix)
    
    return edges, nodes

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def extend_path(path, edges, visited, max_mismatch=1, max_indel=1):
    """Extend the path in the De Bruijn graph allowing mismatches and indels."""
    current = path[-1]
    extended = False
    
    for neighbor in edges[current]:
        if neighbor in visited:
            continue
        mismatches = hamming_distance(current[1:], neighbor[:-1])
        if mismatches <= max_mismatch:
            visited.add(neighbor)
            new_path = path + [neighbor]
            yield new_path
            yield from extend_path(new_path, edges, visited, max_mismatch, max_indel)
            visited.remove(neighbor)
            extended = True

    # Handle insertions
    for ins in product("ACGT", repeat=max_indel):
        inserted = current[1:] + ''.join(ins)
        for neighbor in edges[inserted[:len(current)]]:
            if neighbor in visited:
                continue
            mismatches = hamming_distance(inserted, neighbor[:-1])
            if mismatches <= max_mismatch:
                visited.add(neighbor)
                new_path = path + [neighbor]
                yield new_path
                yield from extend_path(new_path, edges, visited, max_mismatch, max_indel)
                visited.remove(neighbor)
                extended = True

    # Handle deletions
    for i in range(1, max_indel + 1):
        deleted = current[1:-i]
        for neighbor in edges[deleted[:len(current)]]:
            if neighbor in visited:
                continue
            mismatches = hamming_distance(deleted, neighbor[:-1])
            if mismatches <= max_mismatch:
                visited.add(neighbor)
                new_path = path + [neighbor]
                yield new_path
                yield from extend_path(new_path, edges, visited, max_mismatch, max_indel)
                visited.remove(neighbor)
                extended = True

    if not extended:
        yield path

def assemble_sequence(sequences, k, max_mismatch=1, max_indel=1):
    edges, nodes = build_de_bruijn_graph(sequences, k)
    assembled_sequences = []
    
    for start_node in nodes:
        for path in extend_path([start_node], edges, set(), max_mismatch, max_indel):
            sequence = path[0]
            for node in path[1:]:
                sequence += node[-1]
            assembled_sequences.append(sequence)
    
    return assembled_sequences

# 示例用法
sequences = [
    "ACGTACGTGACG",
    "ACGTGACGTACG",
    "GACGTGACGTAC"
]
k = 8
assembled_sequences = assemble_sequence(sequences, k, max_mismatch=1, max_indel=1)
print("Assembled Sequences:")
for seq in assembled_sequences:
    print(seq)

