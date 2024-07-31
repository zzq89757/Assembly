import numpy as np

class OverlapGraph:
    def __init__(self, sequences, min_overlap_length):
        self.sequences = sequences
        self.min_overlap_length = min_overlap_length
        self.graph = self.build_overlap_graph(sequences)

    def calculate_overlap(self, a, b):
        """Calculate the optimal overlap between two sequences considering mismatches and indels."""
        m = len(a)
        n = len(b)

        # Create a scoring matrix
        dp = np.zeros((m + 1, n + 1), dtype=int)
        
        # Initialization
        for i in range(m + 1):
            dp[i][0] = -i  # Cost of indels in 'a'
        for j in range(n + 1):
            dp[0][j] = -j  # Cost of indels in 'b'

        # Fill the scoring matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = dp[i-1][j-1] + (1 if a[i-1] == b[j-1] else -1)  # Match or mismatch
                delete = dp[i-1][j] - 1  # Deletion
                insert = dp[i][j-1] - 1  # Insertion
                dp[i][j] = max(match, delete, insert)

        # Find the maximum overlap
        max_overlap = 0
        overlap_start = 0
        for j in range(n):
            if dp[m][j] > max_overlap:
                max_overlap = dp[m][j]
                overlap_start = j

        return max_overlap, overlap_start

    def build_overlap_graph(self, sequences):
        """Build a graph where edges represent suffix/prefix overlaps considering mismatches and indels."""
        graph = {}
        for a in sequences:
            for b in sequences:
                if a != b:
                    overlap_length, _ = self.calculate_overlap(a, b)
                    if overlap_length >= self.min_overlap_length:
                        graph[(a, b)] = overlap_length
        return graph

    def find_longest_path(self):
        """Find the longest path in the overlap graph."""
        from collections import defaultdict

        dist = defaultdict(int)
        path = defaultdict(list)
        
        for (a, b), length in self.graph.items():
            if dist[a] + length > dist[b]:
                dist[b] = dist[a] + length
                path[b] = path[a] + [b]
        
        end_node = max(dist, key=dist.get)
        return path[end_node]

def assemble_sequences(sequences, min_overlap_length):
    """Assemble sequences using overlap graph and longest path."""
    og = OverlapGraph(sequences, min_overlap_length)
    path = og.find_longest_path()

    if not path:
        return ""

    consensus = path[0]
    for i in range(1, len(path)):
        overlap_length, _ = og.calculate_overlap(path[i-1], path[i])
        consensus += path[i][overlap_length:]

    return consensus

# Example usage:
sequences = [
    "AGTACGCA",
    "TACGCAGA",
    "GCAGAGAT",
    "AGATCGGA",
    # TACGCAGAAGAGAT
]

min_overlap_length = 1
consensus_sequence = assemble_sequences(sequences, min_overlap_length)
print("Consensus sequence:", consensus_sequence)
