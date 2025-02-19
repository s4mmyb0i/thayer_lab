import networkx as nx

def calculate_path_counts(graph):
    path_counts = {3: 0, 4: 0, 5: 0}

    for length in range(3, 6):
        total_paths = 0
        for target in graph.nodes:
            if target != "1":  # Skip the source node
                paths = nx.all_simple_paths(graph, source="1", target=target, cutoff=length)
                total_paths += sum(1 for _ in paths)
        path_counts[length] = total_paths

    return path_counts

G = nx.read_graphml("/home/samvit/Documents/2024/project/all_atom/project_0.8/graphml/P/060.graphml")
print(calculate_path_counts(G))
