from argparse import ArgumentParser
from numpy.typing import NDArray
import numpy as np

def create_partition(mi_matrix: NDArray[np.float64], alpha: float) -> dict[int, list[int]]:
    """
    Creates a partitioning of residues parameterized by a choice of alpha
    (0 <= alpha <= 1) based on their mutual information matrix . Two residues
    are grouped in the same partition if there exists a path of residues such
    that each pair in the path maintians a correlation value greater than or
    equal to alpha.

    Parameters:
        `mi_matrix` - A numpy array of shape (num_residues, num_residues), where each entry is a value between 0 and 1.
        `alpha` - A threshold value which determines the minimum correlation needed between any two adjacent residues in the path
    """
    # Initialize variables
    num_residues = mi_matrix.shape[0]
    ungrouped_residues = list(range(num_residues))
    partitions = dict()

    # Find cluster for each residue
    while len(ungrouped_residues) > 0:
        # Select representative residue for this partition
        representative = ungrouped_residues.pop(0)
        partitions[representative] = [representative]
        unexplored = [representative]

        # Expand community until all nodes within it are explored
        while len(unexplored) > 0:
            # Select an unexplored residue to expand
            selected_residue = unexplored.pop()

            # Find all residues which are highly correlated with the selected residue
            correlated_residues = [resid for resid in ungrouped_residues if mi_matrix[selected_residue][resid] >= alpha]

            # Add these residues to the partition and to the list of residues to explore later
            partitions[representative] += correlated_residues
            unexplored += correlated_residues

            # Mark residues as belonging to some partition
            ungrouped_residues = [resid for resid in ungrouped_residues if resid not in correlated_residues]

    return partitions

def test_function(partition, alpha) -> float:
    """alpha / num partitions"""
    return alpha / len(partition.values())

def modularity(partition, alpha) -> float:
    """
    Reference: Diez, G., Nagel, D., & Stock, G. (2022). Correlation-based
    feature selection to identify functional dynamics in proteins. Journal of
    Chemical Theory and Computation, 18(8), 5079-5088.
    """
    m = sum((len(part) * (len(part) - 1) / 2) for part in partition.values())
    if m == 0:
        return 0
    total = 0
    for part in partition.values():
        n = len(part)
        e_c = sum((data[i][j] for i in part for j in part if i < j))
        k_c = n * (n - 1) / 2
        total += (e_c - (k_c ** 2) / (2 * m))
    return total / (2 * m)

if __name__ == "__main__":
    # Parse cli arguments
    parser = ArgumentParser()

    parser.add_argument(
        "--mifile",
        type=str,
        required=True,
        help="Pairwise mutual information file to be used for partitioning system"
    )

    parser.add_argument(
        "--alpha",
        type=float,
        help="Chosen alpha threshold"
    )

    args = parser.parse_args()

    # Load data
    data = np.loadtxt(args.mifile, dtype=float, delimiter=',')
    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            if data[x][y] < 0:
                data[x][y] = 0

    # If no alpha value is given pick one which optimizes for a specific value
    if args.alpha == None:
        # Find best partition
        alpha_values = tuple(0.001 * i for i in range(1000))
        q = -1

        best_partition = None
        best_alpha = None
        for alpha in alpha_values:
            partition = create_partition(data, alpha)
            curr_q = test_function(partition, alpha)

            if curr_q > q:
                q = curr_q
                best_partition = partition
                best_alpha = alpha

        # Print clusters
        print("clusters:")
        for part in best_partition.values():
            print([x + 1 for x in part])
    else:
        # Calculate chosen partition
        partition = create_partition(data, args.alpha)

        # Print clusters
        print("clusters:")
        for part in partition.values():
            print([x + 1 for x in part])
