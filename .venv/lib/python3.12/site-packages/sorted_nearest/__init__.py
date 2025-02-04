import pkg_resources

from sorted_nearest.src.annotate_clusters import annotate_clusters  # NOQA: F401
from sorted_nearest.src.cluster_by import cluster_by  # NOQA: F401
from sorted_nearest.src.clusters import find_clusters  # NOQA: F401
from sorted_nearest.src.introns import find_introns  # NOQA: F401
from sorted_nearest.src.k_nearest import k_nearest_next_nonoverlapping, k_nearest_previous_nonoverlapping  # NOQA: F401
from sorted_nearest.src.k_nearest_ties import get_all_ties, get_different_ties  # NOQA: F401
from sorted_nearest.src.max_disjoint_intervals import max_disjoint  # NOQA: F401
from sorted_nearest.src.merge_by import merge_by  # NOQA: F401
from sorted_nearest.src.sorted_nearest import (  # NOQA: F401
    nearest_next_nonoverlapping,
    nearest_nonoverlapping,
    nearest_previous_nonoverlapping,
)
from sorted_nearest.src.tiles import maketiles  # NOQA: F401
from sorted_nearest.src.windows import makewindows  # NOQA: F401

__version__ = pkg_resources.get_distribution("sorted_nearest").version
