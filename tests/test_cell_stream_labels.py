import pytest

from malva_client.client import _map_stream_gene_names


@pytest.mark.parametrize(
    ("positions", "features", "requested", "expected"),
    [
        # Current server: positions within the requested panel.
        ([0, 1, 0], ["CD19", "COL1A1"], ["CD19", "COL1A1"], ["CD19", "COL1A1", "CD19"]),
        # Older server: store ordinals, one per requested gene, in panel order.
        ([5394], ["CD19"], ["CD19"], ["CD19"]),
        ([13283, 5394], ["CD19", "COL1A1"], ["CD19", "COL1A1"], ["CD19", "COL1A1"]),
        # Unresolvable indices stay explicit instead of silently wrong.
        ([5394, 9000], [], ["CD19"], ["feature_5394", "feature_9000"]),
    ],
)
def test_map_stream_gene_names(positions, features, requested, expected):
    assert _map_stream_gene_names(positions, features, requested) == expected
