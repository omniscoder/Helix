from helix.schema import SPEC_VERSION, manifest


def test_schema_manifest_stable():
    data = manifest()
    assert data["spec_version"] == SPEC_VERSION
    schemas = data["schemas"]
    assert "viz_alignment_ribbon" in schemas
    assert schemas["viz_alignment_ribbon"]["spec_version"] == SPEC_VERSION
