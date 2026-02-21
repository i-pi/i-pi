import json
from pathlib import Path

import pytest

from ipi.utils.io.inputs.io_json import json_parse_string, xmlnode_to_json
from ipi.utils.io.inputs.io_xml import xml_parse_file


REGTESTS_ROOT = (
    Path(__file__).resolve().parents[3] / "regression_tests" / "tests"
)
REGTEST_INPUTS = sorted(REGTESTS_ROOT.rglob("input.xml"))

if not REGTEST_INPUTS:
    raise RuntimeError(f"No regression input.xml files found under {REGTESTS_ROOT}")


def _norm_text(text):
    return " ".join(str(text).split())


def _normalize_xml_tree(node):
    fields = []
    for key, value in node.fields:
        if hasattr(value, "fields"):
            fields.append((key, _normalize_xml_tree(value)))
        else:
            normalized = _norm_text(value)
            # Ignore indentation/newline-only text nodes.
            if key == "_text" and normalized == "":
                continue
            fields.append((key, normalized))

    return (
        node.name,
        tuple(sorted(node.attribs.items())),
        tuple(fields),
    )


@pytest.mark.parametrize(
    "xml_path",
    REGTEST_INPUTS,
    ids=[str(path.relative_to(REGTESTS_ROOT)) for path in REGTEST_INPUTS],
)
def test_regtest_input_xml_json_roundtrip_semantic_equivalence(xml_path):
    with open(xml_path) as xml_file:
        original = xml_parse_file(xml_file)

    as_json = xmlnode_to_json(original)
    roundtrip = json_parse_string(json.dumps(as_json))

    assert _normalize_xml_tree(roundtrip) == _normalize_xml_tree(original)
