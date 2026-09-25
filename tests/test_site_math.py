"""Equation markup is compiled locally and contains valid native MathML."""
import importlib.util
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('site_math', ROOT / 'scripts/site_math.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_model_equations_compile():
    html = module.typeset((ROOT / 'web/template.html').read_text())
    maths = re.findall(r'<math\b.*?</math>', html, re.S)
    assert len(maths) >= 50
    for source in maths:
        tree = ET.fromstring(source)
        assert not list(tree.iter('{http://www.w3.org/1998/Math/MathML}merror'))
    assert r'\[' not in html
    assert r'\(' not in html
    assert '<script src=' not in html


def test_registry_requires_matching_equation_count():
    import pytest
    with pytest.raises(ValueError):
        module.typeset('<div class="eq">extra</div>')
