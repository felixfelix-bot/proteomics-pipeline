"""
Minimal sanity tests for the proteomics pipeline structure.

These tests verify that required R scripts and config files exist.
They don't require R to be installed — just check file existence.
"""
import os
from pathlib import Path


REPO_ROOT = Path(__file__).parent.parent


def test_r_scripts_exist():
    """All core R analysis scripts should be present."""
    scripts = [
        "R/01_config.R",
        "R/08_targeted_volcanos.R",
        "R/09_flagip_volcano.R",
        "R/10_targeted_venns.R",
        "R/11_targeted_go.R",
        "R/13_chx_crac_analysis.R",
        "R/17_crac_string_network.R",
        "R/18_venn_overflow_examples.R",
        "R/utils.R",
        "R/run_step.R",
    ]
    for s in scripts:
        assert (REPO_ROOT / s).exists(), f"Missing: {s}"


def test_makefile_has_targets():
    """Makefile should contain key targets."""
    makefile = (REPO_ROOT / "Makefile").read_text()
    for target in ["flagip-volcano", "crac-network", "venn-overflow-examples"]:
        assert target in makefile, f"Missing Makefile target: {target}"


def test_config_has_color_palette():
    """Config should define GLOBAL_COLORS."""
    config = (REPO_ROOT / "R/01_config.R").read_text()
    assert "GLOBAL_COLORS" in config
    assert "flag_both" in config
