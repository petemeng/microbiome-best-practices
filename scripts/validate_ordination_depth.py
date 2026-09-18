#!/usr/bin/env python3
"""Check chapter 22's scientific contracts against a completed reader run."""
import argparse
import csv
import json
from pathlib import Path

from PIL import Image
from reader_reproducibility import reader_blocks, script_text


def validate(root, run):
    checks = []
    def check(name, passed):
        checks.append({"check": name, "passed": bool(passed)})
    def table(name):
        with (run / "results/22-ordination" / (name + ".tsv")).open() as stream:
            return list(csv.DictReader(stream, delimiter="\t"))
    source = (root / "chapters/22-ordination-unconstrained.qmd").read_text()
    expected = script_text(source, reader_blocks(source))
    check("companion_matches_executed_chunks",
          (root / "examples/22-ordination-unconstrained.R").read_text() == expected)
    check("same_feature_universe_for_clr", 'clr_transform(comm, pseudocount = 0.5)' in source)
    check("explicit_cailliez_branch", 'correction_name <- if (needs_correction) "Cailliez" else "None"' in source)
    check("removed_false_lingoes_claim", '添加 Lingoes 常数' not in source and 'with a Lingoes additive correction' not in source)
    check("no_broken_backslash_math", '\\frac' not in source and '\\mathrm' not in source)
    check("sampling_design_caveat", 'Sampling-unit independence was not established' in source)
    check("confounding_diagnostic", 'model.matrix(~ Group + Saline, metadata)' in source)
    pcoa = table("pcoa_diagnostic")[0]
    check("no_substantive_negative_eigenvalues", int(pcoa["NegativeEigenvalues"]) == 0)
    check("no_correction_applied", pcoa["Correction"] == "None" and float(pcoa["Constant"]) == 0)
    check("axis1_recomputed", abs(float(pcoa["Axis1Percent"]) - 14.24717) < 0.0001)
    check("axis2_recomputed", abs(float(pcoa["Axis2Percent"]) - 11.65554) < 0.0001)
    nmds = table("nmds_diagnostic")
    check("both_nmds_dimensions", [int(r["Dimensions"]) for r in nmds] == [2, 3])
    check("repeated_best_solutions", all(int(r["RepeatedBestSolutions"]) > 0 for r in nmds))
    check("stress_drops_in_three_dimensions", 0 < float(nmds[1]["Stress"]) < float(nmds[0]["Stress"]) < 1)
    sens = table("sensitivity")
    check("five_controlled_comparisons", len(sens) == 5)
    check("feature_counts", [int(r["Features"]) for r in sens] == [11113, 11113, 5700, 5700, 11113])
    check("geometry_changes_more_than_offset", float(sens[-1]["DistanceSpearman"]) < min(float(r["DistanceSpearman"]) for r in sens[:-1]))
    groups = table("group_tests")
    check("both_groupings_shown", [r["Grouping"] for r in groups] == ["Wetland group", "Soil salinity"])
    check("dispersion_interpretation_matches", float(groups[0]["Dispersion_P"]) < .05 <= float(groups[1]["Dispersion_P"]))
    for stem in ("22-pcoa-bray", "22-nmds-bray", "22-distance-fit", "22-clr-pca"):
        for extension in ("pdf", "png", "tiff"):
            path = run / "figures" / f"{stem}.{extension}"
            check(f"export_{stem}_{extension}", path.is_file() and path.stat().st_size > 1000)
        with Image.open(run / "figures" / f"{stem}.png") as image:
            check(f"resolution_{stem}", min(image.size) >= 1000 and min(image.info.get("dpi", (0, 0))) >= 299)
    return {"status": "passed" if all(c["passed"] for c in checks) else "failed",
            "scope": "chapter 22 scientific and export checks; completed cold run required",
            "checks_total": len(checks), "checks_passed": sum(c["passed"] for c in checks), "checks": checks}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--run-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    report = validate(a.project_root, a.run_root)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "checks"}))
    raise SystemExit(report["status"] != "passed")
