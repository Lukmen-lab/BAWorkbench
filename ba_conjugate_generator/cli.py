"""Command line interface."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter

from rdkit import RDLogger

from .config import GeneratorConfig
from .cores import get_cores
from .exports import export_outputs
from .generator import generate_candidates
from .migration import migrate_previous_output
from .partners import load_partners


def generate_cmd(args: argparse.Namespace) -> None:
    RDLogger.DisableLog("rdApp.*")
    config = GeneratorConfig(mode=args.mode, max_saccharide_dp=args.max_saccharide_dp, allow_dp6_to_dp10=True)
    cores = get_cores(args.cores)
    partners = load_partners(args.partners)
    candidates, invalid, excluded = generate_candidates(cores, partners, config)
    export_outputs(args.output, candidates, invalid, excluded)
    print(json.dumps({"generated": len(candidates), "invalid": len(invalid), "excluded": len(excluded), "output": args.output}, indent=2))


def summarize_cmd(args: argparse.Namespace) -> None:
    with open(args.input, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    print(json.dumps({
        "candidate_count": len(rows),
        "categories": Counter(row["candidate_category"] for row in rows),
        "reaction_channels": Counter(row["reaction_channel"] for row in rows),
        "manual_review_count": sum(row.get("manual_review_required") == "True" for row in rows),
    }, indent=2))


def migrate_cmd(args: argparse.Namespace) -> None:
    RDLogger.DisableLog("rdApp.*")
    summary = migrate_previous_output(args.input, args.output)
    print(json.dumps(summary, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ba-conjugate-generator")
    sub = parser.add_subparsers(required=True)
    gen = sub.add_parser("generate")
    gen.add_argument("--cores", nargs="+", default=["CA", "CDCA", "DCA", "LCA", "UDCA", "7oxoCA", "7oxoCDCA"])
    gen.add_argument("--partners", default=None)
    gen.add_argument("--output", default="outputs/run_001")
    gen.add_argument("--mode", default="broad_covalent", choices=["primary_docking", "broad_covalent", "manual_review"])
    gen.add_argument("--max-saccharide-dp", type=int, default=5)
    gen.set_defaults(func=generate_cmd)
    summ = sub.add_parser("summarize")
    summ.add_argument("--input", required=True)
    summ.set_defaults(func=summarize_cmd)
    migrate = sub.add_parser("migrate-sites", help="Audit and repair site assignments in a previous candidate export")
    migrate.add_argument("--input", required=True)
    migrate.add_argument("--output", required=True)
    migrate.set_defaults(func=migrate_cmd)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
