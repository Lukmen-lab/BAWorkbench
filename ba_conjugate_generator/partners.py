"""Partner registry and normalization."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rdkit import Chem

from .validation import mol_from_smiles


INORGANIC_OR_SALT_LABELS = {
    "potassium phosphate",
    "magnesium sulfate",
    "manganese sulfate",
    "ferrous sulfate",
    "ammonia",
    "hydrogen sulfide",
    "iron ii phosphate",
    "iron iii phosphate",
    "manganese phosphate",
}


@dataclass
class Partner:
    partner_name: str
    normalized_name: str
    synonyms: list[str] = field(default_factory=list)
    smiles: str | None = None
    formula: str | None = None
    source_category: str = "user_supplied"
    functional_group_tags: list[str] = field(default_factory=list)
    allowed_reaction_channels: list[str] = field(default_factory=list)
    preferred_reaction_channels: list[str] = field(default_factory=list)
    plausibility_class: str = "chemically_plausible"
    warning_flags: list[str] = field(default_factory=list)
    notes: str = ""
    original_labels: list[str] = field(default_factory=list)

    @property
    def mol(self) -> Chem.Mol | None:
        return mol_from_smiles(self.smiles, self.partner_name)

    @property
    def is_excluded_context_only(self) -> bool:
        return "exclude_context_only" in self.warning_flags


def _p(
    name: str,
    norm: str,
    smiles: str | None,
    tags: list[str],
    channels: list[str],
    category: str,
    *,
    synonyms: list[str] | None = None,
    plausibility: str = "chemically_plausible",
    warnings: list[str] | None = None,
    notes: str = "",
    labels: list[str] | None = None,
) -> Partner:
    return Partner(
        partner_name=name,
        normalized_name=norm,
        synonyms=synonyms or [],
        smiles=smiles,
        source_category=category,
        functional_group_tags=tags,
        allowed_reaction_channels=channels,
        preferred_reaction_channels=channels[:2],
        plausibility_class=plausibility,
        warning_flags=warnings or [],
        notes=notes,
        original_labels=labels or [name],
    )


def built_in_partners() -> list[Partner]:
    a = "C24_amide"
    e = "C24_ester"
    ac = "BA_O_acyl"
    gly = "BA_O_glycoside"
    react = "reactive_carbonyl"
    return [
        _p("D-glucose", "glucose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", ["sugar", "alcohol"], [e, gly], "saccharide", synonyms=["Glucose"], labels=["Glucose", "D-Glucose"], notes="Anomeric stereochemistry is model-level for glycosides."),
        _p("Galactose", "galactose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", ["sugar", "alcohol"], [e, gly], "saccharide"),
        _p("Ribose", "ribose", "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O", ["sugar", "alcohol"], [e, gly], "saccharide"),
        _p("Lactose", "lactose", "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O", ["disaccharide", "alcohol"], [e, gly], "saccharide", warnings=["large_flexible_warning"]),
        _p("Maltose", "maltose", "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O", ["disaccharide", "alcohol"], [e, gly], "saccharide", warnings=["large_flexible_warning"]),
        _p("Sucrose", "sucrose", "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O", ["disaccharide", "alcohol"], [e, gly], "saccharide", warnings=["large_flexible_warning"]),
        _p("Xylulose-5-phosphate", "xylulose-5-phosphate", "O=P(O)(O)OC[C@@H](O)C(=O)CO", ["phosphate", "sugar", "alcohol"], [e, gly], "saccharide", warnings=["manual_review_required"]),
        _p("High-molecular-weight glucans / Exopolysaccharides", "glucan_exopolysaccharide", None, ["polymer", "sugar"], ["glucan_fragments"], "saccharide", warnings=["polymeric_fragment_model"], notes="Undefined polymer represented only as defined DP fragments."),
        _p("Pyruvate", "pyruvic acid", "CC(=O)C(=O)O", ["carboxylic_acid", "keto"], [ac, react], "organic_acid"),
        _p("Acetyl-phosphate", "acetyl phosphate", "CC(=O)OP(=O)(O)O", ["acyl_phosphate", "carboxylic_acid_like"], [ac], "organic_acid", warnings=["manual_review_required"]),
        _p("L-lactate", "lactic acid", "C[C@H](O)C(=O)O", ["carboxylic_acid", "alcohol"], [e, ac], "organic_acid", labels=["L-Lactate"]),
        _p("D-lactate", "d-lactic acid", "C[C@@H](O)C(=O)O", ["carboxylic_acid", "alcohol"], [e, ac], "organic_acid", labels=["D-Lactate"]),
        _p("Acetate", "acetic acid", "CC(=O)O", ["carboxylic_acid"], [ac], "organic_acid", synonyms=["Sodium acetate"], labels=["Acetate", "Sodium Acetate"], notes="Sodium salt context collapsed to acetate chemistry only."),
        _p("Propionate", "propionic acid", "CCC(=O)O", ["carboxylic_acid"], [ac], "organic_acid"),
        _p("3-Hydroxypropionic acid", "3-hydroxypropionic acid", "OCCC(=O)O", ["carboxylic_acid", "alcohol"], [e, ac], "organic_acid", synonyms=["3-HP"]),
        _p("Ethanol", "ethanol", "CCO", ["alcohol"], [e], "alcohol"),
        _p("Glycerol", "glycerol", "OCC(O)CO", ["polyol"], [e], "alcohol"),
        _p("1,2-Propanediol", "1,2-propanediol", "CC(O)CO", ["diol"], [e], "alcohol", synonyms=["1,2-PDO"]),
        _p("1,3-Propanediol", "1,3-propanediol", "OCCCO", ["diol"], [e], "alcohol", synonyms=["1,3-PDO"]),
        _p("Propanol", "propanol", "CCCO", ["alcohol"], [e], "alcohol"),
        _p("3-Hydroxypropionaldehyde", "3-hydroxypropionaldehyde", "OCCC=O", ["aldehyde", "alcohol"], [e, react], "reactive_carbonyl", synonyms=["3-HPA", "Reuterin"], warnings=["reactive_or_artifact_possible", "manual_review_required"]),
        _p("Propionaldehyde", "propionaldehyde", "CCC=O", ["aldehyde"], [react], "reactive_carbonyl", warnings=["reactive_or_artifact_possible", "manual_review_required"]),
        _p("Acrolein", "acrolein", "C=CC=O", ["aldehyde", "michael_acceptor"], [react], "reactive_carbonyl", warnings=["reactive_or_artifact_possible", "manual_review_required"]),
        _p("L-tryptophan", "tryptophan", "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", ["amine", "carboxylic_acid", "indole"], [a, ac], "amino_acid"),
        _p("L-histidine", "histidine", "N[C@@H](Cc1c[nH]cn1)C(=O)O", ["amine", "carboxylic_acid"], [a, ac], "amino_acid"),
        _p("L-glutamic acid", "glutamic acid", "N[C@@H](CCC(=O)O)C(=O)O", ["amine", "carboxylic_acid"], [a, ac], "amino_acid", synonyms=["L-Glutamate"], labels=["L-Glutamate", "L-Glutamic acid"]),
        _p("Histamine", "histamine", "NCCc1c[nH]cn1", ["amine"], [a], "biogenic_amine"),
        _p("GABA", "gaba", "NCCCC(=O)O", ["amine", "carboxylic_acid"], [a, ac], "amino_acid"),
        _p("Indole-3-pyruvate", "indole-3-pyruvic acid", "O=C(O)C(=O)Cc1c[nH]c2ccccc12", ["carboxylic_acid", "indole", "keto"], [ac, react], "indole_metabolite"),
        _p("Indole-3-lactate", "indole-3-lactic acid", "O=C(O)C(O)Cc1c[nH]c2ccccc12", ["carboxylic_acid", "indole", "alcohol"], [e, ac], "indole_metabolite", synonyms=["ILA"]),
        _p("Indole-3-propionate", "indole-3-propionic acid", "O=C(O)CCc1c[nH]c2ccccc12", ["carboxylic_acid", "indole"], [ac], "indole_metabolite", synonyms=["IPA"]),
        _p("Indole-3-acetate", "indole-3-acetic acid", "O=C(O)Cc1c[nH]c2ccccc12", ["carboxylic_acid", "indole"], [ac], "indole_metabolite", synonyms=["I3A"]),
        _p("Glycine", "glycine", "NCC(=O)O", ["amine", "carboxylic_acid"], [a, ac], "amino_acid"),
        _p("Taurine", "taurine", "NCCS(=O)(=O)O", ["amine", "sulfonic_acid"], [a], "amino_sulfonic_acid", synonyms=["2-aminoethanesulfonic acid"], plausibility="reference_known", notes="Generates the physiological C24 taurine amide (tauro-bile-acid) channel."),
        _p("L-phenylalanine", "phenylalanine", "N[C@@H](Cc1ccccc1)C(=O)O", ["amine", "carboxylic_acid"], [a, ac], "amino_acid"),
        _p("L-proline", "proline", "O=C(O)[C@@H]1CCCN1", ["secondary_amine", "carboxylic_acid"], [a, ac], "amino_acid"),
        _p("L-lysine", "lysine", "NCCCC[C@H](N)C(=O)O", ["amine", "carboxylic_acid"], [a, ac], "amino_acid"),
        _p("Beta-alanine", "beta-alanine", "NCCC(=O)O", ["amine", "carboxylic_acid"], [a, ac], "amino_acid", labels=["Beta-Alanine", "Beta-Alanine (as artifact)"], warnings=["artifact_label_present"]),
        _p("Biotin", "biotin", "O=C(O)CCCCC1SCC2NC(=O)N12", ["carboxylic_acid", "vitamin"], [ac], "vitamin"),
        _p("Pantothenic acid", "pantothenic acid", "CC(C)(CO)C(O)C(=O)NCCC(=O)O", ["carboxylic_acid", "alcohol", "vitamin"], [e, ac], "vitamin", synonyms=["Calcium pantothenate"], labels=["Calcium pantothenate"], notes="Calcium salt context collapsed to pantothenic acid chemistry only."),
        _p("Folic acid", "folic acid", "Nc1nc2ncc(CNc3ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc3)nc2c(=O)[nH]1", ["carboxylic_acid", "amine", "vitamin"], [a, ac], "vitamin", warnings=["large_flexible_warning"]),
        _p("Nicotinic acid", "nicotinic acid", "O=C(O)c1cccnc1", ["carboxylic_acid", "vitamin"], [ac], "vitamin", synonyms=["Niacinamide / Nicotinic acid"]),
        _p("Niacinamide", "niacinamide", "NC(=O)c1cccnc1", ["amide", "vitamin"], [], "vitamin", synonyms=["Niacinamide / Nicotinic acid"], warnings=["no_covalent_channel_in_prototype"]),
        _p("pABA", "para-aminobenzoic acid", "Nc1ccc(C(=O)O)cc1", ["amine", "carboxylic_acid"], [a, ac], "vitamin_related", synonyms=["Para-aminobenzoic acid"]),
        _p("Pyridoxine", "pyridoxine", "Cc1ncc(CO)c(CO)c1O", ["alcohol", "phenol", "vitamin"], [e], "vitamin", synonyms=["Pyridoxine HCl"]),
        _p("Pyridoxal", "pyridoxal", "Cc1ncc(C=O)c(CO)c1O", ["aldehyde", "alcohol", "phenol", "vitamin"], [e, react], "vitamin", synonyms=["Pyridoxal"], warnings=["reactive_or_artifact_possible"]),
        _p("Riboflavin", "riboflavin", "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)CO)c2cc1C", ["alcohol", "vitamin"], [e], "vitamin", warnings=["large_flexible_warning"]),
        _p("Thiamine", "thiamine", "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1", ["alcohol", "vitamin"], [e], "vitamin", synonyms=["Thiamine HCl"], notes="Chloride salt context removed."),
        _p("Vitamin K1", "phylloquinone", "CC1=C(C(=O)c2ccccc2C1=O)CCCC(C)CCCC(C)CCCC(C)C", ["quinone", "vitamin"], [react], "vitamin", warnings=["large_flexible_warning", "manual_review_required"], synonyms=["Phylloquinone"]),
        _p("Adenine", "adenine", "Nc1ncnc2[nH]cnc12", ["nucleobase", "amine"], [a], "nucleobase", warnings=["manual_review_required"]),
        _p("Guanine", "guanine", "Nc1nc2[nH]cnc2c(=O)[nH]1", ["nucleobase", "amine"], [a], "nucleobase", warnings=["manual_review_required"]),
        _p("Thymidine", "thymidine", "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O", ["nucleoside", "alcohol"], [e], "nucleobase"),
        _p("Uracil", "uracil", "O=c1cc[nH]c(=O)[nH]1", ["nucleobase"], [], "nucleobase", warnings=["no_covalent_channel_in_prototype"]),
        _p("Xanthine", "xanthine", "O=c1[nH]c(=O)c2[nH]cnc2[nH]1", ["nucleobase"], [], "nucleobase", labels=["Xanthine", "Xanthine (as artifact)"], warnings=["artifact_label_present", "no_covalent_channel_in_prototype"]),
        _p("5-HMF", "5-hydroxymethylfurfural", "O=Cc1oc(CO)cc1", ["aldehyde", "alcohol"], [e, react], "reactive_carbonyl", synonyms=["5-Hydroxymethylfurfural"], warnings=["reactive_or_artifact_possible", "manual_review_required"]),
        _p("Furfural", "furfural", "O=Cc1ccco1", ["aldehyde"], [react], "reactive_carbonyl", warnings=["reactive_or_artifact_possible", "manual_review_required"]),
        _p("Levulinic acid", "levulinic acid", "CC(=O)CCC(=O)O", ["carboxylic_acid", "keto"], [ac, react], "organic_acid"),
        _p("Formic acid", "formic acid", "O=CO", ["carboxylic_acid"], [ac], "organic_acid"),
        _p("Pyroglutamic acid", "pyroglutamic acid", "O=C(O)[C@@H]1CCC(=O)N1", ["carboxylic_acid", "lactam"], [ac], "organic_acid", synonyms=["5-oxoproline"]),
        _p("Indole", "indole", "c1ccc2[nH]ccc2c1", ["indole"], [], "indole_metabolite", warnings=["no_covalent_channel_in_prototype"]),
        _p("Anthranilic acid", "anthranilic acid", "Nc1ccccc1C(=O)O", ["amine", "carboxylic_acid"], [a, ac], "organic_acid"),
        _p("Tyramine", "tyramine", "NCCc1ccc(O)cc1", ["amine", "phenol"], [a, e], "biogenic_amine"),
        _p("4-Amino-5-aminomethyl-2-methylpyrimidine", "amino-methylpyrimidine", "Cc1ncc(CN)c(N)n1", ["amine"], [a], "vitamin_fragment"),
        _p("4-Methyl-5-hydroxyethylthiazole", "hydroxyethylthiazole", "Cc1ncc(CCO)s1", ["alcohol"], [e], "vitamin_fragment"),
        _p("6-formylpterin", "6-formylpterin", "O=Cc1cnc2nc(N)[nH]c(=O)c2n1", ["aldehyde", "pterine"], [react], "vitamin_fragment", synonyms=["Pterines (e.g., 6-formylpterin)"], warnings=["reactive_or_artifact_possible", "manual_review_required"]),
        _p("p-Aminobenzoylglutamic acid", "p-aminobenzoylglutamic acid", "Nc1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1", ["amine", "carboxylic_acid"], [a, ac], "vitamin_fragment"),
        _p("2,3-Diketogulonic acid", "2,3-diketogulonic acid", "O=C(O)C(=O)C(=O)C(O)C(O)CO", ["carboxylic_acid", "keto", "alcohol"], [e, ac, react], "organic_acid", warnings=["reactive_or_artifact_possible"]),
        _p("Threonic acid", "threonic acid", "O=C(O)C(O)C(O)CO", ["carboxylic_acid", "alcohol"], [e, ac], "organic_acid"),
        _p("Pantoic acid", "pantoic acid", "CC(C)(CO)C(O)C(=O)O", ["carboxylic_acid", "alcohol"], [e, ac], "vitamin_fragment"),
        _p("4-pyridoxic acid", "4-pyridoxic acid", "Cc1ncc(C(=O)O)c(CO)c1O", ["carboxylic_acid", "alcohol", "phenol"], [e, ac], "vitamin_fragment"),
        _p("Hypoxanthine", "hypoxanthine", "O=c1[nH]cnc2[nH]cnc12", ["nucleobase"], [], "nucleobase", warnings=["no_covalent_channel_in_prototype"]),
        _p("Cobalamin", "cobalamin", None, ["vitamin", "metal_complex"], [], "vitamin", warnings=["exclude_context_only", "metal_complex_context", "large_flexible_warning"], notes="Cobalamin is too complex/metal-containing for covalent BA prototype generation."),
        _p("L-Ascorbic acid", "ascorbic acid", "O=C1OC(C(O)=C1O)C(O)CO", ["alcohol", "vitamin"], [e], "vitamin"),
        _p("Melanoidins", "melanoidins", None, ["undefined_polymer"], [], "context_only", warnings=["exclude_context_only"], notes="Undefined mixture/polymer; no exact molecule generated."),
    ] + amino_acid_partners() + excluded_context_partners()


def amino_acid_partners() -> list[Partner]:
    rows = [
        ("L-Alanine", "alanine", "C[C@H](N)C(=O)O"),
        ("L-Arginine", "arginine", "N[C@@H](CCCNC(=N)N)C(=O)O"),
        ("L-Asparagine", "asparagine", "N[C@@H](CC(N)=O)C(=O)O"),
        ("L-Aspartic acid", "aspartic acid", "N[C@@H](CC(=O)O)C(=O)O"),
        ("L-Cysteine", "cysteine", "N[C@@H](CS)C(=O)O"),
        ("L-Glutamine", "glutamine", "N[C@@H](CCC(N)=O)C(=O)O"),
        ("L-Isoleucine", "isoleucine", "CC[C@H](C)[C@H](N)C(=O)O"),
        ("L-Leucine", "leucine", "CC(C)C[C@H](N)C(=O)O"),
        ("L-Methionine", "methionine", "CSCC[C@H](N)C(=O)O"),
        ("L-Serine", "serine", "N[C@@H](CO)C(=O)O"),
        ("L-Threonine", "threonine", "C[C@H](O)[C@H](N)C(=O)O"),
        ("L-Tyrosine", "tyrosine", "N[C@@H](Cc1ccc(O)cc1)C(=O)O"),
        ("L-Valine", "valine", "CC(C)[C@H](N)C(=O)O"),
    ]
    return [_p(n, norm, smi, ["amine", "carboxylic_acid"], ["C24_amide", "BA_O_acyl"], "amino_acid") for n, norm, smi in rows]


def excluded_context_partners() -> list[Partner]:
    names = [
        "Potassium Phosphate (K2HPO4 / KH2PO4)",
        "Magnesium Sulfate (MgSO4)",
        "Manganese Sulfate (MnSO4)",
        "Ferrous Sulfate (FeSO4)",
        "Ammonia (NH3)",
        "Hydrogen sulfide (H2S)",
        "Iron (II) Phosphate [Fe3(PO4)2]",
        "Iron (III) Phosphate [FePO4]",
        "Manganese Phosphate [Mn3(PO4)2]",
    ]
    return [
        _p(n, n.lower().split("(")[0].strip().replace("[fe3", "").replace("[fepo4]", ""), None, ["inorganic_context"], [], "inorganic_context", warnings=["exclude_context_only"], notes="Excluded: inorganic salt, ion-pair, gas, or metal-complex context only.")
        for n in names
    ]


def load_partners(path: str | Path | None = None) -> list[Partner]:
    if path is None:
        return built_in_partners()
    text = Path(path).read_text(encoding="utf-8")
    try:
        raw: list[dict[str, Any]] = json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml  # type: ignore

            raw = yaml.safe_load(text)
        except Exception as exc:
            raise RuntimeError("Partner file must be JSON-compatible YAML or PyYAML must be installed") from exc
    return [Partner(**item) for item in raw]
