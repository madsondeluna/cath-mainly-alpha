#!/usr/bin/env python3
"""
taxonomy_survey.py

Levantamento de taxons representados nos PDBs do conjunto filtrado (~8.2k estruturas).
Parseia headers PDB locais (registros SOURCE) para extrair ORGANISM_TAXID, ORGANISM_SCIENTIFIC.
Para PDBs antigos sem ORGANISM_TAXID, tenta anotar via NCBI Entrez API (busca por nome).

Saidas:
  data_exploration/taxonomy_survey.csv      -- registro por (pdb_id, tax_id, organismo)
  data_exploration/taxonomy_summary.csv     -- contagem de PDBs por taxon
  data_exploration/taxonomy_no_taxid.csv    -- PDBs que ainda ficaram sem taxid apos anotacao
"""

import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Configuracao
# ---------------------------------------------------------------------------
STRUCTURES_DIR = Path("structures")
OUT_DIR = Path("data_exploration")

OUT_SURVEY = OUT_DIR / "taxonomy_survey.csv"
OUT_SUMMARY = OUT_DIR / "taxonomy_summary.csv"
OUT_NO_TAXID = OUT_DIR / "taxonomy_no_taxid.csv"

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_SLEEP = 0.35   # max ~3 req/s sem API key


# ---------------------------------------------------------------------------
# Parser de SOURCE records PDB
# ---------------------------------------------------------------------------
def parse_source_records(pdb_path: Path) -> list[dict]:
    """
    Le os registros SOURCE de um arquivo PDB.
    Retorna lista de dicts por MOL_ID, com ou sem tax_id.
    """
    source_text = []
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if line.startswith("SOURCE"):
                    source_text.append(line[10:].rstrip())
                elif source_text:
                    break
    except (OSError, UnicodeDecodeError):
        return []

    if not source_text:
        return []

    full_text = " ".join(source_text)
    mol_blocks = re.split(r"MOL_ID\s*:", full_text)

    results = []
    seen_taxa: set[int] = set()

    for block in mol_blocks[1:]:
        tax_id_match = re.search(r"ORGANISM_TAXID\s*:\s*(\d+)", block)
        sci_match = re.search(r"ORGANISM_SCIENTIFIC\s*:\s*([^;]+)", block)
        common_match = re.search(r"ORGANISM_COMMON\s*:\s*([^;]+)", block)

        tax_id = int(tax_id_match.group(1)) if tax_id_match else None
        sci_name = sci_match.group(1).strip().title() if sci_match else None
        common_name = common_match.group(1).strip().title() if common_match else None

        # deduplica por tax_id (quando disponivel)
        if tax_id is not None:
            if tax_id in seen_taxa:
                continue
            seen_taxa.add(tax_id)

        results.append(
            {
                "ncbi_taxonomy_id": tax_id,
                "scientific_name": sci_name,
                "common_name": common_name,
            }
        )

    return results


# ---------------------------------------------------------------------------
# Anotacao via NCBI Entrez
# ---------------------------------------------------------------------------
def lookup_taxid_ncbi(scientific_name: str, session: requests.Session) -> int | None:
    """Busca o NCBI taxonomy ID pelo nome cientifico via Entrez esearch."""
    if not scientific_name:
        return None
    params = {
        "db": "taxonomy",
        "term": scientific_name,
        "retmode": "json",
        "retmax": 1,
    }
    try:
        resp = session.get(NCBI_ESEARCH, params=params, timeout=15)
        resp.raise_for_status()
        ids = resp.json()["esearchresult"]["idlist"]
        return int(ids[0]) if ids else None
    except Exception:
        return None


def annotate_missing_taxids(
    missing_records: list[dict],
    session: requests.Session,
) -> list[dict]:
    """
    Para cada registro sem tax_id, tenta buscar pelo scientific_name no NCBI.
    Retorna apenas os que conseguiram ser anotados.
    """
    # agrupa por nome para nao repetir chamadas iguais
    name_to_taxid: dict[str, int | None] = {}
    unique_names = {r["scientific_name"] for r in missing_records if r["scientific_name"]}

    print(f"  Consultando NCBI para {len(unique_names)} nomes...")
    for name in tqdm(sorted(unique_names), desc="  NCBI Entrez lookup"):
        name_to_taxid[name] = lookup_taxid_ncbi(name, session)
        time.sleep(NCBI_SLEEP)

    annotated = []
    for rec in missing_records:
        sci = rec["scientific_name"]
        if sci and name_to_taxid.get(sci) is not None:
            rec = rec.copy()
            rec["ncbi_taxonomy_id"] = name_to_taxid[sci]
            annotated.append(rec)

    return annotated


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    if not STRUCTURES_DIR.exists():
        sys.exit(f"Diretorio nao encontrado: {STRUCTURES_DIR}")
    OUT_DIR.mkdir(exist_ok=True)

    pdb_files = sorted(
        f for f in STRUCTURES_DIR.glob("*.pdb") if not f.name.startswith("._")
    )
    print(f"Arquivos PDB locais: {len(pdb_files):,}")

    records_with_taxid: list[dict] = []
    records_without_taxid: list[dict] = []    # tem nome mas sem taxid
    no_info: list[str] = []                   # sem nome nem taxid

    for pdb_path in tqdm(pdb_files, desc="Parseando PDB headers"):
        pdb_id = pdb_path.stem.upper()
        organisms = parse_source_records(pdb_path)

        if not organisms:
            no_info.append(pdb_id)
            continue

        for org in organisms:
            entry = {"pdb_id": pdb_id, **org}
            if org["ncbi_taxonomy_id"] is not None:
                records_with_taxid.append(entry)
            elif org["scientific_name"]:
                records_without_taxid.append(entry)
            else:
                no_info.append(pdb_id)

    n_with = len({r["pdb_id"] for r in records_with_taxid})
    n_without = len({r["pdb_id"] for r in records_without_taxid})
    print(f"\nCom ORGANISM_TAXID:           {n_with:,} PDBs")
    print(f"Sem taxid, com nome:          {n_without:,} PDBs  -> anotando via NCBI...")
    print(f"Sem taxid e sem nome:         {len(no_info):,} PDBs  (ignorados)")

    # ------------------------------------------------------------------
    # Anotacao via NCBI Entrez para os sem taxid
    # ------------------------------------------------------------------
    annotated = []
    if records_without_taxid:
        session = requests.Session()
        session.headers.update({"User-Agent": "taxonomy_survey/1.0 (research)"})
        annotated = annotate_missing_taxids(records_without_taxid, session)
        print(f"  Anotados com sucesso: {len({r['pdb_id'] for r in annotated})} PDBs")

    # ------------------------------------------------------------------
    # Combina todos os registros
    # ------------------------------------------------------------------
    all_records = records_with_taxid + annotated
    df = pd.DataFrame(all_records)
    df["source"] = "PDB_header"
    df.loc[df["pdb_id"].isin({r["pdb_id"] for r in annotated}), "source"] = "NCBI_lookup"
    df.to_csv(OUT_SURVEY, index=False)
    print(f"\nRegistros totais: {len(df):,}  ->  {OUT_SURVEY}")

    # ------------------------------------------------------------------
    # Sumario por taxon
    # ------------------------------------------------------------------
    pdb_counts = df.groupby("ncbi_taxonomy_id")["pdb_id"].nunique().rename("n_pdbs")
    sci_names = (
        df.groupby(["ncbi_taxonomy_id", "scientific_name"])
        .size()
        .reset_index(name="freq")
        .sort_values("freq", ascending=False)
        .drop_duplicates("ncbi_taxonomy_id")
        .set_index("ncbi_taxonomy_id")["scientific_name"]
    )
    tax_summary = (
        pdb_counts.to_frame()
        .join(sci_names)
        .reset_index()[["ncbi_taxonomy_id", "scientific_name", "n_pdbs"]]
        .sort_values("n_pdbs", ascending=False)
    )
    tax_summary.to_csv(OUT_SUMMARY, index=False)

    # ------------------------------------------------------------------
    # PDBs que ainda ficaram sem taxid
    # ------------------------------------------------------------------
    still_missing_ids = {r["pdb_id"] for r in records_without_taxid} - {r["pdb_id"] for r in annotated}
    if still_missing_ids or no_info:
        missing_df = pd.DataFrame(
            [r for r in records_without_taxid if r["pdb_id"] in still_missing_ids]
            + [{"pdb_id": p, "scientific_name": None} for p in no_info]
        )
        missing_df.to_csv(OUT_NO_TAXID, index=False)

    # ------------------------------------------------------------------
    # Sumario final
    # ------------------------------------------------------------------
    print(f"\nTaxons unicos (NCBI): {len(tax_summary):,}  ->  {OUT_SUMMARY}")
    print(f"\n=== Top 30 Taxons (PDBs unicos) ===")
    print(tax_summary.head(30).to_string(index=False))
    print(f"\n=== Estatisticas gerais ===")
    print(f"  PDBs com taxonomia: {df['pdb_id'].nunique():,}")
    print(f"  Taxa unicos (NCBI): {len(tax_summary):,}")
    print(f"  Sem taxid (final):  {len(still_missing_ids) + len(set(no_info)):,}")


if __name__ == "__main__":
    main()
