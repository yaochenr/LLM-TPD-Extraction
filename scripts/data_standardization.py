"""Data standardization pipeline for the extracted dataset."""

import argparse
import re
import time
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests

HUMAN = 9606

DEFAULT_INPUT_CSV = "data/protacs/llm/protac-data-llm.csv"
DEFAULT_OUTPUT_CSV = "data/AutoTPD-PROTAC.csv"


def _place_after(df: pd.DataFrame, anchor: str, followers: list[str]) -> pd.DataFrame:
    followers = [c for c in followers if c in df.columns]
    if anchor not in df.columns or not followers:
        return df
    cols = [c for c in df.columns if c not in followers]
    idx = cols.index(anchor)
    new_cols = cols[: idx + 1] + followers + cols[idx + 1:]
    return df[new_cols]


# ===========================================================================
# Cell line standardization (Cellosaurus)
# ===========================================================================

def cellosaurus_search(query: str) -> dict | None:
    url = "https://api.cellosaurus.org/search/cell-line"
    try:
        r = requests.get(url, params={
            "q": f"id:{query}",
            "format": "json",
            "fields": "ac,id,sy,ox",
        }, timeout=15)
        if not r.ok:
            return None
        results = r.json().get("Cellosaurus", {}).get("cell-line-list", [])
        if not results:
            return None
        for entry in results:
            names = entry.get("name-list", [])
            rec_name = None
            for n in names:
                if n.get("type") == "identifier":
                    rec_name = n["value"]
                    break
            if rec_name and rec_name.lower() == query.lower():
                acc = None
                for a in entry.get("accession-list", []):
                    if a.get("type") == "primary":
                        acc = a["value"]
                        break
                species = None
                taxonomy_id = None
                sp_list = entry.get("species-list", [])
                if sp_list:
                    species = sp_list[0].get("label")
                    taxonomy_id = int(sp_list[0].get("accession", 0)) or None
                return {"accession": acc, "name": rec_name, "species": species, "taxonomy_id": taxonomy_id}
    except Exception as e:
        print(f"  Cellosaurus error for '{query}': {e}")
    return None


def cellosaurus_fetch(accession: str) -> dict | None:
    url = f"https://api.cellosaurus.org/cell-line/{accession}"
    try:
        r = requests.get(url, params={"format": "json", "fields": "ac,id,sy,ox"}, timeout=15)
        if not r.ok:
            return None
        entry = r.json().get("Cellosaurus", {}).get("cell-line-list", [{}])[0]
        names = entry.get("name-list", [])
        rec_name = next((n["value"] for n in names if n["type"] == "identifier"), None)
        species = None
        taxonomy_id = None
        sp_list = entry.get("species-list", [])
        if sp_list:
            species = sp_list[0].get("label")
            taxonomy_id = int(sp_list[0].get("accession", 0)) or None
        return {"accession": accession, "name": rec_name, "species": species, "taxonomy_id": taxonomy_id}
    except Exception as e:
        print(f"  Cellosaurus fetch error for '{accession}': {e}")
    return None


CL_SEARCH_MAP = {
    "22RV1":                "22Rv1",
    "293T":                 "HEK293T",
    "A375":                 "A-375",
    "A375 xenograft":       "A-375",
    "A549":                 "A-549",
    "Ba/F3 Bcr-AblT315I":  "Ba/F3 BCR-ABL-T315I",
    "BMDMs":                None,
    "CAL-12-T":             "CAL-12T",
    "CAOV4":                "Caov-4",
    "Flip-In 293":          "Flp-In-293",
    "Flp-In T-REx 293":    "Flp-In-T-REx-293",
    "Flp293T GFP-GSPT1(389-499)/mCherry": "Flp-In-293",
    "G2019S LRRK2 MEFs":   "MEF Lrrk2-p.G2019S KI",
    "H1666":                "NCI-H1666",
    "H293T":                "HEK293T",
    "HCC-15":               "HCC15",
    "HCC-364 vr1":          "HCC364",
    "HCT116":               "HCT 116",
    "HCT15":                "HCT 15",
    "HEK-293T (BTK-nLuc)":  "HEK293T",
    "HEK-293T (IKZF1-nLuc)":"HEK293T",
    "HEK-293T (IKZF3-nLuc)":"HEK293T",
    "HEK293 (α1A-AR stably transfected)": "HEK293",
    "HEK293 BRD9-HiBiT/CAS9/FF":  "HEK293",
    "HEK293 BRD9-HiBiT/FF/CAS9":  "HEK293",
    "HEK293-hTau":          "HEK293",
    "HEK293/BRD9-HiBiT":   "HEK293",
    "HepG2":                "Hep-G2",
    "Hs746T":               "Hs 746.T",
    "Insig-silenced HepG2": "Hep-G2",
    "Jurkat E6-1":          "Jurkat E6.1",
    "K562":                 "K-562",
    "K562 HiBiT-tagged":    "K-562",
    "MKN-45":               "MKN45",
    "MM.1S":                "MM1.S",
    "MOLT4":                "MOLT-4",
    "MOLT4 HiBiT-IKZF1":   "MOLT-4",
    "MV-4-11":              "MV4-11",
    "MV4;11":               "MV4-11",
    "NIH-3T3":              "NIH 3T3",
    "Panc Tu-I":            "PancTu-I",
    "R1441C LRRK2 MEFs":   "MEF Lrrk2-p.R1441C KI",
    "Rama 37 (S100A4-transfected)": "Rama 37",
    "SK-MEL-239 C4":        "SK-MEL-239",
    "SK-N-BE(2)-C":         "CVCL_0529",  
    "TMD-8":                "TMD8",
    "TMD8 BTK-GFP/mCh":    "TMD8",
    "U87":                  "U-87MG ATCC",
    "WI38 non-senescent (NCs)": "WI-38",
    "WT LRRK2 MEFs":        "Wt MEFs",
    # --- Molecular glue cell lines ---
    "BxPC3":                "BxPC-3",
    "CWR22Rv1 (CRISPR HiBiT knock-in)": "22Rv1",
    "HEK293T–Cas9":         "HEK293T/Cas9 AAVS1",
    "KBM7":                 "KBM-7",
    "Karpas 422":           "Karpas-422",
    "MOLT4 CK1α-HiBiT":    "MOLT-4",
    "primary human Tregs":  None,
    # --- MolGlueDB baseline cell lines ---
    "BE(2)-C":              "CVCL_0529",  
    "Human PBMCs":          None,
    "Human monocytes (n=3)": None,
    "Human whole blood":    None,
    "JAR":                  "JAR",
    "Mouse primary splenocytes": None,
    "OCI-Ly1":              "OCI-Ly1",
    "OCI-Ly10":             "OCI-Ly10",
    "OCI-Ly11":             "OCI-Ly11",
    "OCILY3":               "OCI-Ly3",
    "PA-1":                 "PA-1",
    "SU-DHL-4":             "SU-DHL-4",
    "hMDM":                 None,
    "primary human erythroblasts": None,
    # --- PROTACDB baseline cell lines ---
    "22Rv1":                "22Rv1",
    "231MFP":               "MDA-MB-231",
    "293FT":                "293FT",
    "4T1":                  "4T1",
    "A-204":                "A-204",
    "A152T neurons":        None,
    "A2780":                "A2780",
    "A431":                 "A-431",
    "A549/Taxol":           "A549-Taxol",
    "AML12":                "AML12",
    "AsPC-1":               "AsPC-1",
    "BT474":                "BT-474",
    "Ba/F3":                "Ba/F3",
    "BaF3 FLT3-ITD":        "Ba/F3",
    "Bel-7402":             "BEL-7402",
    "CA-46":                "CA-46",
    "CAL-12T":              "CAL-12T",
    "CAL-148":              "CAL-148",
    "CAL33":                "CAL 33",
    "CCLP1":                "CC-LP-1",
    "CHL-1":                "CHL-1",
    "COLO 205":             "COLO 205",
    "COLO-205":             "COLO 205",
    "CWR22Rv1":             "22Rv1",
    "Caki-1":               "Caki-1",
    "Calu-1":               "Calu-1",
    "Calu-6":               "Calu-6",
    "Capan-1":              "Capan-1",
    "DB":                   "DB",
    "DOHH2":                "DOHH-2",
    "Dermal fibroblast":    None,
    "EBC-1":                "EBC-1",
    "EOL-1":                "EoL-1",
    "ER-positive breast cancer cells": None,
    "GBM43":                "GBM43",
    "Germ cells":           None,
    "H1650":                "NCI-H1650",
    "H1650R":               "NCI-H1650",
    "H1975":                "NCI-H1975",
    "H23":                  "NCI-H23",
    "H3122":                "NCI-H3122",
    "H3255":                "NCI-H3255",
    "H358":                 "NCI-H358",
    "H838":                 "NCI-H838",
    "HAP1":                 "HAP1",
    "HBL-1":                "HBL-1",
    "HCC-827":              "HCC827",
    "HCC1937":              "HCC1937",
    "HCC827":               "HCC827",
    "HCT-116":              "HCT 116",
    "HD-MB03":              "HD-MB03",
    "HEI-OC1":              "HEI-OC1",
    "HEK293":               "HEK293",
    "HEK293/BRD7-HiBiT":    "HEK293",
    "HEK293A":              "HEK293-A",
    "HEK293T":              "HEK293T",
    "HEL":                  "HEL",
    "HGC-27":               "HGC-27",
    "HGC-27/HUVEC":         "HGC-27",
    "HL-60":                "HL-60",
    "HL60":                 "HL-60",
    "HLE":                  "HLE",
    "HLF":                  "HLF",
    "HT-1080":              "HT-1080",
    "HT-29":                "HT-29",
    "HT1080":               "HT-1080",
    "HT29":                 "HT-29",
    "HUCCT1":               "HuCCT1",
    "HUH-1":                "HuH-1",
    "HeLa":                 "HeLa",
    "Hep3B2.1-7":           "Hep 3B2.1-7",
    "Hs578t":               "Hs 578T",
    "HuH-7":                "HuH-7",
    "Huh7.5":               "Huh-7.5",
    "IgE MM":               None,
    "JeKo-1":               "JeKo-1",
    "Jeko-1":               "JeKo-1",
    "Jurkat":               "Jurkat",
    "KATO III":             "KATO III",
    "KCL22":                "KCL-22",
    "KM12":                 "KM12",
    "KOPN49":               "KOPN-49",
    "KOPT-K1":              "KOPT-K1",
    "KU182":                "KU-182",
    "KU812":                "KU812",
    "KU812 CML":            "KU812",
    "KYSE-270":             "KYSE-270",
    "KYSE520":              "KYSE-520",
    "Karpas 299":           "Karpas-299",
    "Karpas299":            "Karpas-299",
    "Kasumi-1":             "Kasumi-1",
    "Kelly":                "Kelly",
    "L-O2":                 "L-02",
    "LCC2":                 "LCC2",
    "LNCaP":                "LNCaP",
    "LnCaP95":              "LNCaP-95",
    "MCF-7":                "MCF-7",
    "MCF7":                 "MCF-7",
    "MDA-MB-231":           "MDA-MB-231",
    "MDA-MB-436":           "MDA-MB-436",
    "MDA-MB-453":           "MDA-MB-453",
    "MDA-MB-468":           "MDA-MB-468",
    "MDA-Pca-2b":           "MDA-PCa-2b",
    "MEFs":                 "MEF",
    "MHH-CALL-4":           "MHH-CALL-4",
    "MIA PaCa-2":           "MIA PaCa-2",
    "MOLM-13":              "MOLM-13",
    "MOLM-14":              "MOLM-14",
    "MOLM-16":              "MOLM-16",
    "MOLM13":               "MOLM-13",
    "MOLT-4":               "MOLT-4",
    "MPro-eGFP stable cells": None,
    "MV4-11":               "MV4-11",
    "Mino":                 "Mino",
    "Molm-16":              "MOLM-16",
    "Molt-4":               "MOLT-4",
    "Mouse 4935":           None,
    "NAMALWA":              "Namalwa",
    "NB4":                  "NB4",
    "NCI-H1568":            "NCI-H1568",
    "NCI-H1703":            "NCI-H1703",
    "NCI-H1975":            "NCI-H1975",
    "NCI-H2030":            "NCI-H2030",
    "NCI-H2228":            "NCI-H2228",
    "NCI-H23":              "NCI-H23",
    "NCI-H358":             "NCI-H358",
    "NCI-H661":             "NCI-H661",
    "NCI-H838":             "NCI-H838",
    "NGP":                  "NGP",
    "Namalwa":              "Namalwa",
    "OCI-AML2":             "OCI-AML2",
    "OCI-AML3":             "OCI-AML3",
    "OVCAR-5":              "OVCAR-5",
    "OVCAR8":               "OVCAR-8",
    "PA1":                  "PA-1",
    "PANC1":                "PANC-1",
    "PBMC":                 None,
    "PC-3":                 "PC-3",
    "PC3":                  "PC-3",
    "PC3-S1":               "PC-3",
    "PC9":                  "PC-9",
    "PDX SJBALL020589":     None,
    "Pfeiffer":             "Pfeiffer",
    "Primary Cardiomyocytes": None,
    "Primary Sertoli cells": None,
    "RAMOS":                "Ramos",
    "RAW 264.7":            "RAW 264.7",
    "RI-1":                 "RI-1",
    "RKO":                  "RKO",
    "RS4;11":               "RS4;11",
    "Raji":                 "Raji",
    "Ramos":                "Ramos",
    "SH-SY5Y":              "SH-SY5Y",
    "SK-Hep-1":             "SK-HEP-1",
    "SK-LU-1":              "SK-LU-1",
    "SK-MEL-246":           "SK-MEL-246",
    "SK-MEL-28":            "SK-MEL-28",
    "SK-MEL-5":             "SK-MEL-5",
    "SKNO1":                "SKNO-1",
    "SNU-387":              "SNU-387",
    "SNU-398":              "SNU-398",
    "SNU-423":              "SNU-423",
    "SR":                   "SR",
    "SRD15":                "SRD-15",
    "SU-DHL-1":             "SU-DHL-1",
    "SUM149":               "SUM149PT",
    "SUP-M2":               "SUP-M2",
    "SW1573":               "SW1573",
    "SW480":                "SW480",
    "SW620":                "SW620",
    "Sk-Mel-28":            "SK-MEL-28",
    "Snca OE-PFF seeding HEK293T": "HEK293T",
    "T-cell leukemia Jurkat": "Jurkat",
    "T47D":                 "T-47D",
    "THP-1":                "THP-1",
    "THP1":                 "THP-1",
    "TM3":                  "TM3",
    "TMD8":                 "TMD8",
    "U251":                 "U-251MG",
    "U937":                 "U-937",
    "VCaP":                 "VCaP",
    "VCaP AR+":             "VCaP",
    "WI38":                 "WI-38",
    "XLA":                  None,
    "Z138":                 "Z-138",
    "hPBMC":                None,
    "human dermal papilla": None,
    "platelets":            None,
    # Primary cells / non-standard — no Cellosaurus entry
    "4xSNCA iPSC-derived dopaminergic neurons (DANs)": None,
    "EGFR-AXL chimeric cell line":   None,
    "EGFR-MERTK chimeric cell line": None,
    "EGFR-TYRO3 chimeric cell line": None,
    "FTD19-L5-RC6 (tau-A152T)":      None,
    "MGH2046-RC1 (tau-P301L)":       None,
    "PBMCs":                 None,
    "Q175 mouse cortical neurons":   None,
    "human PBMCs":           None,
    "mouse BMDMs (naïve)":   None,
    "mouse primary Sertoli cells (C57BL/6N, 6 dpp)": None,
    "mouse primary germ cells (C57BL/6N, 6 dpp)":    None,
    "murine peritoneal macrophages":  None,
    "primary CLL cells":     None,
    "primary rat neonatal cardiomyocytes": None,
    "tau-P301L iPSC-derived neurons": None,
}


def standardize_cell_lines(df_llm: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cl_list = sorted(
        df_llm["Cell_Line"]
        .dropna()
        .astype(str)
        .str.strip()
        .unique()
    )
    print(f"Unique Cell_Line values: {len(cl_list)}")

    cl_info: dict[str, dict] = {}
    not_found = []

    for i, cl in enumerate(cl_list):
        print(f"[{i+1}/{len(cl_list)}] {cl}", end="")
        search_name = CL_SEARCH_MAP.get(cl, cl)
        if search_name is None:
            not_found.append(cl)
            print(f"  SKIPPED (no standard cell line)")
            continue
        if search_name.startswith("CVCL_"):
            result = cellosaurus_fetch(search_name)
        else:
            result = cellosaurus_search(search_name)
        if result:
            cl_info[cl] = result
            sp = result.get("species", "")
            if result["name"] != cl:
                print(f"  -> {result['name']} ({result['accession']}, {sp})")
            else:
                print(f"  OK ({result['accession']}, {sp})")
        else:
            not_found.append(cl)
            print(f"  NOT FOUND")
        time.sleep(0.15)

    print(f"\n--- Results ---")
    print(f"Found: {len(cl_info)}/{len(cl_list)}")
    print(f"Not found: {len(not_found)}")
    if not_found:
        print(f"\nNot found list:")
        for cl in not_found:
            print(f"  {cl}")

    def _cl_field(cl, field):
        if pd.notna(cl):
            return cl_info.get(cl, {}).get(field)
        return None

    _cl_stripped = df_llm["Cell_Line"].astype("string").str.strip()

    df_llm["Cell_Line_ID"] = _cl_stripped.map(lambda cl: _cl_field(cl, "accession"))
    df_llm["Cell_Line_Species"] = _cl_stripped.map(lambda cl: _cl_field(cl, "species"))
    df_llm["_taxonomy_id"] = _cl_stripped.map(
        lambda cl: _cl_field(cl, "taxonomy_id") or HUMAN
    ).astype(int)

    n_id = df_llm["Cell_Line_ID"].notna().sum()
    n_sp = df_llm["Cell_Line_Species"].notna().sum()
    print(f"\nCell_Line_ID mapped: {n_id}/{len(df_llm)} rows")
    print(f"Cell_Line_Species mapped: {n_sp}/{len(df_llm)} rows")

    df_llm = _place_after(df_llm, "Cell_Line", ["Cell_Line_ID", "Cell_Line_Species"])

    return df_llm, cl_info


def build_taxid_to_species(cl_info: dict) -> dict:
    taxid_to_species = {HUMAN: "Homo sapiens (Human)"}
    for info in cl_info.values():
        tid = info.get("taxonomy_id")
        sp = info.get("species")
        if tid and sp:
            taxid_to_species[tid] = sp
    return taxid_to_species


# ===========================================================================
# Gene / UniProt / mutation lookup helpers
# ===========================================================================

def parse_mutation(mut_str: str | None) -> tuple[str, int, str] | None:
    if not mut_str:
        return None
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', mut_str)
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None


def hgnc_fetch(symbol: str) -> str | None:
    try:
        r = requests.get(
            f"https://rest.genenames.org/fetch/symbol/{quote(symbol)}",
            headers={"Accept": "application/json"}, timeout=10,
        )
        if r.ok:
            docs = r.json().get("response", {}).get("docs", [])
            if docs:
                return docs[0]["symbol"]
    except:
        pass
    return None


def hgnc_search(query: str) -> str | None:
    try:
        r = requests.get(
            f"https://rest.genenames.org/search/{quote(query)}",
            headers={"Accept": "application/json"}, timeout=10,
        )
        if r.ok:
            docs = r.json().get("response", {}).get("docs", [])
            if docs:
                return docs[0]["symbol"]
    except:
        pass
    return None


def uniprot_gene_search(query: str, organism_id: int = HUMAN) -> str | None:
    url = "https://rest.uniprot.org/uniprotkb/search"
    for q in [
        f"(gene:{query}) AND (organism_id:{organism_id}) AND (reviewed:true)",
        f'(protein_name:"{query}") AND (organism_id:{organism_id}) AND (reviewed:true)',
    ]:
        try:
            r = requests.get(url, params={
                "query": q, "format": "tsv",
                "fields": "accession,gene_primary", "size": 1,
            }, headers={"User-Agent": "protac-normalize/1.0"}, timeout=10)
            if r.ok:
                lines = r.text.strip().split("\n")
                if len(lines) > 1:
                    gene = lines[1].split("\t")[1].strip()
                    if gene:
                        return gene
        except:
            pass
    return None


def find_gene_symbol(search_term: str, organism_id: int = HUMAN) -> tuple[str | None, str]:
    if organism_id == HUMAN:
        gene = hgnc_fetch(search_term)
        if gene:
            return gene, f"HGNC fetch '{search_term}'"
        time.sleep(0.1)
        gene = hgnc_search(search_term)
        if gene:
            return gene, f"HGNC search '{search_term}'"
        time.sleep(0.1)

    gene = uniprot_gene_search(search_term, organism_id)
    if gene:
        return gene, f"UniProt({organism_id}) '{search_term}'"

    return None, "NOT FOUND"


def fetch_uniprot_accession_and_seq(gene_symbol: str, organism_id: int = HUMAN) -> tuple[str | None, str | None]:
    if not gene_symbol:
        return None, None
    url = "https://rest.uniprot.org/uniprotkb/search"
    query = f"(gene_exact:{gene_symbol}) AND (organism_id:{organism_id}) AND (reviewed:true)"
    try:
        r = requests.get(url, params={
            "query": query, "format": "tsv",
            "fields": "accession,gene_primary,sequence", "size": 1,
        }, headers={"User-Agent": "protac-normalize/1.0"}, timeout=15)
        if r.ok:
            lines = r.text.strip().split("\n")
            if len(lines) > 1:
                parts = lines[1].split("\t")
                accession = parts[0].strip() if len(parts) > 0 else None
                sequence = parts[2].strip() if len(parts) > 2 else None
                return accession, sequence
    except Exception as e:
        print(f"  UniProt accession+seq error for {gene_symbol}: {e}")
    return None, None


def fetch_uniprot_record_by_accession(accession: str) -> dict | None:
    accession = str(accession).strip()
    if not accession:
        return None
    url = f"https://rest.uniprot.org/uniprotkb/{quote(accession)}.json"
    try:
        r = requests.get(url, headers={"User-Agent": "protac-normalize/1.0"}, timeout=15)
        if not r.ok:
            return None
        data = r.json()
        genes = data.get("genes") or []
        gene = None
        if genes:
            gene = (genes[0].get("geneName") or {}).get("value")
        sequence = (data.get("sequence") or {}).get("value")
        return {
            "gene": gene,
            "uniprot": data.get("primaryAccession") or accession,
            "sequence": sequence,
        }
    except Exception as e:
        print(f"  UniProt accession fetch error for {accession}: {e}")
    return None


_variation_cache: dict[str, list] = {}


def fetch_variation_features(accession: str) -> list:
    if accession in _variation_cache:
        return _variation_cache[accession]
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{accession}"
    try:
        r = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
        if r.ok:
            features = r.json().get("features", [])
            _variation_cache[accession] = features
            return features
    except Exception as e:
        print(f"  Variation API error for {accession}: {e}")
    _variation_cache[accession] = []
    return []


def find_var_id(accession: str, wt_aa: str, position: int, mut_aa: str) -> str | None:
    features = fetch_variation_features(accession)
    for f in features:
        if (str(f.get("begin", "")) == str(position)
                and f.get("wildType") == wt_aa
                and f.get("alternativeSequence") == mut_aa):
            return f.get("ftId")
    return None


def clinvar_check_variant(gene: str, mutation_str: str) -> str | None:
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        r = requests.get(url, params={
            "db": "clinvar",
            "term": f"{gene}[gene] AND {mutation_str}",
            "retmode": "json",
        }, timeout=10)
        if r.ok:
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                return ids[0]
    except:
        pass
    return None


def apply_mutation(sequence: str, wt_aa: str, position: int, mut_aa: str) -> str | None:
    idx = position - 1
    if idx < 0 or idx >= len(sequence):
        print(f"  apply_mutation: position {position} out of range (seq length {len(sequence)})")
        return None
    if sequence[idx] != wt_aa:
        print(f"  apply_mutation: expected {wt_aa} at position {position}, found {sequence[idx]}")
        return None
    return sequence[:idx] + mut_aa + sequence[idx + 1:]


# ===========================================================================
# POI cleaning + enrichment
# ===========================================================================

poi_alias_map = {
    "BCL-XL": "Bcl-xL",
    "BCL-xL": "Bcl-xL",
    "Bcl-xL": "Bcl-xL",
    "Bcl-xl": "Bcl-xL",
    "Brd2": "BRD2",
    "Brd3": "BRD3",
    "Brd4": "BRD4",
}

POI_LOOKUP = {
    "AXL":        {"search": "AXL",      "mutation": None},
    "BCL-2":      {"search": "BCL2",     "mutation": None},
    "BRD2":       {"search": "BRD2",     "mutation": None},
    "BRD3":       {"search": "BRD3",     "mutation": None},
    "BRD4":       {"search": "BRD4",     "mutation": None},
    "BRD4BD1":    {"search": "BRD4",     "mutation": None},
    "BRD4BD2":    {"search": "BRD4",     "mutation": None},
    "BRD7":       {"search": "BRD7",     "mutation": None},
    "BRD9":       {"search": "BRD9",     "mutation": None},
    "BTK":        {"search": "BTK",      "mutation": None},
    "CDK4":       {"search": "CDK4",     "mutation": None},
    "CDK6":       {"search": "CDK6",     "mutation": None},
    "CDK12":      {"search": "CDK12",    "mutation": None},
    "CDK13":      {"search": "CDK13",    "mutation": None},
    "CRBN":       {"search": "CRBN",     "mutation": None},
    "FER":        {"search": "FER",      "mutation": None},
    "GDI2":       {"search": "GDI2",     "mutation": None},
    "GSPT1":      {"search": "GSPT1",    "mutation": None},
    "HDAC1":      {"search": "HDAC1",    "mutation": None},
    "HDAC2":      {"search": "HDAC2",    "mutation": None},
    "HDAC3":      {"search": "HDAC3",    "mutation": None},
    "HDAC4":      {"search": "HDAC4",    "mutation": None},
    "HDAC6":      {"search": "HDAC6",    "mutation": None},
    "HDAC8":      {"search": "HDAC8",    "mutation": None},
    "HMGCR":      {"search": "HMGCR",    "mutation": None},
    "IDO1":       {"search": "IDO1",     "mutation": None},
    "IKZF1":      {"search": "IKZF1",    "mutation": None},
    "IKZF3":      {"search": "IKZF3",    "mutation": None},
    "KEAP1":      {"search": "KEAP1",    "mutation": None},
    "LRRK2":      {"search": "LRRK2",    "mutation": None},
    "MDM2":       {"search": "MDM2",     "mutation": None},
    "MERTK":      {"search": "MERTK",    "mutation": None},
    "MET":        {"search": "MET",      "mutation": None},
    "NR4A1":      {"search": "NR4A1",    "mutation": None},
    "PARP1":      {"search": "PARP1",    "mutation": None},
    "PARP2":      {"search": "PARP2",    "mutation": None},
    "PBRM1":      {"search": "PBRM1",    "mutation": None},
    "REXO4":      {"search": "REXO4",    "mutation": None},
    "RIPK2":      {"search": "RIPK2",    "mutation": None},
    "S100A4":     {"search": "S100A4",   "mutation": None},
    "SGK3":       {"search": "SGK3",     "mutation": None},
    "SLC2A3":     {"search": "SLC2A3",   "mutation": None},
    "SLC38A2":    {"search": "SLC38A2",  "mutation": None},
    "SLC9A1":     {"search": "SLC9A1",   "mutation": None},
    "SMARCA2":    {"search": "SMARCA2",  "mutation": None},
    "SMARCA4":    {"search": "SMARCA4",  "mutation": None},
    "TYRO3":      {"search": "TYRO3",    "mutation": None},
    "WDR5":       {"search": "WDR5",     "mutation": None},
    "XIAP":       {"search": "XIAP",     "mutation": None},

    "AURORA-A":                            {"search": "AURKA",    "mutation": None},
    "Androgen receptor (AR)":              {"search": "Androgen receptor",       "mutation": None},
    "Bcl-xL":                              {"search": "BCL2L1",   "mutation": None},
    "Brd4 long":                           {"search": "BRD4",     "mutation": None},
    "Brd4 short":                          {"search": "BRD4",     "mutation": None},
    "Cdc20":                               {"search": "CDC20",    "mutation": None},
    "Cyclophilin A (CypA)":                {"search": "PPIA",     "mutation": None},
    "EGFRWT":                              {"search": "EGFR",     "mutation": None},
    "ENL":                                 {"search": "MLLT1",    "mutation": None},
    "ENL (WT and mutant)":                 {"search": "MLLT1",    "mutation": None},
    "ERK1/2":                              {"search": "MAPK1",    "mutation": None},
    "FAK":                                 {"search": "PTK2",     "mutation": None},
    "FKBP12":                              {"search": "FKBP1A",   "mutation": None},
    "GSK-3β":                              {"search": "GSK3B",    "mutation": None},
    "HSP90α":                              {"search": "HSP90AA1", "mutation": None},
    "HSP90β":                              {"search": "HSP90AB1", "mutation": None},
    "IMPDH":                               {"search": "IMPDH1",   "mutation": None},
    "MIF":                                 {"search": "MIF",      "mutation": None},
    "PARP-1":                              {"search": "PARP1",    "mutation": None},
    "PDEδ":                                {"search": "PDE6D",    "mutation": None},
    "PIKfyve":                             {"search": "PIKFYVE",  "mutation": None},
    "SHP2":                                {"search": "PTPN11",   "mutation": None},
    "Squalene synthase (SQS; FDFT1)":      {"search": "FDFT1",    "mutation": None},
    "TC-PTP":                              {"search": "PTPN2",    "mutation": None},
    "Tau (P-tau S396)":                    {"search": "MAPT",     "mutation": None},
    "Tau (insoluble fraction)":            {"search": "MAPT",     "mutation": None},
    "Tau (total)":                         {"search": "MAPT",     "mutation": None},
    "tau":                                 {"search": "MAPT",     "mutation": None},
    "c-MET":                               {"search": "MET",      "mutation": None},
    "cIAP1":                               {"search": "BIRC2",    "mutation": None},
    "cIAP2":                               {"search": "BIRC3",    "mutation": None},
    "hCAII":                               {"search": "CA2",      "mutation": None},
    "p38 MAPK":                            {"search": "MAPK14",   "mutation": None},
    "p38α":                                {"search": "MAPK14",   "mutation": None},
    "p38δ":                                {"search": "MAPK13",   "mutation": None},
    "pVHL30":                              {"search": "VHL",      "mutation": None},
    "phosphorylated p38 MAPK (p-p38)":     {"search": "MAPK14",   "mutation": None},
    "α-Synuclein (αSyn)":                  {"search": "SNCA",     "mutation": None},
    "α-synuclein aggregates":              {"search": "SNCA",     "mutation": None},
    "α1A-adrenergic receptor (α1A-AR)":    {"search": "ADRA1A", "mutation": None},

    "BRAF G466V":    {"search": "BRAF", "mutation": "G466V"},
    "BRAF G469A":    {"search": "BRAF", "mutation": "G469A"},
    "BRAF V600E":    {"search": "BRAF", "mutation": "V600E"},
    "BRAF-p61V600E": {"search": "BRAF", "mutation": "V600E"},
    "BRAFV600E":     {"search": "BRAF", "mutation": "V600E"},
    "Bcr-AblT315I":  {"search": "ABL1", "mutation": "T315I"},
    "KRASG12C":      {"search": "KRAS", "mutation": "G12C"},
    # --- PROTACDB baseline targets without provided UniProt accessions ---
    "Alpha-tubulin": {"search": None,  "mutation": None},
    "BCR-ABL":       {"search": "ABL1",    "mutation": None},
    "BRD4-L":        {"search": "BRD4",    "mutation": None},
    "Beta3-tubulin": {"search": "TUBB3",   "mutation": None},
    "EML4-ALK":      {"search": "ALK",     "mutation": None},
    "HSP90":         {"search": "HSP90AA1", "mutation": None},
    "NS3":           {"search": None,      "mutation": None},
    "c-Src":         {"search": "SRC",     "mutation": None},
    "p-p38":         {"search": "MAPK14",  "mutation": None},
    # --- PROTACDB baseline mutation targets ---
    "ALK G1202R":             {"search": "ALK",   "mutation": "G1202R"},
    "AR L702H and T878A":     {"search": "AR",    "mutation": None},
    "AR T878A":               {"search": "AR",    "mutation": "T878A"},
    "BCR-ABL T315I":          {"search": "ABL1",  "mutation": "T315I"},
    "BTK C481A":              {"search": "BTK",   "mutation": "C481A"},
    "BTK C481G":              {"search": "BTK",   "mutation": "C481G"},
    "BTK C481S":              {"search": "BTK",   "mutation": "C481S"},
    "BTK C481T":              {"search": "BTK",   "mutation": "C481T"},
    "BTK C481W":              {"search": "BTK",   "mutation": "C481W"},
    "EGFR Del19":             {"search": "EGFR",  "mutation": None},
    "EGFR Del19/T790M/C797S": {"search": "EGFR",  "mutation": None},
    "EGFR L858R":             {"search": "EGFR",  "mutation": "L858R"},
    "EGFR L858R/T790M":       {"search": "EGFR",  "mutation": None},
    "EML4-ALK C1156Y":        {"search": "ALK",   "mutation": "C1156Y"},
    "EML4-ALK G1202R":        {"search": "ALK",   "mutation": "G1202R"},
    "EML4-ALK L1196M":        {"search": "ALK",   "mutation": "L1196M"},
    "EML4-ALK L1196M/G1202R": {"search": "ALK",   "mutation": None},
    "JAK2 R683G":             {"search": "JAK2",  "mutation": "R683G"},
    "KRAS G12C":              {"search": "KRAS",  "mutation": "G12C"},
    "KRAS G12D":              {"search": "KRAS",  "mutation": "G12D"},
    "LRRK2 G2019S":           {"search": "LRRK2", "mutation": "G2019S"},

    "HCV NS3": {"search": None, "mutation": None},
    # --- Molecular glue targets ---
    "Androgen receptor (f-AR)":             {"search": "Androgen receptor", "mutation": None},
    "BCL6":                                 {"search": "BCL6",    "mutation": None},
    "CK1α":                                 {"search": "CSNK1A1", "mutation": None},
    "Cyclin K (CCNK)":                      {"search": "CCNK",    "mutation": None},
    "cyclin K":                             {"search": "CCNK",    "mutation": None},
    "cyclin K (CCNK)":                      {"search": "CCNK",    "mutation": None},
    "IKZF2":                                {"search": "IKZF2",   "mutation": None},
    "Lin28":                                {"search": "LIN28A",  "mutation": None},
    "SALL4":                                {"search": "SALL4",   "mutation": None},
    # --- MolGlueDB baseline targets ---
    "BRD4-S":                               {"search": "BRD4",    "mutation": None},
    "CDK12-CCNK":                           {"search": None,      "mutation": None},
    "CDO1":                                 {"search": "CDO1",    "mutation": None},
    "GEMIN3":                               {"search": "DDX20",   "mutation": None},
    "GSPT2":                                {"search": "GSPT2",   "mutation": None},
    "IKZF4":                                {"search": "IKZF4",   "mutation": None},
    "Lin28A":                               {"search": "LIN28A",  "mutation": None},
    "Lin28B":                               {"search": "LIN28B",  "mutation": None},
    "NEK7":                                 {"search": "NEK7",    "mutation": None},
    "PDE6D":                                {"search": "PDE6D",   "mutation": None},
    "PLZF":                                 {"search": "ZBTB16",  "mutation": None},
    "RIOK2":                                {"search": "RIOK2",   "mutation": None},
    "WEE1":                                 {"search": "WEE1",    "mutation": None},
    "WIZ":                                  {"search": "WIZ",     "mutation": None},
    "ZBTB16":                               {"search": "ZBTB16",  "mutation": None},
    "ZNF-90":                               {"search": "ZNF90",   "mutation": None},
    "dALK3":                                {"search": "ALK",     "mutation": None},
    "eGFP(-SD36)":                          {"search": None,      "mutation": None},
    "eGFP(-SD40)":                          {"search": None,      "mutation": None},
    # Reporter / fusion constructs
    "FLuc-S4D (Firefly luciferase fused to SALL4 degron)": {"search": None, "mutation": None},
    "eGFP (fused to SD36)":                 {"search": None, "mutation": None},
    "eGFP (fused to SD40)":                 {"search": None, "mutation": None},
    "eGFP (fused to SD40 F18Q)":            {"search": None, "mutation": None},
    "eGFP (fused to SD56)":                 {"search": None, "mutation": None},

}


def _normalize_uniprot_accession(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    s = str(value).strip()
    if not s or s.lower() in {"nan", "none"}:
        return None
    return re.split(r"[;,\s]+", s)[0].strip() or None


def get_lookup_mutation_for_target(target: str | None) -> str | None:
    if not target or pd.isna(target):
        return None
    lookup = POI_LOOKUP.get(str(target))
    if lookup is None:
        return None
    mutation = lookup.get("mutation")
    return mutation if parse_mutation(mutation) else None


def clean_and_enrich_poi_from_uniprot(df_llm: pd.DataFrame, taxid_to_species: dict) -> pd.DataFrame:
    if "Target_Uniprot" in df_llm.columns:
        if "Degradation_Target_Uniprot" not in df_llm.columns:
            df_llm.rename(columns={"Target_Uniprot": "Degradation_Target_Uniprot"}, inplace=True)
        else:
            df_llm["Degradation_Target_Uniprot"] = df_llm["Degradation_Target_Uniprot"].combine_first(
                df_llm["Target_Uniprot"]
            )
            df_llm.drop(columns=["Target_Uniprot"], inplace=True)

    df_llm["Degradation_Target_Uniprot"] = (
        df_llm["Degradation_Target_Uniprot"].map(_normalize_uniprot_accession).astype("string")
    )

    pairs = (
        df_llm[["Degradation_Target", "Degradation_Target_Uniprot"]]
        .dropna(subset=["Degradation_Target", "Degradation_Target_Uniprot"])
        .drop_duplicates()
        .values.tolist()
    )
    pairs.sort(key=lambda x: (x[0], x[1]))
    print(f"Unique (target, UniProt accession) pairs: {len(pairs)}")

    record_cache: dict[str, dict | None] = {}
    poi_info: dict[tuple[str, str], dict] = {}

    for i, (poi, accession) in enumerate(pairs):
        print(f"[{i+1}/{len(pairs)}] {poi}  ({accession})")
        info = {"gene": None, "uniprot": accession, "mutation_id": None, "sequence": None, "source": []}

        if accession not in record_cache:
            record_cache[accession] = fetch_uniprot_record_by_accession(accession)
            time.sleep(0.1)

        record = record_cache[accession]
        if not record:
            print(f"  UniProt: accession not found")
            poi_info[(poi, accession)] = info
            continue

        info["gene"] = record["gene"]
        info["uniprot"] = record["uniprot"]
        info["sequence"] = record["sequence"]
        info["source"].append("UniProt")
        print(f"  Gene: {info['gene'] or 'MISSING'}")
        print(f"  UniProt: {info['uniprot']}  (seq length: {len(info['sequence']) if info['sequence'] else 0})")

        mutation_str = get_lookup_mutation_for_target(poi)
        parsed = parse_mutation(mutation_str)
        if parsed:
            wt, pos, mut = parsed
            var_id = find_var_id(info["uniprot"], wt, pos, mut)
            if var_id:
                info["mutation_id"] = var_id
                info["source"].append("UniProt variation")
                print(f"  Mutation {mutation_str}: {var_id}")
            else:
                print(f"  Mutation {mutation_str}: no UniProt VAR ID found")

            if info["sequence"]:
                mutated = apply_mutation(info["sequence"], wt, pos, mut)
                if mutated:
                    info["sequence"] = mutated
                    print(f"  Sequence: applied {mutation_str} at position {pos}")

        poi_info[(poi, accession)] = info

    fallback_pairs = (
        df_llm[df_llm["Degradation_Target_Uniprot"].isna()][["Degradation_Target", "_taxonomy_id"]]
        .dropna(subset=["Degradation_Target"])
        .drop_duplicates()
        .values.tolist()
    )
    fallback_pairs.sort(key=lambda x: (x[0], x[1]))
    print(f"Unique fallback (target, organism) pairs without UniProt accession: {len(fallback_pairs)}")

    fallback_info: dict[tuple[str, int], dict] = {}

    for i, (poi, org_id) in enumerate(fallback_pairs):
        org_id = int(org_id)
        print(f"[fallback {i+1}/{len(fallback_pairs)}] {poi}  ({taxid_to_species.get(org_id, str(org_id))})")
        info = {"gene": None, "uniprot": None, "mutation_id": None, "sequence": None, "source": []}

        lookup = POI_LOOKUP.get(poi)
        if lookup is None:
            print(f"  WARNING: '{poi}' not in POI_LOOKUP, skipping")
            fallback_info[(poi, org_id)] = info
            continue

        search_term = lookup["search"]
        mutation_str = lookup["mutation"]

        if search_term is None:
            print(f"  (no gene)")
            fallback_info[(poi, org_id)] = info
            continue

        gene, gene_src = find_gene_symbol(search_term, org_id)
        info["gene"] = gene
        if gene:
            src_tag = gene_src.split("(")[0].split(" ")[0]
            info["source"].append(src_tag)
            print(f"  Gene: {gene}  [{gene_src}]")
        else:
            print(f"  Gene: MISSING  [{gene_src}]")
            fallback_info[(poi, org_id)] = info
            continue

        accession, sequence = fetch_uniprot_accession_and_seq(gene, org_id)
        info["uniprot"] = accession
        info["sequence"] = sequence
        if accession:
            if "UniProt" not in info["source"]:
                info["source"].append("UniProt")
            print(f"  UniProt: {accession}  (seq length: {len(sequence) if sequence else 0})")
        else:
            print(f"  UniProt: not found for gene {gene} (organism {org_id})")
            fallback_info[(poi, org_id)] = info
            continue

        parsed = parse_mutation(mutation_str)
        if parsed:
            wt, pos, mut = parsed
            if org_id == HUMAN:
                clinvar_id = clinvar_check_variant(gene, mutation_str)
                if clinvar_id:
                    info["source"].append("ClinVar")
                    print(f"  ClinVar: variant confirmed (ID={clinvar_id})")

            var_id = find_var_id(accession, wt, pos, mut)
            if var_id:
                info["mutation_id"] = var_id
                print(f"  Mutation {mutation_str}: {var_id}")
            else:
                print(f"  Mutation {mutation_str}: no UniProt VAR ID found")

            if sequence:
                mutated = apply_mutation(sequence, wt, pos, mut)
                if mutated:
                    info["sequence"] = mutated
                    print(f"  Sequence: applied {mutation_str} at position {pos}")

        fallback_info[(poi, org_id)] = info
        time.sleep(0.1)

    def _get_field(row, field):
        accession = _normalize_uniprot_accession(row.get("Degradation_Target_Uniprot"))
        if accession:
            key = (row["Degradation_Target"], accession)
            return poi_info.get(key, {}).get(field)
        key = (row["Degradation_Target"], int(row["_taxonomy_id"]))
        return fallback_info.get(key, {}).get(field)

    df_llm["Degradation_Target_Gene"] = df_llm.apply(lambda r: _get_field(r, "gene"), axis=1)
    df_llm["Degradation_Target_Uniprot"] = df_llm.apply(lambda r: _get_field(r, "uniprot"), axis=1)
    df_llm["Degradation_Target_Uniprot_MutationID"] = df_llm.apply(lambda r: _get_field(r, "mutation_id"), axis=1)
    df_llm["Degradation_Target_Sequence"] = df_llm.apply(lambda r: _get_field(r, "sequence"), axis=1)

    df_llm.drop(columns=["POI_raw", "POI_Name"], inplace=True, errors="ignore")

    mapping_records = [
        {
            "Degradation_Target": k[0],
            "UniProt": v["uniprot"],
            "Gene": v["gene"],
            "MutationID": v["mutation_id"],
            "Seq_Length": len(v["sequence"]) if v["sequence"] else None,
            "Source": "/".join(v["source"]) if v["source"] else "",
        }
        for k, v in sorted(poi_info.items())
    ]
    mapping_records.extend([
        {
            "Degradation_Target": k[0],
            "UniProt": v["uniprot"],
            "Gene": v["gene"],
            "MutationID": v["mutation_id"],
            "Seq_Length": len(v["sequence"]) if v["sequence"] else None,
            "Source": "/".join(v["source"]) if v["source"] else "",
        }
        for k, v in sorted(fallback_info.items())
    ])
    mapping_df = pd.DataFrame(mapping_records)

    n_missing_accession = df_llm["Degradation_Target_Uniprot"].isna().sum()
    n_gene = mapping_df["Gene"].notna().sum() if not mapping_df.empty else 0
    n_mut = mapping_df["MutationID"].notna().sum() if not mapping_df.empty else 0
    print(f"\nRows without Degradation_Target_Uniprot: {n_missing_accession}/{len(df_llm)}")
    print(f"Gene symbols found: {n_gene}/{len(mapping_df)}")
    print(f"Mutation IDs found: {n_mut}/{len(mapping_df)}")

    df_llm = _place_after(df_llm, "Degradation_Target", [
        "Degradation_Target_Gene", "Degradation_Target_Uniprot",
        "Degradation_Target_Uniprot_MutationID", "Degradation_Target_Sequence",
    ])

    return df_llm


def clean_and_enrich_poi(df_llm: pd.DataFrame, taxid_to_species: dict) -> pd.DataFrame:
    df_llm["Degradation_Target"] = df_llm["Degradation_Target"].astype("string").str.strip()
    df_llm["Degradation_Target"] = df_llm["Degradation_Target"].replace(poi_alias_map)

    if "Target_Uniprot" in df_llm.columns or "Degradation_Target_Uniprot" in df_llm.columns:
        print("Using provided target UniProt accessions for POI enrichment")
        return clean_and_enrich_poi_from_uniprot(df_llm, taxid_to_species)

    pairs = (
        df_llm[["Degradation_Target", "_taxonomy_id"]]
        .dropna(subset=["Degradation_Target"])
        .drop_duplicates()
        .values.tolist()
    )
    pairs.sort(key=lambda x: (x[0], x[1]))
    print(f"Unique (target, organism) pairs: {len(pairs)}")

    poi_info: dict[tuple[str, int], dict] = {}

    for i, (poi, org_id) in enumerate(pairs):
        org_id = int(org_id)
        species = taxid_to_species.get(org_id, str(org_id))
        print(f"[{i+1}/{len(pairs)}] {poi}  ({species})")
        info = {"gene": None, "uniprot": None, "mutation_id": None, "sequence": None, "source": []}

        lookup = POI_LOOKUP.get(poi)
        if lookup is None:
            print(f"  WARNING: '{poi}' not in POI_LOOKUP, skipping")
            poi_info[(poi, org_id)] = info
            continue

        search_term = lookup["search"]
        mutation_str = lookup["mutation"]

        if search_term is None:
            print(f"  (no gene)")
            poi_info[(poi, org_id)] = info
            continue

        # Step 1: Gene symbol (HGNC for human, UniProt for others)
        gene, gene_src = find_gene_symbol(search_term, org_id)
        info["gene"] = gene
        if gene:
            src_tag = gene_src.split("(")[0].split(" ")[0]
            info["source"].append(src_tag)
            print(f"  Gene: {gene}  [{gene_src}]")
        else:
            print(f"  Gene: MISSING  [{gene_src}]")
            poi_info[(poi, org_id)] = info
            continue

        accession, sequence = fetch_uniprot_accession_and_seq(gene, org_id)
        info["uniprot"] = accession
        info["sequence"] = sequence
        if accession:
            if "UniProt" not in info["source"]:
                info["source"].append("UniProt")
            print(f"  UniProt: {accession}  (seq length: {len(sequence) if sequence else 0})")
        else:
            print(f"  UniProt: not found for gene {gene} (organism {org_id})")
            poi_info[(poi, org_id)] = info
            continue

        parsed = parse_mutation(mutation_str)
        if not parsed:
            poi_info[(poi, org_id)] = info
            continue

        wt, pos, mut = parsed

        if org_id == HUMAN:
            clinvar_id = clinvar_check_variant(gene, mutation_str)
            if clinvar_id:
                info["source"].append("ClinVar")
                print(f"  ClinVar: variant confirmed (ID={clinvar_id})")

        var_id = find_var_id(accession, wt, pos, mut)
        if var_id:
            info["mutation_id"] = var_id
            print(f"  Mutation {mutation_str}: {var_id}")
        else:
            print(f"  Mutation {mutation_str}: no UniProt VAR ID found")

        if sequence:
            mutated = apply_mutation(sequence, wt, pos, mut)
            if mutated:
                info["sequence"] = mutated
                print(f"  Sequence: applied {mutation_str} at position {pos}")

        poi_info[(poi, org_id)] = info
        time.sleep(0.1)

    def _get_field(row, field):
        key = (row["Degradation_Target"], int(row["_taxonomy_id"]))
        return poi_info.get(key, {}).get(field)

    df_llm["Degradation_Target_Gene"] = df_llm.apply(lambda r: _get_field(r, "gene"), axis=1)
    df_llm["Degradation_Target_Uniprot"] = df_llm.apply(lambda r: _get_field(r, "uniprot"), axis=1)
    df_llm["Degradation_Target_Uniprot_MutationID"] = df_llm.apply(lambda r: _get_field(r, "mutation_id"), axis=1)
    df_llm["Degradation_Target_Sequence"] = df_llm.apply(lambda r: _get_field(r, "sequence"), axis=1)

    df_llm.drop(columns=["POI_raw", "POI_Name"], inplace=True, errors="ignore")

    mapping_df = pd.DataFrame([
        {
            "Degradation_Target": k[0],
            "Organism": taxid_to_species.get(k[1], str(k[1])),
            "Gene": v["gene"],
            "UniProt": v["uniprot"],
            "MutationID": v["mutation_id"],
            "Seq_Length": len(v["sequence"]) if v["sequence"] else None,
            "Source": "/".join(v["source"]) if v["source"] else "",
        }
        for k, v in sorted(poi_info.items())
    ])

    n_gene = mapping_df["Gene"].notna().sum()
    n_uniprot = mapping_df["UniProt"].notna().sum()
    n_mut = mapping_df["MutationID"].notna().sum()
    n_nonhuman = (mapping_df["Organism"] != "Human").sum()
    print(f"\nUnique (target, organism) pairs: {len(mapping_df)}")
    print(f"  Non-human pairs: {n_nonhuman}")
    print(f"Gene symbols found: {n_gene}/{len(mapping_df)}")
    print(f"UniProt accessions found: {n_uniprot}/{len(mapping_df)}")
    print(f"Mutation IDs found: {n_mut}/{len(mapping_df)}")

    df_llm = _place_after(df_llm, "Degradation_Target", [
        "Degradation_Target_Gene", "Degradation_Target_Uniprot",
        "Degradation_Target_Uniprot_MutationID", "Degradation_Target_Sequence",
    ])

    return df_llm


# ===========================================================================
# E3-Ligase (Recruiter) cleaning + enrichment
# ===========================================================================

e3_alias_map = {
    "CEREBLON":           "CRBN",
    "Cereblon":           "CRBN",
    "cereblon":           "CRBN",
    "CRBN (cereblon)":    "CRBN",
    "CRBN (pomalidomide)":"CRBN",
    "Cereblon (CRBN)":    "CRBN",
    "Unknown":            None,
}

E3_LOOKUP = {
    "CRBN":   "CRBN",
    "VHL":    "VHL",
    "IAP":    "BIRC2",
    "cIAP1":  "BIRC2",
    "MDM2":   "MDM2",
    "DCAF1":  "DCAF1",
    "HRD1":   "SYVN1",
    "KLHL20": "KLHL20",
    # --- Molecular glue recruiters ---
    "CRL2KLHDC3":                  "KLHDC3",
    "DCAF16":                      "DCAF16",
    "DCAF16 (CRL4DCAF16)":         "DCAF16",
    "DDB1 (CRL4-RBX1 E3 ligase)":  "DDB1",
    "DDB1 (CRL4B)":                "DDB1",
    "DDB1-CUL4-RBX1":              "DDB1",
    "DDB1–CUL4–RBX1":              "DDB1",
    "RNF126":                      "RNF126",
    "SIAH1":                       "SIAH1",
    # --- PROTACDB baseline recruiters ---
    "FEM1B":    "FEM1B",
    "Keap1":    "KEAP1",
    "XIAP":     "XIAP",
}


def clean_and_enrich_e3(df_llm: pd.DataFrame, taxid_to_species: dict) -> pd.DataFrame:
    df_llm["Recruiter"] = df_llm["Recruiter"].astype("string").str.strip()
    df_llm["Recruiter"] = df_llm["Recruiter"].replace(e3_alias_map)

    e3_after = sorted(df_llm["Recruiter"].dropna().unique())
    print(f"Unique Recruiter values: {len(e3_after)}")

    e3_pairs = (
        df_llm[["Recruiter", "_taxonomy_id"]]
        .dropna(subset=["Recruiter"])
        .drop_duplicates()
        .values.tolist()
    )
    e3_pairs.sort(key=lambda x: (x[0], x[1]))
    print(f"Unique (recruiter, organism) pairs: {len(e3_pairs)}")

    e3_info: dict[tuple[str, int], dict] = {}

    for i, (e3, org_id) in enumerate(e3_pairs):
        org_id = int(org_id)
        species = taxid_to_species.get(org_id, str(org_id))
        print(f"[{i+1}/{len(e3_pairs)}] {e3}  ({species})")
        info = {"gene": None, "uniprot": None, "sequence": None, "source": []}

        search_term = E3_LOOKUP.get(e3)
        if search_term is None:
            print(f"  WARNING: '{e3}' not in E3_LOOKUP, skipping")
            e3_info[(e3, org_id)] = info
            continue

        gene, gene_src = find_gene_symbol(search_term, org_id)
        info["gene"] = gene
        if gene:
            src_tag = gene_src.split("(")[0].split(" ")[0]
            info["source"].append(src_tag)
            print(f"  Gene: {gene}  [{gene_src}]")
        else:
            print(f"  Gene: MISSING  [{gene_src}]")
            e3_info[(e3, org_id)] = info
            continue

        accession, sequence = fetch_uniprot_accession_and_seq(gene, org_id)
        info["uniprot"] = accession
        info["sequence"] = sequence
        if accession:
            if "UniProt" not in info["source"]:
                info["source"].append("UniProt")
            print(f"  UniProt: {accession}  (seq length: {len(sequence) if sequence else 0})")
        else:
            print(f"  UniProt: not found for gene {gene} (organism {org_id})")

        e3_info[(e3, org_id)] = info
        time.sleep(0.1)

    def _get_e3_field(row, field):
        key = (row["Recruiter"], int(row["_taxonomy_id"]))
        return e3_info.get(key, {}).get(field)

    df_llm["Recruiter_Gene"] = df_llm.apply(lambda r: _get_e3_field(r, "gene"), axis=1)
    df_llm["Recruiter_Uniprot"] = df_llm.apply(lambda r: _get_e3_field(r, "uniprot"), axis=1)
    df_llm["Recruiter_Sequence"] = df_llm.apply(lambda r: _get_e3_field(r, "sequence"), axis=1)

    df_llm.drop(columns=["_taxonomy_id"], inplace=True)

    e3_mapping_df = pd.DataFrame([
        {
            "Recruiter": k[0],
            "Organism": taxid_to_species.get(k[1], str(k[1])),
            "Gene": v["gene"],
            "UniProt": v["uniprot"],
            "Seq_Length": len(v["sequence"]) if v["sequence"] else None,
            "Source": "/".join(v["source"]) if v["source"] else "",
        }
        for k, v in sorted(e3_info.items())
    ])

    n_gene = e3_mapping_df["Gene"].notna().sum()
    n_uniprot = e3_mapping_df["UniProt"].notna().sum()
    n_nonhuman = (e3_mapping_df["Organism"] != "Human").sum()
    print(f"\nUnique (recruiter, organism) pairs: {len(e3_mapping_df)}")
    print(f"  Non-human pairs: {n_nonhuman}")
    print(f"Gene symbols found: {n_gene}/{len(e3_mapping_df)}")
    print(f"UniProt accessions found: {n_uniprot}/{len(e3_mapping_df)}")

    df_llm = _place_after(df_llm, "Recruiter", [
        "Recruiter_Gene", "Recruiter_Uniprot", "Recruiter_Sequence",
    ])

    return df_llm


# ===========================================================================
# DC50 / Dmax / Dmax_conc standardization
# ===========================================================================

_RAW_COLS = ["DC50_raw", "DC50_h_raw", "Dmax_raw", "Dmax_h_raw", "Dmax_conc_raw"]
_INSERT_AFTER = {
    "DC50_h": ["DC50_operator", "DC50_range_min", "DC50_range_max"],
    "Dmax_h": ["Dmax_operator", "Dmax_range_min", "Dmax_range_max"],
    "Dmax_conc": [
        "Dmax_conc_operator", "Dmax_conc_range_min",
        "Dmax_conc_range_max", "Dmax_conc_units",
    ],
}


def reorder_standardized_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.drop(columns=[c for c in _RAW_COLS if c in df.columns])
    moved = {c for cols in _INSERT_AFTER.values() for c in cols}
    base_cols = [c for c in df.columns if c not in moved]
    new_cols = []
    for col in base_cols:
        new_cols.append(col)
        for follower in _INSERT_AFTER.get(col, []):
            if follower in df.columns:
                new_cols.append(follower)
    for col in df.columns:
        if col not in new_cols:
            new_cols.append(col)
    return df[new_cols]


_MONTH_TO_NUMBER = {
    "jan": 1,
    "feb": 2,
    "mar": 3,
    "apr": 4,
    "may": 5,
    "jun": 6,
    "jul": 7,
    "aug": 8,
    "sep": 9,
    "oct": 10,
    "nov": 11,
    "dec": 12,
}


def _excel_date_range_bounds(s: str) -> tuple[float, float] | None:
    m = re.match(r"^([0-9]{1,2})-([A-Za-z]{3})$", s.strip())
    if not m:
        return None
    hi = float(m.group(1))
    lo = _MONTH_TO_NUMBER.get(m.group(2).lower())
    if lo is None:
        return None
    return float(lo), hi


def parse_dc50(raw: str) -> dict:
    if pd.isna(raw):
        return {"value": None, "operator": None, "range_min": None, "range_max": None}

    s = str(raw).strip()
    if not s:
        return {"value": None, "operator": None, "range_min": None, "range_max": None}
    if s.strip("()").lower() in {"n/a", "na", "nan", "none"}:
        return {"value": None, "operator": None, "range_min": None, "range_max": None}

    s = s.replace(",", "")

    s = (
        s.replace("∼", "~")
        .replace("≈", "~")
        .replace("≃", "~")
        .replace("≅", "~")
        .replace("～", "~")
        .replace("≥", ">=")
        .replace("≤", "<=")
        .replace("–", "-")
    )
    s = re.sub(r"^DC_?50\s*[:=]?\s*", "", s, flags=re.IGNORECASE).strip()
    s = re.sub(r"\s*nM\s*$", "", s, flags=re.IGNORECASE).strip()

    excel_date_range = _excel_date_range_bounds(s)
    if excel_date_range:
        lo, hi = excel_date_range
        return {"value": None, "operator": None, "range_min": lo, "range_max": hi}

    m = re.match(r'^([0-9.]+)\s*(?:-|~)\s*([0-9.]+)$', s)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        return {"value": None, "operator": None, "range_min": lo, "range_max": hi}

    m = re.match(r'^(<=|<)\s*([0-9.]+)$', s)
    if m:
        val = float(m.group(2))
        return {"value": val, "operator": m.group(1), "range_min": None, "range_max": None}

    m = re.match(r'^(>=|>)\s*([0-9.]+)$', s)
    if m:
        val = float(m.group(2))
        return {"value": val, "operator": m.group(1), "range_min": None, "range_max": None}

    m = re.match(r'^~\s*([0-9.]+)$', s)
    if m:
        val = float(m.group(1))
        return {"value": val, "operator": "~", "range_min": None, "range_max": None}

    m = re.match(r'^([0-9.]+)\s*±\s*([0-9.]+)$', s)
    if m:
        center, err = float(m.group(1)), float(m.group(2))
        return {"value": center, "operator": None, "range_min": center - err, "range_max": center + err}

    m = re.match(r'^([0-9.]+)$', s)
    if m:
        val = float(m.group(1))
        return {"value": val, "operator": None, "range_min": None, "range_max": None}

    print(f"  WARNING: could not parse DC50 value: {raw!r}")
    return {"value": None, "operator": None, "range_min": None, "range_max": None}


def parse_dc50_h(raw) -> int | float | None:
    if pd.isna(raw):
        return None

    s = str(raw).strip()
    if not s:
        return None
    if s.strip("()").lower() in {"n/a", "na", "nan", "none"}:
        return None
    if s.lower() == "overnight":
        return 16

    s = (
        s.replace("≥", ">=")
        .replace("≤", "<=")
        .replace("≈", "~")
        .replace("∼", "~")
        .replace("–", "-")
    )
    s = re.sub(r"^(?:>=|<=|>|<|~)\s*", "", s).strip()
    s = re.sub(r"\s*h(?:ours?)?\s*$", "", s, flags=re.IGNORECASE).strip()

    excel_date_range = _excel_date_range_bounds(s)
    if excel_date_range:
        lo, hi = excel_date_range
        mid = (lo + hi) / 2
        return int(mid) if mid.is_integer() else mid

    m = re.match(r"^([0-9.]+)\s*-\s*([0-9.]+)$", s)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        mid = (lo + hi) / 2
        return int(mid) if mid.is_integer() else mid

    m = re.match(r"^([0-9.]+)$", s)
    if m:
        val = float(m.group(1))
        return int(val) if val.is_integer() else val

    print(f"  WARNING: could not parse DC50_h value: {raw!r}")
    return None


def _normalize_hour(v):
    if v is None or pd.isna(v):
        return pd.NA
    val = float(v)
    return int(val) if val.is_integer() else val


def _normalize_numeric(v):
    if v is None or pd.isna(v):
        return pd.NA
    val = float(v)
    return int(val) if val.is_integer() else val


def _normalize_numeric_series(series: pd.Series) -> pd.Series:
    normalized = pd.Series(
        [_normalize_numeric(v) for v in series],
        index=series.index,
        dtype=object,
    )
    has_fraction = normalized.dropna().apply(lambda v: not float(v).is_integer()).any()
    return normalized if has_fraction else normalized.astype("Int64")


def standardize_dc50(df_llm: pd.DataFrame, input_csv: str) -> pd.DataFrame:
    dc50_raw = pd.read_csv(input_csv, usecols=["DC50"])["DC50"]
    parsed = dc50_raw.apply(parse_dc50).apply(pd.Series)

    df_llm["DC50"] = _normalize_numeric_series(parsed["value"])
    df_llm["DC50_operator"] = parsed["operator"]
    df_llm["DC50_range_min"] = _normalize_numeric_series(parsed["range_min"])
    df_llm["DC50_range_max"] = _normalize_numeric_series(parsed["range_max"])

    df_llm["DC50_units"] = (
        df_llm["DC50_units"]
        .astype("string")
        .str.strip()
        .replace({"nm": "nM", "nmol/L": "nM"})
    )

    dc50_h_raw = pd.read_csv(input_csv, usecols=["DC50_h"], dtype=str)["DC50_h"]
    dc50_h_parsed = dc50_h_raw.apply(parse_dc50_h)

    dc50_h_normalized = dc50_h_parsed.map(_normalize_hour)
    has_fraction = dc50_h_normalized.dropna().apply(lambda v: not float(v).is_integer()).any()
    if has_fraction:
        df_llm["DC50_h"] = dc50_h_normalized.astype(object)
    else:
        df_llm["DC50_h"] = dc50_h_normalized.astype("Int64")

    df_llm = reorder_standardized_columns(df_llm)
    return df_llm


def parse_dmax(raw) -> dict:
    if pd.isna(raw):
        return {"value": None, "operator": None, "range_min": None, "range_max": None}

    s = str(raw).strip()
    if not s or s == "-":
        return {"value": None, "operator": None, "range_min": None, "range_max": None}
    if s.strip("()").lower() in {"n/a", "na", "nan", "none"}:
        return {"value": None, "operator": None, "range_min": None, "range_max": None}

    s = (
        s.replace("∼", "~")
        .replace("≈", "~")
        .replace("≃", "~")
        .replace("≅", "~")
        .replace("～", "~")
        .replace("≥", ">=")
        .replace("≤", "<=")
        .replace("–", "-")
    )

    if re.match(r"^-\s*[0-9.]+$", s):
        return {"value": None, "operator": None, "range_min": None, "range_max": None}

    m = re.match(r'^(<=|<)\s*([0-9.]+)$', s)
    if m:
        val = float(m.group(2))
        return {"value": val, "operator": m.group(1), "range_min": None, "range_max": None}

    m = re.match(r'^(>=|>)\s*([0-9.]+)$', s)
    if m:
        val = float(m.group(2))
        return {"value": val, "operator": m.group(1), "range_min": None, "range_max": None}

    m = re.match(r'^~\s*([0-9.]+)$', s)
    if m:
        val = float(m.group(1))
        return {"value": val, "operator": "~", "range_min": None, "range_max": None}

    m = re.match(r'^([0-9.]+)\s*-\s*([0-9.]+)$', s)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        return {"value": None, "operator": None, "range_min": lo, "range_max": hi}

    m = re.match(r'^([0-9.]+)\s*±\s*([0-9.]+)$', s)
    if m:
        center, err = float(m.group(1)), float(m.group(2))
        return {"value": center, "operator": None, "range_min": center - err, "range_max": center + err}

    m = re.match(r'^([0-9.]+)$', s)
    if m:
        val = float(m.group(1))
        return {"value": val, "operator": None, "range_min": None, "range_max": None}

    print(f"  WARNING: could not parse Dmax value: {raw!r}")
    return {"value": None, "operator": None, "range_min": None, "range_max": None}


def parse_dmax_h(raw) -> int | float | None:
    if pd.isna(raw):
        return None

    s = str(raw).strip()
    if not s:
        return None

    s = re.sub(r"\s*h\s*$", "", s, flags=re.IGNORECASE).strip()
    s = s.replace("–", "-")

    m = re.match(r"^([0-9.]+)\s*-\s*([0-9.]+)$", s)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        val = (lo + hi) / 2
        return int(val) if val.is_integer() else val

    m = re.match(r"^([0-9.]+)$", s)
    if m:
        val = float(m.group(1))
        return int(val) if val.is_integer() else val

    print(f"  WARNING: could not parse Dmax_h value: {raw!r}")
    return None


def standardize_dmax(df_llm: pd.DataFrame, input_csv: str) -> pd.DataFrame:
    dmax_raw = pd.read_csv(input_csv, usecols=["Dmax"])["Dmax"]
    parsed = dmax_raw.apply(parse_dmax).apply(pd.Series)

    df_llm["Dmax"] = _normalize_numeric_series(parsed["value"])
    df_llm["Dmax_operator"] = parsed["operator"]
    df_llm["Dmax_range_min"] = _normalize_numeric_series(parsed["range_min"])
    df_llm["Dmax_range_max"] = _normalize_numeric_series(parsed["range_max"])

    dmax_h_raw = pd.read_csv(input_csv, usecols=["Dmax_h"], dtype=str)["Dmax_h"]
    dmax_h_parsed = dmax_h_raw.apply(parse_dmax_h)

    dmax_h_normalized = dmax_h_parsed.map(_normalize_hour)
    has_fraction = dmax_h_normalized.dropna().apply(lambda v: not float(v).is_integer()).any()
    if has_fraction:
        df_llm["Dmax_h"] = dmax_h_normalized.astype(object)
    else:
        df_llm["Dmax_h"] = dmax_h_normalized.astype("Int64")

    df_llm = reorder_standardized_columns(df_llm)
    return df_llm


_EMPTY = {"value": None, "operator": None, "range_min": None, "range_max": None}


def _conc_as_nM(amount: float, unit: str) -> int | float | None:
    u = unit.strip().lower().replace("µ", "μ")
    if u in {"nm", "nmol/l"}:
        factor = 1.0
    elif u in {"μm", "um", "μmol/l", "umol/l"}:
        factor = 1000.0
    else:
        return None
    val = amount * factor
    return int(val) if float(val).is_integer() else float(val)


def parse_dmax_conc(raw) -> dict:
    if pd.isna(raw):
        return _EMPTY.copy()

    s = str(raw).strip()
    if not s or re.search(r"mg/kg", s, flags=re.IGNORECASE):
        return _EMPTY.copy()

    s = s.replace("–", "-").replace("≥", ">=").strip()

    m = re.match(
        r"^([0-9.]+)\s*(\S+)\s+and\s+([0-9.]+)\s*(\S+)$",
        s,
        flags=re.IGNORECASE,
    )
    if m:
        lo = _conc_as_nM(float(m.group(1)), m.group(2))
        hi = _conc_as_nM(float(m.group(3)), m.group(4))
        if lo is None or hi is None:
            print(f"  WARNING: could not parse Dmax_conc multi-value: {raw!r}")
            return _EMPTY.copy()
        return {"value": None, "operator": None, "range_min": lo, "range_max": hi}

    m = re.match(r"^([0-9.]+)\s*-\s*([0-9.]+)\s*(\S+)$", s)
    if m:
        lo = _conc_as_nM(float(m.group(1)), m.group(3))
        hi = _conc_as_nM(float(m.group(2)), m.group(3))
        if lo is None or hi is None:
            print(f"  WARNING: could not parse Dmax_conc range: {raw!r}")
            return _EMPTY.copy()
        mid = (float(lo) + float(hi)) / 2
        mid = int(mid) if mid.is_integer() else mid
        return {"value": mid, "operator": None, "range_min": lo, "range_max": hi}

    for pattern, op in [
        (r"^>=\s*([0-9.]+)\s*(\S+)$", ">="),
        (r"^>\s*([0-9.]+)\s*(\S+)$", ">"),
        (r"^~\s*([0-9.]+)\s*(\S+)$", "~"),
    ]:
        m = re.match(pattern, s)
        if m:
            val = _conc_as_nM(float(m.group(1)), m.group(2))
            if val is None:
                print(f"  WARNING: could not parse Dmax_conc value: {raw!r}")
                return _EMPTY.copy()
            return {"value": val, "operator": op, "range_min": None, "range_max": None}

    m = re.match(r"^approaching\s+([0-9.]+)\s*(\S+)$", s, flags=re.IGNORECASE)
    if m:
        val = _conc_as_nM(float(m.group(1)), m.group(2))
        if val is None:
            print(f"  WARNING: could not parse Dmax_conc value: {raw!r}")
            return _EMPTY.copy()
        return {"value": val, "operator": "~", "range_min": None, "range_max": None}

    m = re.match(r"^([0-9.]+)\s*(\S+)$", s)
    if m:
        val = _conc_as_nM(float(m.group(1)), m.group(2))
        if val is None:
            print(f"  WARNING: could not parse Dmax_conc value: {raw!r}")
            return _EMPTY.copy()
        return {"value": val, "operator": None, "range_min": None, "range_max": None}

    print(f"  WARNING: could not parse Dmax_conc value: {raw!r}")
    return _EMPTY.copy()


def _normalize_conc(v):
    if v is None or pd.isna(v):
        return pd.NA
    val = float(v)
    return int(val) if val.is_integer() else val


def standardize_dmax_conc(df_llm: pd.DataFrame, input_csv: str) -> pd.DataFrame:
    dmax_conc_raw = pd.read_csv(input_csv, usecols=["Dmax_conc"], dtype=str)["Dmax_conc"]
    parsed = dmax_conc_raw.apply(parse_dmax_conc).apply(pd.Series)

    for col in ["value", "range_min", "range_max"]:
        parsed[col] = parsed[col].map(_normalize_conc)

    conc_values = pd.concat([parsed["value"], parsed["range_min"], parsed["range_max"]]).dropna()
    has_fraction = conc_values.apply(lambda v: not float(v).is_integer()).any()
    conc_dtype = object if has_fraction else "Int64"

    df_llm["Dmax_conc"] = parsed["value"].astype(conc_dtype)
    df_llm["Dmax_conc_operator"] = parsed["operator"]
    df_llm["Dmax_conc_range_min"] = parsed["range_min"].astype(conc_dtype)
    df_llm["Dmax_conc_range_max"] = parsed["range_max"].astype(conc_dtype)

    has_conc = (
        df_llm["Dmax_conc"].notna()
        | df_llm["Dmax_conc_range_min"].notna()
        | df_llm["Dmax_conc_range_max"].notna()
    )
    df_llm["Dmax_conc_units"] = pd.NA
    df_llm.loc[has_conc, "Dmax_conc_units"] = "nM"

    df_llm = reorder_standardized_columns(df_llm)
    return df_llm


# ===========================================================================
# Orchestration
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Standardize the extracted dataset."
    )
    parser.add_argument("--input_csv", type=str, default=DEFAULT_INPUT_CSV,
                        help="Path to the extracted CSV.")
    parser.add_argument("--output_csv", type=str, default=DEFAULT_OUTPUT_CSV,
                        help="Path to write the standardized CSV.")
    args = parser.parse_args()

    input_csv = args.input_csv
    output_csv = args.output_csv

    if not Path(input_csv).exists():
        print(f"Error: input CSV does not exist - {input_csv}")
        raise SystemExit(1)

    print(f"{'='*80}")
    print(f"Data standardization pipeline")
    print(f"  Input:  {input_csv}")
    print(f"  Output: {output_csv}")
    print(f"{'='*80}\n")

    df_llm = pd.read_csv(input_csv)

    print("\n--- Stage 1/4: Cell line standardization ---")
    df_llm, cl_info = standardize_cell_lines(df_llm)
    taxid_to_species = build_taxid_to_species(cl_info)

    print("\n--- Stage 2/4: Degradation target enrichment ---")
    df_llm = clean_and_enrich_poi(df_llm, taxid_to_species)

    print("\n--- Stage 3/4: Recruiter enrichment ---")
    df_llm = clean_and_enrich_e3(df_llm, taxid_to_species)

    print("\n--- Stage 4/4: DC50 / Dmax / Dmax_conc standardization ---")
    df_llm = standardize_dc50(df_llm, input_csv)
    df_llm = standardize_dmax(df_llm, input_csv)
    df_llm = standardize_dmax_conc(df_llm, input_csv)

    df_llm = _place_after(df_llm, "IUPAC_Name", ["IUPAC_Name_Source"])

    out_path = Path(output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_llm.to_csv(out_path, index=False)

    print(f"\n{'='*80}")
    print(f"Standardization complete: {len(df_llm)} rows, {len(df_llm.columns)} columns")
    print(f"Saved to: {out_path}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
