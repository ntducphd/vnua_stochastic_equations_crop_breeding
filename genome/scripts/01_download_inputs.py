"""Fetch the raw inputs of 07_population_differentiation.py and 13_a1_sh2_interval.py into data/raw/
and check each against the size and MD5 in data/MANIFEST_inputs_md5.tsv.

This is the input part of the project's analysis_genome/scripts/01_download_data.py, restricted to
the files these two scripts read:
  07  B73v5_to_B73v4.tsv and Zm00001d.2.full_gene_data.txt.gz (MaizeGDB), the Hu et al. 2021 GBS VCF
      and its accession groups (Supplementary Data 4), and the Qu et al. 2026 array tables;
  13  the Pan-Zea v4 pan-gene table, the GFF3 annotation and the genome FASTA (with MaizeGDB's .md5)
      of each of the 17 assemblies listed in data/curated/genomes.csv and hybrid_parents.csv.
MaizeGDB files are downloaded from the URL in the manifest (curl when available; MaizeGDB's CDN
resets urllib connections). The three population-genetic files are not on MaizeGDB: give the folder
that holds them under the paths shown in the manifest ('project:<path>') with --local, or place
them in data/raw/popgen/ by hand (sources in README.md). Files already present are not fetched
again. The genome FASTAs are about 10 GB in all.

Usage:
  python -I scripts/01_download_inputs.py --list          # show what would be fetched; no network access
  python -I scripts/01_download_inputs.py [--local DIR]   # fetch, copy and verify
"""
import argparse
import hashlib
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = ROOT / "data" / "MANIFEST_inputs_md5.tsv"


def read_manifest():
    rows = []
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or line.startswith("path\t") or not line.strip():
            continue
        path, size, md5, _check, url = line.split("\t")
        rows.append(dict(path=path, bytes=int(size), md5=md5, url=url))
    return rows


def md5sum(path, block=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"  downloading {url}", flush=True)
    if shutil.which("curl"):
        subprocess.run(["curl", "-fsSL", "--retry", "5", "--retry-all-errors", "--retry-delay", "10", "-o", str(tmp), url],
                       check=True)
    else:
        req = urllib.request.Request(url, headers={"User-Agent": "curl/8.9.1"})
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
            while chunk := r.read(1 << 22):
                fh.write(chunk)
    tmp.replace(dest)


def published_md5(md5_file, fasta_name):
    """MaizeGDB .md5 files list '<md5>  <file name>' lines; return the FASTA's MD5 if listed."""
    if not md5_file.exists():
        return None
    for line in md5_file.read_text(errors="ignore").splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[-1].lstrip("*").endswith(fasta_name):
            return parts[0].lower()
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--list", action="store_true", help="print the inputs and exit without network access")
    ap.add_argument("--local", type=Path, help="folder holding the 'project:' files under their manifest paths")
    args = ap.parse_args()
    rows = read_manifest()
    if args.list:
        total = sum(r["bytes"] for r in rows)
        for r in rows:
            state = "present" if (RAW / r["path"]).exists() else "missing"
            print(f"{state:8} {r['bytes']:>13,}  {r['path']}  <- {r['url']}")
        print(f"{len(rows)} files, {total / 1e9:.2f} GB")
        return 0
    bad = 0
    for r in rows:
        dest = RAW / r["path"]
        if not dest.exists():
            if r["url"].startswith("https://"):
                fetch(r["url"], dest)
            elif args.local and (args.local / r["url"].split(":", 1)[1]).exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(args.local / r["url"].split(":", 1)[1], dest)
            else:
                print(f"   missing  {r['path']}: place it there by hand (README.md gives its source)")
                bad += 1
                continue
        digest = md5sum(dest)
        status = "OK" if digest == r["md5"] else "MD5 differs from the manifest"
        if r["path"].startswith("genomes/") and r["path"].endswith(".fa.gz"):
            ref = published_md5(dest.with_name(dest.name.replace(".fa.gz", ".md5")), dest.name)
            status += "; MaizeGDB .md5 " + ("not listed" if ref is None else ("OK" if ref == digest else "MISMATCH"))
        bad += digest != r["md5"]
        print(f"{status:>10}  {r['path']}")
    print("all inputs present and identical to the manifest" if not bad else f"{bad} file(s) missing or different")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
