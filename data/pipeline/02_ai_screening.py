"""
02_ai_screening.py -- rule-based keyword classifier for title/abstract screening (advisory suggestions).

The classifier was written by the AI assistant (Claude, Anthropic). Its first version (v1) required both a domain
term and a method term for Include and agreed poorly (Cohen's kappa = 0.180, 45.3% agreement) with the assistant's
own reading of 75 records (25 per v1 suggestion class), made without seeing the classifier's suggestions. v2 adds a
STRONG_GENETICS_METHODS carve-out (Wright-Fisher, chemical master equation, Moran model, coalescent, branching
process, Potts model, diffusion approximation count as sufficient genetics-domain evidence) and, re-scored on the
same 75 records, reached kappa = 0.365 (70.7% agreement). That value is in-sample (the same records were used for
tuning) and both readings come from the same assistant. The record-level calibration readings were not kept; only
the summary table survives (workbook sheet Calibration_Validation), so the calibration cannot be re-run.

The suggestions are advisory. pipeline/03_merge_screening_decisions.py applies the confirmation pass, made by the AI
assistant at the first author's instruction, that produces the final decisions (manuscript Section 2.2).

Usage: python pipeline/02_ai_screening.py   (run from data/)
"""
import pandas as pd

IN_PATH = "pipeline_output/01_merged_deduplicated.csv"
OUT_PATH = "pipeline_output/02_merged_with_ai_screening.csv"

POSITIVE_TERMS = [
    "population genetics", "quantitative genetics", "genomic selection",
    "marker-assisted", "plant breeding", "animal breeding", "crop breeding",
    "allele frequenc", "heritability", "genetic drift", "genetic variance",
    "molecular evolution", "evolutionary genetics", "quantitative trait loci",
    "genome-wide association", "genomic prediction", "gene regulatory network",
    "gene expression", "population genomics", "mutation rate",
    "fixation probability", "selection coefficient", "fitness landscape",
    "breeding value", "pedigree", "inbreeding", "heterozygos", "recombination",
    "linkage disequilibrium", "quantitative trait", "genome-wide", "allele",
    "gene network", "transcription factor", "gene circuit", "genotype",
    "phenotype", "locus", "loci", "hered", "genetic variation",
    "evolutionary biology", "genetic architecture", "trait evolution",
    "fixation", "deleterious", "speciation", "transposable element",
    "kimura", "coevolutionary", "biochemical switch", "mrna",
    "protein network", "genetic circuit", "competence", "altruism",
    "demographic", "site frequency spectrum", "neutrality test",
    "population dynamics", "evolutionary game", "epistasis", "haplotype",
]

METHOD_TERMS = [
    "stochastic differential equation", "chemical master equation",
    "wright-fisher", "wright–fisher", "diffusion approximation",
    "fokker-planck", "fokker–planck", "hamilton-jacobi-bellman",
    "mean field game", "ising model", "potts model", "master equation",
    "branching process", "score-based generative", "neural sde",
    "stochastic optimal control", "stochastic process", "markov chain",
    "moran model", "birth-death process", "diffusion process", "coalescent",
]

STRONG_GENETICS_METHODS = {
    "wright-fisher", "wright–fisher", "chemical master equation",
    "moran model", "coalescent", "diffusion approximation",
    "branching process", "potts model",
}

NEGATIVE_TERMS = [
    "genetic algorithm", "evolutionary algorithm", "particle swarm",
    "epidemic", "covid-19", "sars-cov", "pandemic", "financial market",
    "stock price", "option pricing", "image segmentation", "computer vision",
    "autonomous vehicle", "wireless network", "traffic flow", "power grid",
    "supply chain", "quantum",
]


def count_hits(text, terms):
    return [t for t in terms if t in text]


def classify(row):
    if pd.isna(row["abstract"]) or not str(row["abstract"]).strip():
        return pd.Series({
            "ai_suggested_decision": "Maybe",
            "ai_suggested_reason": "insufficient_information",
            "ai_confidence": "n/a",
            "ai_flags": "missing_abstract",
        })

    text = (str(row["title"]) + " " + str(row["abstract"])).lower()
    pos_hits = count_hits(text, POSITIVE_TERMS)
    method_hits = count_hits(text, METHOD_TERMS)
    neg_hits = count_hits(text, NEGATIVE_TERMS)
    strong_method_hits = [m for m in method_hits if m in STRONG_GENETICS_METHODS]

    flags = []
    if pos_hits:
        flags.append("pos:" + ",".join(sorted(set(pos_hits))[:5]))
    if method_hits:
        flags.append("method:" + ",".join(sorted(set(method_hits))[:5]))
    if neg_hits:
        flags.append("neg:" + ",".join(sorted(set(neg_hits))))
    flag_str = " | ".join(flags)

    if neg_hits and not pos_hits and not strong_method_hits:
        return pd.Series({"ai_suggested_decision": "Exclude", "ai_suggested_reason": "keyword_coincidence",
                           "ai_confidence": "high", "ai_flags": flag_str})
    if neg_hits and (pos_hits or strong_method_hits):
        return pd.Series({"ai_suggested_decision": "Maybe", "ai_suggested_reason": "none",
                           "ai_confidence": "low", "ai_flags": flag_str})
    if strong_method_hits and pos_hits:
        return pd.Series({"ai_suggested_decision": "Include", "ai_suggested_reason": "none",
                           "ai_confidence": "high", "ai_flags": flag_str})
    if strong_method_hits and not pos_hits:
        return pd.Series({"ai_suggested_decision": "Include", "ai_suggested_reason": "none",
                           "ai_confidence": "medium", "ai_flags": flag_str})
    if len(pos_hits) >= 2 and method_hits:
        return pd.Series({"ai_suggested_decision": "Include", "ai_suggested_reason": "none",
                           "ai_confidence": "high", "ai_flags": flag_str})
    if len(pos_hits) == 1 and method_hits:
        return pd.Series({"ai_suggested_decision": "Include", "ai_suggested_reason": "none",
                           "ai_confidence": "medium", "ai_flags": flag_str})
    if not pos_hits and not method_hits:
        return pd.Series({"ai_suggested_decision": "Exclude", "ai_suggested_reason": "not_genetics_domain",
                           "ai_confidence": "medium", "ai_flags": flag_str})
    if len(pos_hits) == 1 and not method_hits:
        return pd.Series({"ai_suggested_decision": "Maybe", "ai_suggested_reason": "none",
                           "ai_confidence": "low", "ai_flags": flag_str})
    return pd.Series({"ai_suggested_decision": "Maybe", "ai_suggested_reason": "none",
                       "ai_confidence": "low", "ai_flags": flag_str})


df = pd.read_csv(IN_PATH, encoding="utf-8-sig")
result = df.apply(classify, axis=1)
out = pd.concat([df, result], axis=1)
out.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")

print("Classifier suggestions (version 2, deduplicated records):")
print(out["ai_suggested_decision"].value_counts())
print()
print(f"Saved {OUT_PATH}")
