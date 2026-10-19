"""
02_screen_records.py -- title-and-abstract screening by explicit, versioned rules.

Every deduplicated record receives a decision (Include, Maybe or Exclude), the exclusion reason, the rule that
decided it and the terms that matched. The rules implement the eligibility criteria in
../protocol/eligibility_criteria.md; the term lists and thresholds are in pipeline/screening_rules.yaml, whose
SHA-256 is written next to each decision. The same records and the same rules file always give the same
decisions: nothing is learned from the data and no service is called.

Inputs
  pipeline_output/01_merged_deduplicated.csv   title, abstract and document type of the deduplicated records
  raw/scopus/*.csv, raw/WoS/*.xls              author keywords and Scopus Index Keywords, joined as in step 13
                                               (Web of Science Keywords Plus are not used: Clarivate derives them
                                               from the titles of cited references)
  pipeline/screening_rules.yaml                term lists and thresholds
Output
  pipeline_output/03_final_decisions.csv       the records with final_decision, final_exclusion_reason,
                                               final_notes (rule and reason), screen_rule, the scores, the
                                               matched terms, rules_version and rules_sha256

Scores. A term counts once, with the weight of the field it occurs in (title and keywords 2, abstract 1).
  domain       strong genetics, evolution, breeding and gene-regulation terms (criterion I4)
  domain_bio   the same without the model names that other fields borrow from population genetics
               (Wright-Fisher, coalescent, genetic drift ...), which do not show a biological subject
  domain_weak  number of weak biological context terms (gene, cell, evolution ...)
  off_domain   terms of other subjects (criterion E1), including computing senses masked out of the text
               (a genetic algorithm, a language model), plus context terms (epidemiology, climate ...) at
               off_domain_context_weight; context terms never exclude a record through the title rule
  stochastic   a named stochastic model (I3), or general stochastic terms with a weighted score of 2 or more

Decision: the first rule that applies.
  E6   document type excluded                                          -> Exclude wrong_document_type
  E5   no abstract                                                     -> Maybe   insufficient_information
  E3   no stochastic term                                              -> Exclude not_stochastic_framework
       (Maybe when domain_bio >= population_genetics_without_model: a population-genetic subject whose
       model is not named in the abstract)
  E1a  no domain term and no weak biological term                      -> Exclude not_genetics_domain
  E1c  an off-domain subject term in the title, none of the strong
       domain terms there                                             -> Exclude not_genetics_domain
       (keyword_coincidence when the title term is a computing sense: a genetic algorithm, a language model)
  E1b  off_domain >= off_domain_dominates and off_domain >= domain_bio  -> Exclude not_genetics_domain
       (Maybe when domain_bio >= off_domain_keeps_maybe; keyword_coincidence when computing senses make up
       half the off-domain score or more)
  I    stochastic and one of: domain >= domain_strong_for_include;
       domain >= 1 with domain_weak >= domain_weak_support and a named model;
       domain_weak >= domain_weak_alone with a named model and no off-domain term
                                                                       -> Include
       (Maybe instead when the abstract is shorter than min_abstract_words, unless domain >= 2 with a named
       model: title and a short abstract are then sufficient, criterion I5)
  E2   anything else (passing mention, or the scope or model unclear)  -> Maybe

Usage: python -I pipeline/02_screen_records.py [output.csv]   (run from data/)
"""
import glob
import hashlib
import re
import sys

import numpy as np
import pandas as pd
import yaml

RULES_PATH = "pipeline/screening_rules.yaml"
IN_PATH = "pipeline_output/01_merged_deduplicated.csv"
OUT_PATH = "pipeline_output/03_final_decisions.csv"


# --------------------------------------------------------------------------------------------- text and matching
def text(x):
    return "" if x is None or (isinstance(x, float) and np.isnan(x)) else str(x)


DASHES = dict.fromkeys(map(ord, "‐‑‒–—−"), "-")


def normalise(s):
    """Lower case, every dash as '-', single spaces."""
    return re.sub(r"\s+", " ", text(s).translate(DASHES).lower()).strip()


def term_regex(term, case_sensitive=False):
    """Whole-word pattern. It may not start inside a word or a hyphenated compound ('non-genetic' does not match
    'genetic') but may end at a hyphen ('Wright-Fisher-like' matches 'Wright-Fisher'). A trailing * matches any
    word ending; '-' and ' ' inside a term match either."""
    t = term if case_sensitive else term.lower()
    stem, wild = (t[:-1], True) if t.endswith("*") else (t, False)
    body = r"[- ]".join(re.escape(part) for part in re.split(r"[- ]", stem))
    return re.compile(r"(?<![\w-])" + body + (r"\w*" if wild else r"(?!\w)"))


class Lexicon:
    def __init__(self, terms, case_sensitive=()):
        self.cs = set(case_sensitive)
        self.terms = [(t, term_regex(t)) for t in terms] + [(t, term_regex(t, True)) for t in case_sensitive]

    def hits(self, low, raw):
        return [t for t, rx in self.terms if rx.search(raw if t in self.cs else low)]


def load_rules(path=RULES_PATH):
    blob = open(path, "rb").read()
    rules = yaml.safe_load(blob)
    rules["sha256"] = hashlib.sha256(blob).hexdigest()
    rules["mask_rx"] = [term_regex(m) for m in rules["mask"]]
    rules["mask_off_rx"] = [(m, term_regex(m)) for m in rules["mask_off_domain"]]
    rules["model_terms"] = set(rules["domain_model_terms"])
    rules["lex"] = {
        "stoch_strong": Lexicon(rules["stochastic_strong"], rules.get("case_sensitive_stochastic", [])),
        "stoch_weak": Lexicon(rules["stochastic_weak"]),
        "dom_strong": Lexicon(rules["domain_strong"], rules.get("case_sensitive_domain", [])),
        "dom_weak": Lexicon(rules["domain_weak"]),
        "off": Lexicon(rules["off_domain_subject"]),
        "off_ctx": Lexicon(rules["off_domain_context"]),
    }
    return rules


def field_hits(record, rules):
    """{lexicon: {term: weight}}: each term with the largest field weight at which it occurs."""
    out = {k: {} for k in list(rules["lex"]) + ["off_masked"]}
    for field, weight in rules["fields"].items():
        raw = re.sub(r"\s+", " ", text(record.get(field, "")).translate(DASHES))
        low = raw.lower()
        for term, rx in rules["mask_off_rx"]:                  # computing senses: off-domain evidence, then removed
            if rx.search(low):
                out["off_masked"][term] = max(out["off_masked"].get(term, 0), weight)
                low, raw = rx.sub(" ", low), re.sub(rx.pattern, " ", raw, flags=re.I)
        for rx in rules["mask_rx"]:
            low, raw = rx.sub(" ", low), re.sub(rx.pattern, " ", raw, flags=re.I)
        for name, lex in rules["lex"].items():
            for t in lex.hits(low, raw):
                out[name][t] = max(out[name].get(t, 0), weight)
    return out


# --------------------------------------------------------------------------------------------- decision
def decide(record, rules):
    lv = rules["levels"]
    title_weight = rules["fields"]["title"]
    h = field_hits(record, rules)
    strong, weak = h["stoch_strong"], h["stoch_weak"]
    cw = rules["off_domain_context_weight"]
    off_terms = {**h["off"], **h["off_masked"], **{t: w * cw for t, w in h["off_ctx"].items()}}
    dom = sum(h["dom_strong"].values())
    dom_bio = sum(w for t, w in h["dom_strong"].items() if t not in rules["model_terms"])
    dom_weak = len(h["dom_weak"])
    off = sum(off_terms.values())
    stochastic = bool(strong) or sum(weak.values()) >= 2
    title = normalise(record.get("title", ""))
    title_off = [t for t, rx in rules["lex"]["off"].terms if rx.search(title)] + \
                [t for t, rx in rules["mask_off_rx"] if rx.search(title)]
    title_dom = [t for t, rx in rules["lex"]["dom_strong"].terms
                 if rx.search(record.get("title", "") if t in rules["lex"]["dom_strong"].cs else title)]
    abstract_words = len(text(record.get("abstract", "")).split())
    doc_type = normalise(record.get("doc_type", ""))

    def result(decision, reason, rule, why):
        terms = lambda d: "; ".join(sorted(d))
        return {
            "final_decision": decision, "final_exclusion_reason": reason, "screen_rule": rule,
            "final_notes": f"{rule}: {why}",
            "domain_score": dom, "domain_bio_score": dom_bio, "domain_weak_terms": dom_weak, "off_domain_score": off,
            "stochastic_terms": terms(strong) or terms(weak), "domain_terms": terms(h["dom_strong"]),
            "context_terms": terms(h["dom_weak"]), "off_domain_terms": terms(off_terms),
        }

    if any(re.search(r"(?<![a-z])" + re.escape(t) + r"(?![a-z])", doc_type) for t in rules["excluded_document_types"]):
        return result("Exclude", "wrong_document_type", "E6", f"document type '{record.get('doc_type')}'")
    if abstract_words == 0:
        return result("Maybe", "insufficient_information", "E5", "no abstract")
    if not strong and not weak:
        if dom_bio >= lv["population_genetics_without_model"]:
            return result("Maybe", "", "E3?", f"no stochastic term, but a population-genetic subject (domain {dom_bio})")
        return result("Exclude", "not_stochastic_framework", "E3", "no stochastic term")
    if dom == 0 and dom_weak == 0:
        return result("Exclude", "not_genetics_domain", "E1a", "no genetics, evolution, breeding or gene-regulation term")
    masked = set(h["off_masked"])
    if title_off and not title_dom:
        reason = "keyword_coincidence" if set(title_off) <= masked else "not_genetics_domain"
        return result("Exclude", reason, "E1c", f"title names another subject ({', '.join(title_off)})")
    if off >= lv["off_domain_dominates"] and off >= dom_bio:
        if dom_bio >= lv["off_domain_keeps_maybe"]:
            return result("Maybe", "", "E1b?", f"off-domain {off} against biological domain {dom_bio}")
        reason = "keyword_coincidence" if sum(h["off_masked"].values()) * 2 >= off else "not_genetics_domain"
        return result("Exclude", reason, "E1b", f"off-domain {off} against biological domain {dom_bio}")

    scope = (dom >= lv["domain_strong_for_include"]
             or (dom >= 1 and dom_weak >= lv["domain_weak_support"] and strong)
             or (dom_weak >= lv["domain_weak_alone"] and strong and off == 0))
    if stochastic and scope:
        if abstract_words < rules["min_abstract_words"] and not (dom >= 2 and strong):
            return result("Maybe", "insufficient_information", "E5",
                          f"abstract of {abstract_words} words with limited evidence")
        return result("Include", "", "I", f"domain {dom} (biological {dom_bio}, context {dom_weak}); "
                                         f"stochastic: {'named model' if strong else 'general terms'}")
    if scope:
        return result("Maybe", "", "E2s", "the stochastic model is not named (a single general term)")
    return result("Maybe", "", "E2d", f"domain evidence below the threshold (domain {dom}, context {dom_weak})")


# --------------------------------------------------------------------------------------------- records
def norm_title(t):
    return re.sub(r"[^a-z0-9]", "", text(t).lower())


def norm_doi(d):
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", text(d).strip().lower())


def year4(y):
    m = re.search(r"(\d{4})", text(y))
    return m.group(1) if m else ""


def keyword_table():
    """Author keywords and Scopus Index Keywords of each deduplicated work (same union-find as steps 01 and 13)."""
    scopus = pd.read_csv(glob.glob("raw/scopus/*.csv")[0], encoding="utf-8-sig")
    wos = pd.read_excel(glob.glob("raw/WoS/*.xls")[0])
    raw = pd.concat([
        pd.DataFrame({"title": scopus["Title"], "year": scopus["Year"], "doi": scopus["DOI"],
                      "kw": scopus["Author Keywords"].fillna("") + "; " + scopus["Index Keywords"].fillna("")}),
        pd.DataFrame({"title": wos["Article Title"], "year": wos["Publication Year"], "doi": wos["DOI"],
                      "kw": wos["Author Keywords"].fillna("")}),
    ], ignore_index=True)
    raw["tk"] = raw["title"].map(norm_title) + "|" + raw["year"].map(year4)
    raw["dk"] = raw["doi"].map(norm_doi)
    parent = list(range(len(raw)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for key, usable in (("dk", raw["dk"] != ""), ("tk", raw["title"].map(norm_title) != "")):
        groups = {}
        for i in raw.index[usable]:
            groups.setdefault(raw.at[i, key], []).append(i)
        for members in groups.values():
            for j in members[1:]:
                ri, rj = find(members[0]), find(j)
                if ri != rj:
                    parent[rj] = ri
    raw["group"] = [find(i) for i in raw.index]
    kw = raw.groupby("group")["kw"].apply(lambda s: "; ".join(x.strip("; ") for x in s if x.strip("; ")))
    return dict(zip(raw["tk"], raw["group"])), kw


def screen(records, rules):
    decided = pd.DataFrame([decide(r, rules) for r in records.to_dict("records")], index=records.index)
    out = pd.concat([records, decided], axis=1)
    out["reviewer_source"] = "screening_rules"
    out["rules_version"] = rules["version"]
    out["rules_sha256"] = rules["sha256"]
    return out


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else OUT_PATH
    rules = load_rules()
    records = pd.read_csv(IN_PATH, encoding="utf-8-sig")
    tk_group, kw = keyword_table()
    groups = (records["title"].map(norm_title) + "|" + records["year"].map(year4)).map(tk_group)
    if groups.isna().any() or not groups.is_unique:
        sys.exit("screened records do not match the raw exports one to one")
    records["keywords"] = groups.map(kw).fillna("")
    out = screen(records, rules)
    out.to_csv(out_path, index=False, encoding="utf-8-sig")

    print(f"rules {rules['version']}  sha256 {rules['sha256'][:16]}")
    print(out["final_decision"].value_counts().to_string())
    print(out.loc[out["final_decision"] == "Exclude", "final_exclusion_reason"].value_counts().to_string())
    print(out["screen_rule"].value_counts().to_string())
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
