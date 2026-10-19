"""
test_screen_records.py -- checks that each screening rule decides the case it is written for.

Usage: python -I pipeline/test_screen_records.py   (run from data/; prints PASS/FAIL per case, exits 1 on a failure)
"""
import importlib.util
import sys

spec = importlib.util.spec_from_file_location("screen", "pipeline/02_screen_records.py")
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)
RULES = screen.load_rules()

FILLER = " ".join(["word"] * 60)     # makes an abstract long enough for the length rule


def rec(title, abstract, keywords="", doc_type="Article"):
    return {"title": title, "abstract": abstract + " " + FILLER, "keywords": keywords, "doc_type": doc_type}


CASES = [
    ("wright-fisher population genetics", rec("Fixation of a beneficial allele in the Wright-Fisher model",
     "We study allele frequency dynamics under genetic drift and selection with a diffusion approximation."), "Include", "I"),
    ("named model in a hyphenated compound", rec("Wright-Fisher-like models with constant population size",
     "A Galton-Watson construction gives the fixation probability of a mutant allele."), "Include", "I"),
    ("stochastic gene expression", rec("Noise in a gene regulatory network",
     "A chemical master equation describes stochastic gene expression of a toggle switch."), "Include", "I"),
    ("erratum", rec("Erratum to: allele fixation", "Correction of a figure.", doc_type="Erratum"), "Exclude", "E6"),
    ("no abstract", {"title": "Allele fixation", "abstract": "", "keywords": "", "doc_type": "Article"}, "Maybe", "E5"),
    ("deterministic plant physiology", rec("Sugar signalling in cucumber genotypes",
     "Expression of genes under abscisic acid treatment was measured by qPCR."), "Exclude", "E3"),
    ("pure mathematics", rec("Duality of spatial birth-death processes in Hilbert space",
     "We prove a duality relation for Brownian motion on a manifold."), "Exclude", "E1a"),
    ("language change with borrowed genetics model", rec("A stochastic model of language change",
     "We adapt the Wright-Fisher model of population genetics to linguistic variants and word frequency."), "Exclude", "E1c"),
    ("genetic algorithm in the title", rec("Schema dynamics in parallel genetic algorithms",
     "Genetic algorithms with island models are analysed with a Markov chain and a diffusion approximation of genetic drift."),
     "Exclude", "E1c"),
    ("epidemiology with phylodynamics stays in", rec("Phylodynamic inference of an epidemic",
     "A birth-death master equation gives the likelihood of a pathogen phylogeny; mutation rates and the effective population size are estimated from genomes."),
     "Include", "I"),
    ("coalescence of particles is not the coalescent", rec("Stochastic models for coalescence of particles",
     "Coalescence and coagulation of polymer particles are described by a master equation."), "Exclude", "E1a"),
    ("passing mention", rec("Random walks on graphs",
     "A random walk on a lattice; such walks also appear in gene studies."), "Maybe", "E2d"),
]


def main():
    failures = 0
    for name, record, decision, rule in CASES:
        got = screen.decide(record, RULES)
        ok = got["final_decision"] == decision and got["screen_rule"] == rule
        failures += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {got['final_decision']} ({got['screen_rule']})"
              + ("" if ok else f"; expected {decision} ({rule}); {got['final_notes']}"))
    print(f"{len(CASES) - failures}/{len(CASES)} passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
