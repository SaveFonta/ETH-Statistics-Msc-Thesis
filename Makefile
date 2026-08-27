.PHONY: all thesis pull-output push-scripts archive-manuscript appendix clean clean-all help

# Default target: build the thesis
all: thesis

# ---------------------------------------------------------------------------
# The ETH SfS thesis (MasterThesisSfSVersionSep25/). This is the live document.
# ---------------------------------------------------------------------------
# Two stages on purpose: knit2pdf runs LaTeX too few times to populate the
# table of contents, list of figures and list of tables, so knitr and latexmk
# are invoked separately and latexmk iterates until the references settle.
thesis:
	cd MasterThesisSfSVersionSep25 && \
	Rscript -e "knitr::knit('MasterThesisSfS.Rnw')" && \
	latexmk -pdf -interaction=nonstopmode MasterThesisSfS.tex

# ---------------------------------------------------------------------------
# Remote results. All simulation and MCMC work runs on doksum, but the knitr
# chunks read ../Output/*.RDS, so those results must be present locally for
# "make thesis" to compile here. Run pull-output after any script is rerun
# remotely, otherwise the local build fails on a missing RDS or silently uses
# a stale one.
# ---------------------------------------------------------------------------
REMOTE := doksum:~/Thesis-EUII

pull-output:
	rsync -av $(REMOTE)/Output/ Output/

push-scripts:
	rsync -av scripts/ $(REMOTE)/scripts/

# ---------------------------------------------------------------------------
# Superseded article version, kept for reference only (archive/).
# Its content has been migrated into the thesis; do not develop it further.
# ---------------------------------------------------------------------------
archive-manuscript:
	cd archive/manuscript-article && \
	Rscript -e "knitr::knit2pdf('manuscript.Rnw')"

# Compile only the old standalone appendix (fast, no knitr)
appendix:
	cd archive/manuscript-article && Rscript -e "tinytex::latexmk('appendix_only.tex')"

# ---------------------------------------------------------------------------
clean:
	cd MasterThesisSfSVersionSep25 && rm -f *.aux *.log *.out *.toc *.lof *.lot \
	    *.bbl *.blg *.fls *.fdb_latexmk *.synctex.gz MasterThesisSfS.tex \
	    Chapters/*.aux Chapters/*.log Chapters/Methodology.tex Chapters/Applications.tex
	-@cd archive/manuscript-article 2>/dev/null && rm -f *.aux *.log *.out *.toc *.synctex.gz manuscript.tex

clean-all: clean
	cd MasterThesisSfSVersionSep25 && rm -rf figure cache MasterThesisSfS.pdf
	-@cd archive/manuscript-article 2>/dev/null && rm -rf figure cache manuscript.pdf

help:
	@echo "Thesis build targets:"
	@echo "  make thesis             - Build the ETH SfS thesis (knitr + latexmk)"
	@echo "  make pull-output        - Fetch Output/*.RDS from doksum (needed to build)"
	@echo "  make push-scripts       - Send scripts/ to doksum before running them"
	@echo "  make archive-manuscript - Rebuild the superseded article version"
	@echo "  make appendix           - Compile the old standalone appendix"
	@echo "  make clean              - Remove build artifacts (keeps PDFs)"
	@echo "  make clean-all          - Remove build artifacts including PDFs"
	@echo "  make help               - Show this message"
