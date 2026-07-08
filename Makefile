.PHONY: all clean help reports final manuscript appendix

# Default target
all: reports final

# Build all reports
reports: reports/First_report.pdf reports/How.RBesT.works.inside.pdf reports/Bayesian_group.pdf

# Build final manuscript
final: Final/manuscript.pdf

# Individual report targets
reports/First_report.pdf: reports/First_report.Rnw
	cd reports && Rscript -e "knitr::knit2pdf('First_report.Rnw')" && cd ..

reports/How.RBesT.works.inside.pdf: reports/How.RBesT.works.inside.Rnw
	cd reports && Rscript -e "knitr::knit2pdf('How.RBesT.works.inside.Rnw')" && cd ..

reports/Bayesian_group.pdf: reports/Bayesian_group.Rnw
	cd reports && Rscript -e "knitr::knit2pdf('Bayesian_group.Rnw')" && cd ..

# Manuscript target (incremental: rebuilds when the Rnw or any \input'd file changes)
Final/manuscript.pdf: Final/manuscript.Rnw Final/newCommands.tex Final/appendix.tex Final/Thesis.lit.rev.bib
	cd Final && Rscript -e "knitr::knit2pdf('manuscript.Rnw')" && cd ..

# Force a full recompile, even when make thinks nothing is out of date
manuscript:
	$(MAKE) -B Final/manuscript.pdf

# Compile only the appendix (fast, standalone; no knitr, no main manuscript)
appendix:
	cd Final && Rscript -e "tinytex::latexmk('appendix_only.tex')"

# Clean build artifacts (but keep PDFs for reference)
clean:
	cd reports && rm -f *.tex *.log *.aux *.out *.toc *.synctex.gz
	cd Final && rm -f *.tex *.log *.aux *.out *.toc *.synctex.gz

# Clean everything including PDFs
clean-all: clean
	cd reports && rm -f *.pdf
	cd Final && rm -f *.pdf

# Help target
help:
	@echo "Thesis Build Targets:"
	@echo "  make all              - Build all reports and manuscript (default)"
	@echo "  make reports          - Build all reports in reports/"
	@echo "  make final            - Build manuscript in Final/"
	@echo "  make manuscript       - Force a full recompile of the manuscript"
	@echo "  make appendix         - Compile only the appendix (standalone, fast)"
	@echo "  make clean            - Remove build artifacts (keeps PDFs)"
	@echo "  make clean-all        - Remove all build artifacts including PDFs"
	@echo "  make help             - Show this help message"
