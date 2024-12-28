# Makefile to generate standard and inverted color diagrams

# Define variables
TEX=composition_diagram.tex
STANDARD_PDF=composition_diagram_light.pdf
INVERTED_PDF=composition_diagram_dark.pdf
STANDARD_SVG=composition_diagram_light.svg
INVERTED_SVG=composition_diagram_dark.svg

.PHONY: all clean allclean

all: $(STANDARD_SVG) $(INVERTED_SVG) clean

# Compile standard colors
$(STANDARD_PDF): $(TEX)
	@echo "Compiling standard color diagram..."
	pdflatex -interaction=nonstopmode -jobname=composition_diagram_light $(TEX)

# Compile inverted colors
$(INVERTED_PDF): $(TEX)
	@echo "Compiling inverted color diagram..."
	pdflatex -interaction=nonstopmode -jobname=composition_diagram_dark "\def\invertcolors{} \input{$(TEX)}"

# Convert PDFs to SVGs using pdf2svg
$(STANDARD_SVG): $(STANDARD_PDF)
	@echo "Converting standard PDF to SVG..."
	pdf2svg $(STANDARD_PDF) $(STANDARD_SVG)

$(INVERTED_SVG): $(INVERTED_PDF)
	@echo "Converting inverted PDF to SVG..."
	pdf2svg $(INVERTED_PDF) $(INVERTED_SVG)

clean:
	@echo "Cleaning up extra generated files..."
	rm -f composition_diagram*.aux composition_diagram*.log composition_diagram*.out
	rm -f composition_diagram*.pdf composition_diagram*.fdb_latexmk composition_diagram*.fls
	rm -f composition_diagram*.synctex.gz

allclean: clean
	@echo "Cleaning up generated SVGs..."
	rm -f composition_diagram*.svg
