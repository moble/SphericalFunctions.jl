# Makefile to generate standard and inverted color diagrams

# Define variables
TEX=composition_diagram.tex
STANDARD_DVI=composition_diagram_light.dvi
INVERTED_DVI=composition_diagram_dark.dvi
STANDARD_PDF=composition_diagram_light.pdf
INVERTED_PDF=composition_diagram_dark.pdf
STANDARD_SVG=composition_diagram_light.svg
INVERTED_SVG=composition_diagram_dark.svg

.PHONY: all clean allclean

all: $(STANDARD_SVG) $(INVERTED_SVG) clean

$(STANDARD_DVI): $(TEX)
	@echo "Compiling standard-color diagram to DVI..."
	lualatex -output-format=dvi -jobname=composition_diagram_light $(TEX)

$(INVERTED_DVI): $(TEX)
	@echo "Compiling inverted-color diagram to DVI..."
	lualatex -output-format=dvi -jobname=composition_diagram_dark "\def\invertcolors{} \input{$(TEX)}"

$(STANDARD_SVG): $(STANDARD_DVI)
	@echo "Converting standard DVI to SVG..."
	dvisvgm --libgs=/opt/homebrew/lib/libgs.dylib -o $(STANDARD_SVG) $(STANDARD_DVI)

$(INVERTED_SVG): $(INVERTED_DVI)
	@echo "Converting inverted DVI to SVG..."
	dvisvgm --libgs=/opt/homebrew/lib/libgs.dylib -o $(INVERTED_SVG) $(INVERTED_DVI)

# # Compile standard colors
# $(STANDARD_PDF): $(TEX)
# 	@echo "Compiling standard color diagram..."
# 	xelatex -interaction=nonstopmode -jobname=composition_diagram_light $(TEX)

# # Compile inverted colors
# $(INVERTED_PDF): $(TEX)
# 	@echo "Compiling inverted color diagram..."
# 	xelatex -interaction=nonstopmode -jobname=composition_diagram_dark "\def\invertcolors{} \input{$(TEX)}"

# # Convert PDFs to SVGs using pdf2svg
# $(STANDARD_SVG): $(STANDARD_PDF)
# 	@echo "Converting standard PDF to SVG..."
# 	pdf2svg $(STANDARD_PDF) $(STANDARD_SVG)

# $(INVERTED_SVG): $(INVERTED_PDF)
# 	@echo "Converting inverted PDF to SVG..."
# 	pdf2svg $(INVERTED_PDF) $(INVERTED_SVG)

clean:
	@echo "Cleaning up extra generated files..."
	rm -f composition_diagram*.aux composition_diagram*.log composition_diagram*.out
	rm -f composition_diagram*.fdb_latexmk composition_diagram*.fls
	rm -f composition_diagram*.synctex.gz
	rm -f composition_diagram*.dvi
	rm -f composition_diagram*.pdf

allclean: clean
	@echo "Cleaning up generated SVGs..."
	rm -f composition_diagram*.svg
