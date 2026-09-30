# Define variables
TEX=composition_diagram.tex
STANDARD_DVI=composition_diagram.dvi
STANDARD_SVG=composition_diagram.svg

.PHONY: all clean allclean

all: $(STANDARD_SVG) clean

$(STANDARD_DVI): $(TEX)
	@echo "Compiling diagram to DVI..."
	lualatex -output-format=dvi $(TEX)

$(STANDARD_SVG): $(STANDARD_DVI)
	@echo "Converting DVI to SVG..."
	dvisvgm --libgs=/opt/homebrew/lib/libgs.dylib -o $(STANDARD_SVG) $(STANDARD_DVI)

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
