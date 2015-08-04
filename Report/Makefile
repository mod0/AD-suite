DOCNAME=report
DVIPSFLAGS = -f -tletter -Ppdf -G0 -O 0.0in,0.0in
SOURCES = $(wildcard *.tex) $(wildcard figures/*.pdf) $(wildcard figures/*.jpg) $(wildcard *.bib) $(wildcard *.c) $(shell find ../ -name 'README.tex')

all: $(DOCNAME).pdf 

$(DOCNAME).pdf: $(SOURCES) $(DOCNAME).bbl
	-pdflatex $(LATEXFLAGS) $(DOCNAME) 
	-bibtex $(DOCNAME) 
	-pdflatex $(LATEXFLAGS) $(DOCNAME) 
	-pdflatex $(LATEXFLAGS) $(DOCNAME)

$(DOCNAME).bbl: $(wildcard *.bib)
	-pdflatex $(LATEXFLAGS) $(DOCNAME)
	-bibtex $(DOCNAME) 
	-pdflatex $(LATEXFLAGS) $(DOCNAME)

final: all

clean:
	$(RM) $(DOCNAME).bbl $(DOCNAME).aux $(DOCNAME).dvi $(DOCNAME).blg $(DOCNAME).log $(DOCNAME).ps $(DOCNAME).pdf $(DOCNAME).toc $(DOCNAME).out *~

.PHONY: final
