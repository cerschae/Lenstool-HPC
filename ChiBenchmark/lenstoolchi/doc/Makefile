# Makefile for the doc directory and subdirectories

all: clean lenstool.ps lenstool.pdf

SOURCES = lenstool.tex run.tex inputfile.tex datafile.tex getting_started.tex example_par.tex 

lenstool.dvi: $(SOURCES)
	latex lenstool
	latex lenstool

lenstool.ps: lenstool.dvi
	dvips -f -Pcmz -o $@ $^

lenstool.pdf: lenstool.dvi
	dvipdf $^ $@ 

clean: 
	rm -f *.aux *.toc *.log *.dvi lenstool.ps
