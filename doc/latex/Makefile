# Makefile for the help files

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


.PHONY: all distclean clean aux
all: manual.pdf

aux:
	PYTHONPATH=../../:../scripts/:$(PYTHONPATH) python3 ../scripts/create_latex_manual.py
	rm -f *.pyc

manual.pdf: aux
	pdflatex manual
	bibtex manual
	pdflatex manual
	pdflatex manual

clean:
	bash -c "rm -rf input_docs manual.{aux,bbl,blg,brf,idx,log,lof,out,toc}"

distclean: clean
	rm -f manual.pdf manual.xml
