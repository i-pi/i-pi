# Makefile for the help files
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http.//www.gnu.org/licenses/>.

.PHONY : all distclean clean aux
all : manual.pdf

aux:
	python create_man.py
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
