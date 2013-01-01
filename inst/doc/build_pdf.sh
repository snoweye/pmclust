#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex pmclust-guide.Rnw
bibtex pmclust-guide
pdflatex pmclust-guide.Rnw
pdflatex pmclust-guide.Rnw
pdflatex pmclust-guide.Rnw
rm *.aux *.bbl *.blg *.log *.out *.toc
