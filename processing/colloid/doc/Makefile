all: manual.dvi

source = manual.tex \
		 introduction.tex \
		 normalmode.tex \
	     ellipsoid_normalmode.tex

manual.dvi: $(source)
	latex manual

pdf: manual.dvi
	dvipdfm manual
