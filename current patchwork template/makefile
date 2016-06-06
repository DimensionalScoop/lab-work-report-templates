all: python Abgabe.pdf

spelling-mistakes: tex/content/*.tex
	cd tex/content/ && bash ../../spell-checker.sh

Abgabe.pdf: tex/main.tex tex/preamble.tex tex/content/* plots

Abgabe.pdf: FORCE | build
	cd tex/ && \
	pdflatex --output-directory "../build/" main.tex
	cp "./build/main.pdf" "Abgabe.pdf"

python: python/*.py
	cd python/ && \
	python main.py


clean:
	rm -f Abgabe.pdf
	rm -f spelling-mistakes
	rm -f build/*

FORCE:

.PHONY: all clean
