all: python protocol.pdf

protocol.pdf: folders

python: folders

protocol.pdf: FORCE | build
	cd latex-template/ && \
	pdflatex --output-directory "../build/tex" main.tex
	cp "./build/tex/main.pdf" "protocol.pdf"

python: *.py python_lib/*.py
	python main.py


init: purge
	make folders
	cp -r template/text text
	cp -r template/python/* .

clean:
	rm -f protocol.pdf
	rm -rf build

purge: clean
	echo "Do you really want to destroy all work and start over fresh? Press enter to continue."
	read -n 1 no
	rm -rf text
	rm -f *.py

folders:
	mkdir -p build
	mkdir -p build/tables
	mkdir -p build/plots
	mkdir -p build/tex

FORCE:

.PHONY: all clean
