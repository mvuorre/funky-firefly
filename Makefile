all: renv docs

renv: renv.lock
	Rscript -e "renv::restore()"

docs: index.qmd
	quarto render

clean:
	rm -rf *.pdf *.typ *.png *_cache/ *_files/ docs/

clean-hard:
	rm -rf *.pdf *.typ *.png *_cache/ *_files/ docs/ cache/

.PHONY: all clean clean-hard
