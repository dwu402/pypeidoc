src_name = pypei_variants
SHELL=/bin/bash

all:
	pandoc --from markdown --to latex ${src_name}.md --output ${src_name}.pdf

.PHONY: clean
clean:
	rm --force ${src_name}.pdf
	rm --force --dir ./tex2pdf.*/