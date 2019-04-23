all:  docs install

docs:
	R -e "devtools::document()"
vig:
	R -e "devtools::build_vignettes()"

build:
	(cd ..; R CMD build --no-build-vignettes TrenaProjectHG38)

install:
	(cd ..; R CMD INSTALL --no-test-load TrenaProjectHG38)

check:
	(cd ..; R CMD check `ls -t TrenaProjectHG38_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t TrenaProjectHG38_* | head -1`)


test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

