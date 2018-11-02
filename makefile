all:  docs install

docs:
	R -e "devtools::document()"
vig:
	R -e "devtools::build_vignettes()"

build:
	(cd ..; R CMD build --no-build-vignettes TrenaProject)

install:
	(cd ..; R CMD INSTALL TrenaProject)

check:
	(cd ..; R CMD check `ls -t TrenaProject_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t TrenaProject_* | head -1`)


test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

