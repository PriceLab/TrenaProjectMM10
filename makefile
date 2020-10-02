all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes trenaProjectMM10)

install:
	(cd ..; R CMD INSTALL --no-test-load trenaProjectMM10)

check:
	(cd ..; R CMD check `ls -t trenaProjectMM10) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

