all:
	make README.md

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif


README.md: README.Rmd
	"$(R_HOME)/bin/R" --vanilla -e "knitr::knit('README.Rmd');"

.PHONY: REAMDE.Rmd all
