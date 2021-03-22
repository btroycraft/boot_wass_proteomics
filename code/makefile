.PHONY: all

EMAIL=
SERVER=

R = R --vanilla --file=rscripts/$(1) > rout/$(1).out
JULIA = julia --file=rscripts/$(1) > rout/$(1).out
MAIL = $(if $(EMAIL), mail -s $(1) ${EMAIL} <<< "$(1) has finished running on ${SERVER}.")

all:
clear_logs:
	rm rout/*.out

rscripts/plots.R: rscripts/gen_data.R rscripts/funcs.R

pbetti/pbetti_mean.csv: rscripts/pbetti_mean.R
	$(call R,pbetti_mean.R)
	$(call MAIL,pbetti_mean.R)

plots/plot_bootstrap_comp.pdf: rscripts/plots.R
	$(call R,plots.R)
