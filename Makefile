.PHONY: yml

yml:
	conda env export > conda-environment.yml
