.PHONY: test

test:
	@# See https://github.com/pytest-dev/pytest/issues/7650 -- this can help
	@# avoid weird issues due to multiple pytest versions
	python3 -m pytest notebooks/
