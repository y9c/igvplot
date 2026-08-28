# igvplot developer tasks
.PHONY: install test lint build gallery clean release-help

VENV := .venv
PY := $(VENV)/bin/python

install:            ## create venv + editable install
	python -m venv $(VENV)
	$(PY) -m pip install -U pip
	$(PY) -m pip install -e ".[dev]"

test:               ## run the test suite (warnings-as-errors via pyproject filterwarnings)
	$(PY) -m pytest -q

lint:               ## ruff check
	ruff check igvplot scripts examples

build:              ## build wheel + sdist
	$(PY) -m build

gallery:            ## regenerate the README example images
	$(PY) examples/generate_gallery.py

clean:              ## remove build artifacts
	rm -rf build dist *.egg-info .pytest_cache

release-help:       ## show release commands
	@echo "python scripts/bump_version.py --minor --tag && git push && git push --tags"
