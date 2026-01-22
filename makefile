##################################################################
# MultiNEAs Makefile
##################################################################

.PHONY: help install install-dev test clean build upload upload-test docs push release
RELMODE=release
PYTHON ?= python3
COMMIT_MSG ?= chore: sync tracked changes

help:
	@echo "MultiNEAs Development Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  install      - Install the package"
	@echo "  install-dev  - Install package in development mode with dev dependencies"
	@echo "  test         - Run tests with pytest"
	@echo "  clean        - Remove build artifacts and cache files"
	@echo "  build        - Build distribution packages"
	@echo "  upload       - Upload package to PyPI"
	@echo "  upload-test  - Upload package to TestPyPI"
	@echo "  docs         - Build documentation (installs docs requirements)"
	@echo "  push         - Commit (all files) and push current branch"
	@echo "  release      - Release a new version (usage: make release RELMODE=release VERSION=x.y.z)"

install:
	pip install .

install-dev:
	pip install -e .
	pip install -r requirements-dev.txt

test:
	pytest

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf src/*.egg-info
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete

build: clean
	$(PYTHON) -m build

upload: build
	$(PYTHON) -m twine upload dist/*

upload-test: build
	$(PYTHON) -m twine upload --repository testpypi dist/*

docs:
	$(PYTHON) -m pip install -r docs/requirements.txt
	cd docs && $(PYTHON) -m sphinx.cmd.build -M html "." "_build"

push:
	@echo "Committing tracked changes (if any)..."
	@if ! git diff --quiet || ! git diff --cached --quiet || [ -n "$$(git status --porcelain)" ]; then \
		git add . && \
		git commit -m "$(COMMIT_MSG)"; \
	else \
		echo "Working tree is clean (tracked files); nothing to commit."; \
	fi
	@echo "Pushing current branch..."
	@git push -u origin HEAD

# Example: make release RELMODE=release VERSION=0.2.0.2
release: push
	@test -n "$(VERSION)" || (echo "ERROR: VERSION is required. Example: make release RELMODE=release VERSION=0.2.0" && exit 1)
	@echo "Releasing a new version..."
	@bash bin/release.sh $(RELMODE) $(VERSION)
