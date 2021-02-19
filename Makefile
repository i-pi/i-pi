.PHONY: lint, pretty
pretty:
	black ./
lint: 
	flake8 --select=F --ignore= --ignore=F403,F405 --per-file-ignores=**__init__.py:F401 ./
	black --check ./
