.PHONY: lint, pretty
lint: 
	flake8 --select=F --ignore= --ignore=F403,F405 --per-file-ignores=**__init__.py:F401 ./
	black --check ./
pretty:
	black ./

