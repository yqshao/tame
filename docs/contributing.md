# Contributing

Contributions in the form of issues, bug reports and extension of
functionalities are welcome. The development process of TAME is tracked with
[Github Issues](https://github.com/teoroo-cmc/tame) and the long-term roadmap
 can be found at [Github
Project](https://github.com/orgs/teoroo-cmc/projects/2). If you would like to
request or implement new features, you are suggested to open an issue to discuss
it or contact the core developer (Yunqi Shao). Below are some resources if you
would like to help with the development of TAME.

**Code of conduct**: be respectful to each other.

## Getting started

The easiest way to try TAME is with its CodeSpace with pre-installed
dependencies ready for testing and development, the development environment can
be used in the browser without any installation, or with your own editor. Below
are the instructions for those wish to set up a development environment
**locally**.

### Environment preparation

To clone the repository and install the necessary dependencies:

``` bash
git clone
# it's recommended to separate your development environments, 
# here we use the most readily available venv
python venv ~/
pip install
```

You can verify the installation by runninng the unit test:

```bash
pytest 
```

### Editor setup with VS Code

When contributing to the code please conform to the coding conventions (below)
to minimize friction. Those are automatically set up if you edit the source code
using the VS Code editor. The configuration includes:

- Code formatting for Python and Markdown;
- Local service for documentation building and preview.

To make sure that VSCode recognize the virtual environment you created, open the
command palette (Ctrl-Shift-P or F1); type and select
`Python: Select Interpreter` and find your environment. Doing so makes sure that
whenever you open the project, the correct Python is used.

## Coding conventions

### Dependencies

Basic requirements (numpy, click, etc) are specified in `setup.py`, version
numbers could be given as the `>=old.ver.that.works`, and both old and new
versions should be tested in the CI process.

Development and documentation building requirements are fixed with exact version
numbers. Versions may be changed if new features becomes desirable, bugs are
found, or infrequently between version changes.

### Formatting

Python code, both in source and documentation/notebooks, should be formatted
with [black](https://github.com/psf/black) by default. This will be checked upon
commit and pull requests. Markdown for documentation should be formatted
with the [markdownlint](https://github.com/DavidAnson/markdownlint).

### Documentation

The Python classes and functions should be documentated according to the [Google
stlye
guide](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings).
New recipes should be supplied with corresponding references and equations
whenever necessary.

### Versioning

Versioning is determined according to the [PEP
440](https://peps.python.org/pep-0440/) convention. Specifically in the form:
`X.Y.Z(.suffix)`, where X is the major release, Y is the minor release and Z is
the micro release. X is reserved for major milestones in the project roadmap, Y
for the implementation of new features, and Z for the small bug fixes.

The main branch should be named as X.Y.Z.dev, after the next release number.
Developmental or pre releases may be tagged for internal usage (in connection
with publication). `CHANGELOG.md` should be updated between version changes.
