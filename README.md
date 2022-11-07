# Solvers project
Systematically explore various solvers for 3d simulations

Testing and exploration ground for Solvers
Generate the site https://mwiesenberger.github.io/solvers/matrix-free.html

## Install
In order to generate the static website we use
[jupyter-book](https://jupyterbook.org).

```bash
# install jupyter-book
pip install -U jupyter-book
# clone this repository
git clone https://github.com/mwiesenberger/solvers
# build the book
jupyter-book build path/to/solvers
# we use ghp-import to publish changes on the github pages
pip install ghp-import
```
In order to locally generate the simulation data you will need the
[Feltor](https://github.com/feltor-dev/feltor) code repository.  Follow the
quick-start guide to install.  It is recommended to keep Feltor and this
repository next to each other.  If you prefer not to, you need to set the
`FELTOR_PATH` environment variable in order for the `Makefile` to
compile the code.

## Usage
Build the book with
```bash
jupyter-book build path/to/solvers
```
To publish changes after the book was built
```bash
cd path/to/solvers
ghp-import -n -f -p -o _build/html
```

## Author
Matthias Wiesenberger
