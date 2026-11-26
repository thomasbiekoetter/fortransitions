# Build and install the Python package pyfortransitions.
#
#   make wrapper            build the shared library and pip-install the
#                           package (editable) into the active Python
#   make wrapper FC=ifx     the same with the Intel compiler
#   make lib                only build the shared library
#   make clean              uninstall and remove all build products
#
# The pixi tasks (pixi.toml) call these targets as well.

# make predefines FC = f77 and PC = pc, so ?= would not apply to them;
# replace only those built-in defaults (values from the command line or
# the environment are kept).
ifeq ($(origin FC),default)
FC := gfortran
endif
ifeq ($(origin PC),default)
PC := python
endif
PIP_FLAGS ?=

PKG := pyfortransitions

# The libraries find each other through the runpath $ORIGIN (their own
# directory), so no patchelf step is needed. Inside a conda/pixi
# environment, the environment's lib/ is added so that its libgfortran
# is used instead of the system one.
RPATH := -Wl,-rpath,$$ORIGIN
ifdef CONDA_PREFIX
RPATH += -Wl,-rpath,$(CONDA_PREFIX)/lib
endif

.PHONY: wrapper lib clean

wrapper: lib
	$(PC) -m pip install -e . $(PIP_FLAGS)

# Always pass --profile: without it, fpm applies none of the flags in
# fpm.toml (for ifx, this includes the required -fp-model=precise).
lib:
	fpm install --profile release --compiler '$(FC)' --flag '-fPIC' \
		--link-flag '$(RPATH)' --prefix ./$(PKG)/

clean:
	-$(PC) -m pip uninstall -y $(PKG)
	rm -rf build $(PKG)/lib $(PKG)/include *.egg-info
