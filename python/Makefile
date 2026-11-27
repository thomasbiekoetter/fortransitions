# Build and install the Python package fortransitions.
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

PKG := fortransitions

# The libraries find each other through the runpath $ORIGIN (their own
# directory), so no patchelf step is needed. Inside a conda/pixi
# environment, the environment's lib/ is added so that its libgfortran
# is used instead of the system one.
RPATH := -Wl,-rpath,$$ORIGIN
ifeq ($(notdir $(FC)),lfortran)
# LFortran runs its link command through a shell, which would expand
# $ORIGIN as an (empty) variable. Escape it twice: fpm removes one level
# of backslashes, LFortran's shell the other.
RPATH := -Wl,-rpath,\\$$ORIGIN
endif
ifdef CONDA_PREFIX
RPATH += -Wl,-rpath,$(CONDA_PREFIX)/lib
endif

# One fpm build directory per toolchain (pixi environment or system
# compiler). fpm does not relink when only the link flags change, so with
# a shared build directory, switching between e.g. pixi and a system
# Python would install a library linked for the other environment.
ifdef CONDA_PREFIX
BUILD_DIR := build/$(notdir $(FC))-env-$(notdir $(CONDA_PREFIX))
else
BUILD_DIR := build/$(notdir $(FC))-system
endif

.PHONY: wrapper lib clean

# Strict editable mode: the package is reached through a plain sys.path
# entry (a link tree in build/), not through an import hook. With the
# default hook, `import fortransitions` in the directory that contains the
# repository would find the repository folder (also named fortransitions)
# as an empty namespace package instead.
wrapper: lib
	$(PC) -m pip install -e . --config-settings editable_mode=strict $(PIP_FLAGS)

# Always pass --profile: without it, fpm applies none of the flags in
# fpm.toml (for ifx, this includes the required -fp-model=precise).
lib:
	rm -rf $(PKG)/lib $(PKG)/include
	fpm install --profile release --compiler '$(FC)' --flag '-fPIC' \
		--link-flag '$(RPATH)' --build-dir '$(BUILD_DIR)' --prefix ./$(PKG)/

clean:
	-$(PC) -m pip uninstall -y $(PKG)
	rm -rf build $(PKG)/lib $(PKG)/include *.egg-info
