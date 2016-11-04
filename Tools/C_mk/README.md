The build system requires GNU `make` >= 3.81. 

Typically an application will have its own `GNUmakefile`.  (See
`Tutorials/HelloWorld_C/` for a simple example, or
`Tutorials/AMR_Adv_C_v2/Exec/SingleVortex/` for a lightly more
complicated example.)  `Make.defs` is included near the beginning, and
`Make.rules` is included in the end.  Depending the need,
`GNUmakefile` includes a number of
`$(BOXLIB_HOME)/Src/xxx/Make.package`, where `xxx` is `C_BaseLib` etc.
These `Make.package` files add sources to the make system.  An
application also has its own `Make.package`.  The make variables for
the sources are:

* `CEXE_sources` for C++ source files with `.cpp` extension.
* `CEXE_headers` for C++ headers with `.h` or `.H` extension.
* `cEXE_sources` for C source files with `.c` extension.
* `cEXE_headers` for C headers with `.h` or `.H` extension.
* `f90EXE_sources` for free format Fortran sources with `.f90` extension.
* `F90EXE_sources` for free format Fortran sources with `.F90` extension requiring preprocessing.
* `fEXE_sources` for fixed format Fortran sources with `.f` extension.
* `FEXE_sources` for fixed format Fortran sources with `.F` extension requiring preprocessing.

Typically one types `make` to build an executable.  Other useful
commands include:

* `make clean` removes the executable, `.o` files, and the resulting files of preprocessing.
* `make realclean` removes all files generated by the make system.
* `make help` shows the rules for compilation.
* `make print-xxx` shows the value of `xxx`.  This is very useful for
  debugging.
* `make file_locations` shows the path of each file.
* `make tags` and `make TAGS` generate tag files using `ctags` and `etags`, respectively.

The `Make.defs` includes the following files in the listed order:

* `Make.machines`: 

* `comps/xxx.mak`:

* `packages/Make.xxx`:

* `sites/Make.xxx`:

* `Make.local`:


VPATH_LOCATIONS INCLUDE_LOCATIONS libraries...

See User's Guide for more information on various build options (e.g.,
`USE_MPI`). 