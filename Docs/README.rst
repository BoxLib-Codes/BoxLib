Building BoxLib Documentation
=============================

This directory contains, in addition to the BoxLib users'
gude, the framework for building documentation that leverages 
both embedded in-source documentation as well as additional
notes and examples.

The mechanics are that the source code is processed with Doxygen, 
which in addition to outputting it's own html with a complete 
description of documented functions, class hierarchy diagrams, etc,
outputs XML that can be referenced by the primary Sphinx-based 
documentation by using Breathe. The Sphinx documentation includes
additional examples and working notes on a curated subset of the 
functions available.

The build process is therefore two steps::

    Doxygen doxyfile
    cd sphinx_doc
    make html

Dependencies
------------


#. `Doxygen <http://www.doxygen.org>`_.
#. `Sphinx <http://sphinx-doc.org>`_.
#. `Breathe <http://breathe.readthedocs.io>`_.


