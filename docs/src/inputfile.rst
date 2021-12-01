Input File
==========

SHAW's input file uses the yaml format,
and is organized into sections: *general*, *io*, *source*, *material*.

At a high-level, it looks something like:

.. code-block:: yaml

  general:
   # ...

  io:
   # ...

  source:
   # ...

  material:
   # ...


.. Attention::

  A valid input file must contain a **single** instance of **each** of these sections.

Follow the following links to see the details for each section:

.. toctree::
    :maxdepth: 2

    inputfile_description

A full template of the input file is available here.
You can use it as a starting point.

.. toctree::

    inputfile_template
