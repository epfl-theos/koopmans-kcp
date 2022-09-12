============
koopmans-kcp
============

| |GH Actions| |GPL License|

An implementation of Koopmans functionals in ``Quantum ESPRESSO v4.1`` with full orbital optimization.

Instead of using this as a stand-alone code, it is strongly recommended you install and use the larger `koopmans package <https://github.com/epfl-theos/koopmans>`_.

Installation
------------

Installation is a two-step procedure: first you must run ``configure``:

.. code-block:: bash

   ./configure MPIF90=<mpif90>

where ``<mpif90>`` should be replaced by the name of your chosen MPI Fortran90 compiler e.g. ``MPIF90=mpiifort``. The code should automatically detect and link the requisite libraries. It it does not, edit ``make.sys`` as desired.

Then you compile the code itself:

.. code-block:: bash

   make kcp

Running
-------
Calculations can be run with the command

.. code-block:: bash

   kcp.x <seed>.cpi

where <seed>.cpi is the ``kcp`` input file. However, it is strongly recommended that you do not run these commands directly but instead via the ``koopmans`` package. See https://koopmans-functionals.org for more details.

Contact
-------
Written and maintained by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna (2020-)

For help and feedback email edward.linscott@gmail.com

.. |GH Actions| image:: https://img.shields.io/github/workflow/status/epfl-theos/koopmans-kcp/Make/master?label=master&logo=github
   :target: https://github.com/epfl-theos/koopmans-kcp/actions?query=branch%3Amaster
.. |GPL License| image:: https://img.shields.io/badge/license-GPL-blue
   :target: https://github.com/epfl-theos/koopmans-kcp/blob/master/LICENSE

