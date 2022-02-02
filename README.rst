========================
MAPQ plugin
========================

This plugin provide wrappers around several programs of `MapQ <https://github.com/gregdp/mapq>`_.

+------------------+------------------+
| stable: |stable| | devel: | |devel| |
+------------------+------------------+

.. |stable| image:: http://scipion-test.cnb.csic.es:9980/badges/eman2_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/eman2_sdevel.svg


Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version (not avalaible yet)

.. code-block::

    scipion installp -p scipion-em-mapq

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-mapq.git

    * install

    .. code-block::

        scipion installp -p path_to_scipion-em-mapq --devel

MapQ binaries will be installed automatically with the plugin.

* Default installation path assumed is ``software/em/chimera-1.16``, if you want to change it, set *MAPQ_HOME* in ``scipion.conf`` file pointing to the folder where Chimera with MapQ is installed.

To check the installation, simply run one of the following Scipion tests:

.. code-block::

   TODO

A complete list of tests can also be seen by executing ``scipion test --show --grep mapq``

Supported versions
------------------

* 1.16

Protocols
---------

* compute q-scores

References
----------

1. Pintilie, G., Zhang, K., Su, Z. et al. Measurement of atom resolvability in cryo-EM maps with Q-scores. Nat Methods 17, 328–334 (2020). https://doi.org/10.1038/s41592-020-0731-1
