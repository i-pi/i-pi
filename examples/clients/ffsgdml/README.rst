*******************
sGDML Documentation
*******************

This is a highly optimized implementation of the recently proposed symmetric gradient domain machine learning (sGDML) force field model [#gdml]_ [#sgdml]_. It is able to faithfully reproduce detailed global potential energy surfaces (PES) for small- and medium-sized molecules from a limited number of user-provided reference calculations [#analysis]_. 

We provide a set of Python routines to reconstruct and evaluate custom sGDML force fields [#soft]_. A user-friendly command-line interface offers assistance through the complete process of model creation, in an effort to make this novel machine learning approach accessible to broad practitioners.

It's easy to get going!
-----------------------

.. _installation:

Installation
============


Stable release (recommended)
----------------------------

Most systems come with the default package manager for Python ``pip`` already preinstalled. Installing ``sgdml`` is as simple as calling:

.. code-block:: console

	$ pip install sgdml
	
The ``sgdml`` command-line interface can now be used from anywhere on the system:

.. code-block:: console

	$ sgdml train datasets/npz/ethanol.npz 200 1000 5000

Similarly, the ``sgdml`` Python API is available to other programs across the whole system:

.. code:: python

	from sgdml.predict import GDMLPredict
	gdml_predict = GDMLPredict(<model>)
  
  
Download pretrained sGDML models
----------------------------

To run i-Pi using the sGDML force field on pretained models, you can download them via: 

.. code-block:: bash

	$ sgdml-get model
  
Generating you own model is a straightforward procedure with the sGDML package. 
A detailed description of how to do it can be found in: http://sgdml.org/doc/


Code development
----------------

The sGDML code is developed through our GitHub repository: https://github.com/stefanch/sGDML

Use i-Pi with sGDML
----------------

The sGDML force field can be deployed by just specifying the name of the sGDML model in the input.xml file (e.g. benzene.DFT.PBE-TS.npz):

.. code-block:: xml

  ...
  <ffsgdml name='sgdml' pbc='False'>
    <sGDML_model> benzene.DFT.PBE-TS.npz </sGDML_model>
  </ffsgdml>
  ...
  <system>
    ...
    <forces>
      <force forcefield='sgdml'> </force>
    </forces>
    ...
  </system>
  ...

Citing
======

Please cite GDML and sGDML as follows:

.. [#gdml] Chmiela, S., Tkatchenko, A., Sauceda, H. E., Poltavsky, Igor, Schütt, K. T., Müller, K.-R. (2017). `Machine Learning of Accurate Energy-conserving Molecular Force Fields <http://advances.sciencemag.org/content/3/5/e1603015>`_. *Sci. Adv.*, **3(5)**, e1603015.
.. [#sgdml] Chmiela, S., Sauceda, H. E., Müller, K.-R., Tkatchenko, A. (2018). `Towards Exact Molecular Dynamics Simulations with Machine-Learned Force Fields <https://www.nature.com/articles/s41467-018-06169-2>`_. *Nat. Commun.*, **9(1)**, 3887.
.. [#soft] Chmiela, S., Sauceda, H. E., Poltavsky, Igor, Müller, K.-R., Tkatchenko, A. (2019). `sGDML: Constructing Accurate and Data Efficient Molecular Force Fields Using Machine Learning <https://doi.org/10.1016/j.cpc.2019.02.007>`_. *Comput. Phys. Commun.*, **240**, 38-45.
.. [#analysis] Sauceda, H. E., Chmiela, S., Poltavsky, Igor, Müller, K.-R., Tkatchenko, A. (2019). `Molecular Force Fields with Gradient-Domain Machine Learning: Construction and Application to Dynamics of Small Molecules with Coupled Cluster Forces <https://doi.org/10.1016/j.cpc.2019.02.007>`_. *J. Chem. Phys.*, **150**, 114102.