========================================
NCMW (Nasal community modeling workflow)
========================================
A metabolic modeling workflow for communities of organisms in an nasal medium. This implements functionality to analyze the interaction of metabolic models within a community. 


Installation
============

Please clone the repository using ``git clone URL`` . The package and all the dependencies can be installed via pip. Just enter:
    ``pip install ncmw``
when you are in the same directory as the pacakage. This will install the python package  ``ncmw`` but also console scripts ``ncmw_setup`` , ``ncmw_analysis`` and ``ncmw_community`` .
These should be added to your PATH, to try if this worked use::
    ncmw_setup --help
    ncmw_analysis --help
    ncmw_community --help

These scripts use `hydra <https://hydra.cc/docs/intro/>`__, thus we can overwrite all default values which are organized within ``data/hydra``. More on overwritting arguments in `scripts`_

*NOTE* Automatically adding scripts to the PATH may not work if the python environment is not setuped correctly. Alternatively you can call the scripts by  ``python scripts/ncmw_setup``


Overview
========

The package is organised into five parts
    * **setup_models** - which setup metabolic models to the given default configuration and medium e.g. typically SNM3.
    * **analysis** - which contains all function for single model analysis .
    * **community** - which contains all code for community analysis and implemented community models.
    * **visualization** - which contains all code to produce figures.
    
The **data** folder can be used to modify default values
    * **configs** contains the default settings for the maximum lower and upper bound of reactions, which is -1000 and 1000 .
    * **hydra** contains all default parameters for the scripts, which can be overwritten in the console. 
    * **medium** contains the default *SNM3 medium* .
    * **models** contains all models that are used by the scripts. 
You can change the default values if needed or add additional models to the workflow.


Quick Start
===========

Scripts
-------
.. _scripts:
asdf


Results
=======
öjölklkj