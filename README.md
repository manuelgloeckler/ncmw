[![DOI](https://zenodo.org/badge/377459650.svg)](https://zenodo.org/badge/latestdoi/323376678)

NCMW (Nasal community modeling workflow)
========================================

A metabolic modeling workflow for communities of organisms in an nasal
medium. This implements functionality to analyze the interaction of
metabolic models within a community.

The full documentation can be found on:
<https://manuelgloeckler.github.io/ncmw/>

Installation
------------

Please clone the repository using `git clone URL` . The package and all the dependencies can be installed via pip. Just enter:
:   `pip install -e ncmw`

when you are in the same directory as the pacakage. This will install
the python package `ncmw` but also console scripts `ncmw_setup` ,
`ncmw_analysis`, `ncmw_community` and `ncmw`, which calls the latter
three scripts in sequential order.

These should be added to your PATH, to try if this worked use::
:   ncmw --help ncmw\_setup --help ncmw\_analysis --help ncmw\_community
    --help

These scripts use [hydra](https://hydra.cc/docs/intro/), thus we can
overwrite all default values which are organized within `data/hydra`.
More on overwritting arguments in scripts\_

*NOTE* Automatically adding scripts to the PATH may not work if the
python environment is not setuped correctly. Alternatively you can call
the scripts by `python scripts/ncmw_setup`

Overview
--------

The package is organised into five parts
:   -   **setup\_models** - which setup metabolic models to the given
        default configuration and medium e.g. typically SNM3.
    -   **analysis** - which contains all function for single model
        analysis .
    -   **community** - which contains all code for community analysis
        and implemented community models.
    -   **visualization** - which contains all code to produce figures.

The **data** folder can be used to modify default values
:   -   **configs** contains the default settings for the maximum lower
        and upper bound of reactions, which is -1000 and 1000 .
    -   **hydra** contains all default parameters for the scripts, which
        can be overwritten in the console.
    -   **medium** contains the default *SNM3 medium* .
    -   **models** contains all models that are used by the scripts.

You can change the default values if needed or add additional models to
the workflow.

Quick Start
-----------

### Scripts

The scripts will automatically produce some results, that is:

    -   **ncmw\_setup**:
        :   -   Will set the default bounds as specified in
                data/configs, as well as the medium as specified in
                data/medium.
            -   Will gapfill the model with reactions or the medium with
                metabolites, such that all models obtain growth on the
                specified medium

    -   **ncmw\_analysis**
        :   -   Will perform flux variability analysis and visualize
                results on all exchange reactions.
            -   Will analyze the uptake/sekretion overlap between models
            -   Will compute the similarity of models, based on number
                of shared metabolites/reactions.
            -   Will compute the **COMPM medium**, which is the medium
                in which all models are able to obtain their maximal
                biomass rate.

    -   **ncmw\_community**
        :   -   Will create several kinds of community models.
            -   Will compute the **COOPM medium**, which is the smallest
                medium such that the community achives 10% of the
                maximal biomass rate, which induces cooperation.
            -   Will visualize the observed interactions between models.
            -   Will investigate the dependence of community weight and
                observed growth.

