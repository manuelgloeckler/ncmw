# NCMW (Nasal community modeling workflow)

A metabolic modeling workflow for communities of organisms in an nasal medium. Please clone the repository at into a place you want. The package and all the dependencies can be installed via pip. Just enter:
```ruby
    pip install ncmw
```

The package is organised into five parts:
* Setup - which setup metabolic models to the given default configuration and medium e.g. typically SNM3
* Analysis - which contains all function for single model analysis 
* Community - which contains all code for community analysis and implemented community models
* Visualization - which contains all code to produce figures
* Workflows - which contain scripts that is runable from the command line.

The data folder can be used as input, but also contains all default configurations, managed by hydra.
