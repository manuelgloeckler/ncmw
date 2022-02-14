=================
Script parameters
=================

In this section, we will describe the possible parameters for the workflow scripts. Further we will quickly introduce how to change parameters.

The default parameters are listed in `data/hydra` within .yml files. To run the whole workflow just open your terminal of choice and write

.. code-block:: ruby

   Some Ruby code.


It will lsit all parameters at the beginning, you can verify that these are exactly the on listed within the files. Thus one way to modify the parameters is to simply modify these files!

Yet, in some cases we may not want to permanetly change a default value, but just for one run. In this case we can overwrite the parameters within the terminal. For example to create a new project we can just modify the command using

.. code-block::

   ncmw name="other_name"

Note that `name` is one of the general parameters and thus can be set directly. Lets also consider how to change parmeters from setup or community. Lets say we want to disable "fastcc" because our models are already curated and further do not want to produce results based on the BagOfReactions Community model, then we can use 

.. code-block::

   ncmw setup.fastcc=false community.bag_of_reactions_model=false

We can also directly add some new parameters, this maybe usefull to extend the "names" i.e 

.. code-block::

   ncmw '+visualization.names={"id":"name"}'

General parameters for any of the scripts:

* **name**: The name of the project. This will affect the folder name within the result folder, the default name is "default_project_name"
* **seed**: The random seed, which is used by any random operation (default = 1).
* **eps**: That's a cutoff growth on which we perform e.g. gap filling (default= 0.01)
* **run_setup**: Wheather to run setup (default true)
* **run_analysis**: Wheather to run analysis (default true)
* **run_community**: Wheather to run community (default true)
* **overwrite_all**: Wheather to overwrite all files (by default the scripts skip tasks which were already completed) (default false)

Now there are some parameters for the setup:

* **fastcc**: Wheather to improve model quality with FastCC (default true)
* **set_bounds_and_medium**: Wheather to set bound and medium as specified in "data" (default true)
* **memote_evaluation**: Wheather to run a memote evaluation. This may take a while! (default true)
* **memote_solver_time_out**: A time out for memote (it can get stuck some times...) (default 10 minutes)
* **gapfill**: Which type of gap filling should be performed, we use: "medium" or "model" to either extend the medium or the model (default "medium")
* **configs**: The name of the config file (default default.json)
* **medium**: The name of the medium file (default snm3.json)
* **models**: The folder name of the models (default models)

Next, there are some parameters for analysis:

* **fva_fraction**: The fraction of optimal growth that must be achieved in FVA. (default 1.)
* **sekretion_uptake**: Wheather to compute secretion and uptake (default true)
* **require_setup**: Wheather we require that setup was run before (default true)
* **check_transport**: Wheater we should check the transport reactions (default true)
* **check_uptake_sekretion**: Wheater we should check the uptakes and secretion (default true)

Now also the community has several parameters which can be set:

* **bag_of_reactions_model**: Wheater we should perform all experiments with the BagOfReaction model (default true)
* **shuttle_reactions_model**: Wheater we should perform all experiments with the ShuttleReaction model (default true)
* **shuttle_reaction_params**:

  * shared_reactions: The set of reactions which are shared (experimental, not tested, default None -> all)
  * shuttle_reaction_lower_bound: The lower bound of shuttle reactions (default -50)
  * shuttle_reaction_upper_bound: The upper bound of shuttle reactions (default is 1000)
  * shuttle_regularization: Wheater to regularize shuttles to make them growth depndent (0 growth -> 0 input flux, default true)
  * shuttle_regularization_factor: Shuttle regularization factor, the higher the less influence (default 10000)
* **main_weights**: The weights considered for each experiment. You can either specify it by yourself as list i.e [1,1,1,1,1]. Equivalently to this you can use "ones". We also implement "fair" which tries to choose weights antiproportional to their individual growth values. (default "ones")
* **models_folder**: The folder for models that should be considered. If "setup" is specified it takes them from a previous setup run! (default "setup")
* **save_models**: Wheater to save models using pickle! (default true)
* **save_as_sbml**: Wheater to save models as sbml (xml) (default true)
* **coopm_params**:

  * fraction: The fraction of maximum biomass rate that the community have to attain minimally (default 0.1)
  * enforce_survival: Wheater to enforece the survival of **all** members (default  0.01)
* **medium_strict_subset_of_default**: Wheater the medium must be a strict subset of default medium i.e. excluding previous addition i.e. by gapfilling (default false)
* **medium**:

  * default: Perform all experiments on default medium (default true)
  * compm: Perform all experiments on COMPM medium (requires analysis run previous to that) (default true)
  * coopm: Perform all experiments on COOPM medium (requires COMPM) (default true)
* **optimize_params**:

  * enforce_survival: If for each optimzation run we should enforece the survival of 
* **compute_community_fva**: Wheater to compute FVA results for the community model (default true)
* **community_fva_params**:

  * fraction_of_optimum: Which fraction from MBR used in FBA. (default 0.9)
  * processes: Number of processors involved (default 10)
* **cooperative_tradeoff**: Wheather to compute results with cooperative_tradeoff (default true).
* **cooperative_tradeoff_params**:

  * alpha: The tradeoff value used (defualt 0.9)
* **pairwise_growth**: Wheater to compute all pairwise growth relationships (default true).
* **pairwise_growth_params**:

  * h: Number of discretiztaion steps for weights (default 100)
* **compute_infer_weights**: Wheater to infer weights (default true)
* **infer_weights**:

  * simulations_for_different_weights: Simulation used within the simulation-based inference procedure (default 2000)
  * observed_individual_biomass: Observed individual biomass values. As default we use "balanced" i.e. each species should have the same biomass value. If you have any custom inference goal you can set it here by passing a list of values.
  * medium: On which medium we do inference (default is on default medium)
  * enforce_survival: Wheater we use enforce survival constraints (default 0.0)
  * competitive_tradeoff: Wheater we use competitive_tradeoff (default false)
  * competitive_tradeoff_alpha: Parameterse used in competitive tradeoff (default 0.9)

Last but not least we add parameters for the visualization:

* **scaled_medium_growth_plot**:

  * min_scale: The minimal medium scale value (default 0.1)
  * max_scale: The maximum medium scale value (default 100)
  * evaluations:The evaluations (discretization steps) (default 1000)
* **cmap**: A color map used to choose colors!
* **names**: Here you can pass a map that maps model_id -> name. This will only change the name within figures!

  * CarveMe_SNM_gapfilled_model: M. catarrhalis
  * DP_83VPs_KB5: D. pigrum
  * iYS854: S. aureus
  * Staphylococcus_epidermidis_ATCC_12228: S. epidermidis
  * himodel: H. influenzae
  * Aba: A. baumannii
  * iDPM21RW: D. pigrum
  * MODEL1507180054: K. pneumoniae
  * Hin: H. influenzae
  * Slu: S. lugdunensis
