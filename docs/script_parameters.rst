=================
Script parameters
=================

In this section we will describe the possibe parameters for the workflow scripts.

General parameters for any of the scripts:
* name: The name of the project. This will affect the folder name within the result folder, the default name is "default_project_name"
* seed: The random seed, which is used by any random operation (default = 1).
* eps: Thats a cutoff growth on which we perform e.g. gapfilling (default= 0.01)
* run_setup: Wheather to run setup (default true)
* run_analysis: Wheather to run analysis (default true)
* run_community: Wheather to run community (default true)
* overwrite_all: Wheather to overwrite all files (by default the scripts skip tasks which were already completed) (default false)

Now there are some parameters for the setup:

* fastcc: Wheather to improve model quality with FastCC (default true)
* set_bounds_and_medium: Wheather to set bound and medium as specified in "data" (default true)
* memote_evaluation: Wheather to run a memote evaluation. This may take a while! (default true)
* memote_solver_time_out: A time out for memote (it can get stuck some times...) (default 10 minutes)
* gapfill: Which type of gapfilling should be performed, we use: "medium" or "model" to either extend the medium or the model (default "medium")
* configs: The name of the config file (default default.json)
* medium: The name of the medium file (default snm3.json)
* models: The folder name of the models (default models)

Next there are some parameters for analysis:
* fva_fraction: The fraction of optimal growth that must be achived in FVA. (default 1.)
* sekretion_uptake: Wheather to compute sekretion and uptake (default true)
* require_setup: Wheather we require that setup was run before (default true)
* check_transport: Wheater we shoud check the transport reactions (default true)
* check_uptake_sekretion: Wheater we should check the uptakes and sekretion (default true)

Now also the community has several parameters which can be set:

* bag_of_reactions_model: Wheater we should perform all experiments with the BagOfReaction model (default true)
* shuttle_reactions_model: Wheater we should perform all experiments with the ShuttleReaction model (default true)
* shuttle_reaction_params:
  * shared_reactions: The set of reactions which are shared (experimental, not tested, default None -> all)
  * shuttle_reaction_lower_bound: The lower bound of shuttle reactions (default -50)
  * shuttle_reaction_upper_bound: The upper bound of shuttle reactions (default is 1000)
  * shuttle_regularization: Wheater to regularize shuttles to make them growth depndent (0 growth -> 0 input flux, default true)
  * shuttle_regularization_factor: Shuttle regularization factor, the higher the less influence (default 10000)
* main_weights: The weights considered for each experiment. You can either specify it by yourself as list i.e [1,1,1,1,1]. Equivalently to this you can use "ones". We also implement "fair" which tries to choose weights antiproportional to their individual growth values. (default "ones")
* models_folder: The folder for models that should be considered. If "setup" is specified it takes them from a previous setup run! (default "setup")
* save_models: Wheater to save models using pickle! (default true)
* save_as_sbml: Wheater to save models as sbml (xml) (default true)
* coopm_params:
  * fraction: The fraction of maximum biomass rate that the community have to attain minimally (default 0.1)
  * enforce_survival: Wheater to enforece the survival of **all** members (default  0.01)
* medium_strict_subset_of_default: Wheater the medium must be a strict subset of default medium i.e. excluding previous addition i.e. by gapfilling (default false)
* medium:
  * default: Perform all experiments on default medium (default true)
  * compm: Perform all experiments on COMPM medium (requires analysis run previous to that) (default true)
  * coopm: Perform all experiments on COOPM medium (requires COMPM) (default true)
* optimize_params:
  * enforce_survival: If for each optimzation run we should enforece the survival of 
* compute_community_fva: true
* community_fva_params:
  * fraction_of_optimum: 0.9
  * processes: 10
* cooperative_tradeoff: true
* cooperative_tradeoff_params:
  * alpha: 0.9
* pairwise_growth: true
* pairwise_growth_params:
  * h: 100
* compute_infer_weights: true
* infer_weights:
  * simulations_for_different_weights: 2000
  * observed_individual_biomass: balanced
  * medium: default
  * enforce_survival: 0.0
  * competitive_tradeoff: false
  * competitive_tradeoff_alpha: 0.9

Last but not least we add parameters for the visualization:
* scaled_medium_growth_plot:
  * min_scale: 0.1
  * max_scale: 100
  * evaluations: 1000
* cmap: tab20
* names:
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
