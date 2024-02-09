These are scripts used for modeling in the paper:
A Scube2-Shh feedback loop links morphogen release and spread to morphogen signaling to enable scale invariant patterning of the ventral neural tube
Zachary M. Collins, Anna Cha, Albert Qin, Kana Ishimatsu, Tony Y.C. Tsai, Ian A. Swinburne, Pulin Li, Sean G. Megason
bioRxiv 469239; doi: https://doi.org/10.1101/469239

These files are written by Albert Qin and Pulin Li

here's what the various scripts do:

+metrics: specifies various metrics 

+models: the actual model equations are here; model is a base class from which HS (first-order degradation) and HS_self_enhanced (self-enhanced degradation) inherit

+utils: miscellaneous stuff I downloaded to make progress bars and animations and the like

run_model: a function that takes in a model and parameters as input and spits out the dynamics of all species

gridsearch_final: run this one to run *all* the gridsearches (the run results for all the searches are stored in the folder "gridsearch_results")

diffusion_vs_release: generates the figure that compares diffusion to diffusion + release for different values of hedgehog flux

dot_plots: generates the dot_plots

phase_plot: generates the D on kappa phase plot

degradation_phase_plot: generates a gamma on kappa phase plot

draw_SS_gradients: draws six plots: the ligand concentration on absolute / relative position for a model that does both release+diffusion, release only, and nothing (I used this to generate the example plots of good, medium, bad scaling)

plot_scaling_dynamics: generates three gifs that plot the ligand, scube, and scube expression dynamics
