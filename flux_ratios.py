from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.core.Model import optimize
from copy import deepcopy
from draw_flux import DrawFlux
from html_writer import HtmlWriter
from analysis_toolbox import model_summary
from cobra.flux_analysis.single_deletion import single_reaction_deletion_fba 

model = create_cobra_model_from_sbml_file('data/ecoli_core.xml')

ko_list = ['PGM', 'ICL','G6PDH2r', 'PGL', 'MALS'] # ecoli with pfk (delta GPM, and not hexa)
knockouts = map(model.reactions.get_by_id, ko_list) 
model.remove_reactions(knockouts)

glucose_uptake = model.reactions.get_by_id('EX_glc_e')
pyruvate_uptake = model.reactions.get_by_id('EX_pyr_e')
xylose_uptake = model.reactions.get_by_id('EX_xyl-D_e')

glucose_uptake.lower_bound = 0
pyruvate_uptake.lower_bound = -10
xylose_uptake.lower_bound = -.001