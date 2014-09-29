from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from draw_flux import DrawFlux
from html_writer import HtmlWriter
from analysis_toolbox import model_summary
from copy import deepcopy
from cobra.flux_analysis.single_deletion import single_reaction_deletion_fba 
from cobra.flux_analysis.double_deletion import double_reaction_deletion_fba
from cobra.core import Reaction
from models import knockin_reactions,init_wt_model, knockout_reactions
from optknock import OptKnock

model = init_wt_model('full', carbon_sources={})

ki = 'PRK+RBC,RED'     # add prk+rubisco and a "free" redox reaction 
ko = 'PGM,ICL,G6PDH2r,PFK' # tetra deletions
ko = ko + ',PFL,POR5' # "pyruvate formate lyase" and pyruvate synthase 
ko = ko + ',EDA,PGCD,MGSA,GLYCK,GLYCK2,DRPA' # gpm breakers:
                                            # ED pathway
                                            # serine pathway 
                                            # methylglyoxal
                                            # glycerate kinase1
                                            # glycerate kinase2
                                            # deoxyribose-phosphate aldolase (Alternate Carbon Metabolism)
    
                                        # deoxyribose-phosphate aldolase, 2-methylcitrate dehydratase

knockin_reactions(model,ki)
knockout_reactions(model,ko)

optimize_minimal_flux(model)
print "solution: %.01f 1/h" % model.solution.f
html = HtmlWriter('res/report.html')
df = DrawFlux('EcoliMetabolism.svg')
df.ToSVG(model, model.solution, html)
model_summary(model, model.solution, html)


background_gr_dict, background_status_dict = single_reaction_deletion_fba(model, model.reactions)




#non_essential_reactions = [r for r,v in background_gr_dict.iteritems() if v>1e-5 and r not in ko.split(',')]
#
#optimize_minimal_flux(model)
#print "solution: %.01f 1/h" % model.solution.f
#
#growth_rate_dict, status_dict = single_reaction_deletion_fba(model, non_essential_reactions)
#essential_reactions = {r:model.reactions.get_by_id(r).subsystem for r,v in growth_rate_dict.iteritems() if v<1e-5}

