from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.core.Model import optimize
from copy import deepcopy
from draw_flux import DrawFlux
from html_writer import HtmlWriter
from analysis_toolbox import model_summary
from cobra.flux_analysis.single_deletion import single_reaction_deletion_fba 

model = create_cobra_model_from_sbml_file('data/iJO1366.xml')

ko_list = ['PGM', 'ICL','G6PDH2r', 'PGL', 'MALS'] # ecoli with pfk (delta GPM, and not hexa)
knockouts = map(model.reactions.get_by_id, ko_list) 
model.remove_reactions(knockouts)

eat_glucose = model.reactions.get_by_id('EX_glc_e')
eat_glucose.lower_bound = 0

eat_pyruvate = model.reactions.get_by_id('EX_pyr_e')
eat_pyruvate.lower_bound = -10

growth_rate_dict, status_dict = single_reaction_deletion_fba(model, model.reactions)

potential_bypasses = set([k for k,v in growth_rate_dict.iteritems() if v==0])

model.add_reactions(knockouts)
growth_rate_dict, status_dict = single_reaction_deletion_fba(model, model.reactions)
essential_genes = set([k for k,v in growth_rate_dict.iteritems() if v==0])

bypasses = potential_bypasses - (potential_bypasses & essential_genes)
bypasses = map(model.reactions.get_by_id, bypasses)
genes = map(lambda x: x.genes, bypasses)
reactions = map(lambda x: x.name, bypasses)

gene_dictionary = {row[0:5]:row[84:].split(';')[0].strip() 
                            for row in open("data/all_ecoli_genes.txt", 'r')}

gene_names = [map(lambda x: gene_dictionary[x.id], x) for x in genes]
output = zip(reactions, gene_names)


for r in bypasses:
    old_lb = r.lower_bound
    old_ub = r.upper_bound
    r.lower_bound = 0
    r.upper_bound = 0
    solution = optimize(model)

    html = HtmlWriter('res/pyruvate_leaks_%s.html' %r.name)
    df = DrawFlux('EcoliMetabolism.svg')
    df.ToSVG(model, solution, html)
    model_summary(model, solution, html)

    html.close
    
    r.lower_bound = old_lb
    r.upper_bound = old_ub


