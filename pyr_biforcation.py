from cobra.io.sbml import create_cobra_model_from_sbml_file
from html_writer import HtmlWriter
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from draw_flux import DrawFlux
from analysis_toolbox import model_summary
from models import *
import pandas as pd
from cobra.core import Reaction, Metabolite, Formula

def draw_flux(html_fname):
    html = HtmlWriter(html_fname)
    df = DrawFlux('CentralMetabolism.svg')
    df.ToSVG(model, model.solution, html)
    model_summary(model, model.solution, html)
    html.close

def add_metabolite(model, id, formula='', name='', compartment='C'):
        met = Metabolite(id=id, formula=formula, 
                         name=name, compartment=compartment)
        model.add_metabolites([met])
        
def add_reaction(model, id, name, sparse,
                 lower_bound=0, upper_bound=1000):
    """
        Adds a reaction to the model
    """
    # convert the sparse representation using the metabolites in the model
    for key in sparse.keys():
        if key not in map(lambda x: x.id, model.metabolites):
            raise Exception("cannot find the metabolite %s in the model" % key)

    s = dict([(model.metabolites.get_by_id(key), val)
              for key, val in sparse.iteritems()])
    reaction = Reaction(name)
    reaction.id = id
    reaction.add_metabolites(s)
    reaction.lower_bound = lower_bound
    reaction.upper_bound = upper_bound
    model.add_reactions([reaction])
    return reaction

def add_metabolite_exchange(model, metabolite, lower_bound=-10, upper_bound=0):
    try:
        met = model.metabolites.get_by_id(metabolite + '_c')
    except AttributeError:
        raise KeyError('Model does not have a metabolite with ID: ' + metabolite)
    
    add_metabolite(model, metabolite + '_e', 'E')
    add_reaction(model, metabolite + '_transport', met.name + ' permease',
                 {metabolite + '_c' : -1, metabolite + '_e' : 1}, -1000, 1000)
                 
    add_reaction(model, 'EX_' + metabolite + '_e', met.name + ' exchange',
                 {metabolite + '_e' : -1}, lower_bound, upper_bound)
    
    
model = create_cobra_model_from_sbml_file('data/ecoli_core_model.xml')
mals = model.reactions.get_by_id('MALS')
mals.upper_bound = 0
# the core model has these annoying '_b' metabolites that are used as
# 'ghost' metabolites that balance the exchange reactions. they should
# be ignored in the mass-balance equation and therefore the best way to
# deal with them is to remove them from all the reactions

for m in model.metabolites:
    if m.id.endswith('_b'):
        for r in m.reactions:
            coeff = r.get_coefficient(m)
            r.add_metabolites({m : -coeff})

eat_glucose = model.reactions.get_by_id('EX_glc_e')
eat_glucose.lower_bound = 0

eat_pyruvate = model.reactions.get_by_id('EX_pyr_e')
eat_pyruvate.lower_bound = -10
#eat_pyruvate.upper_bound = 1

accoa = model.metabolites.get_by_id('accoa_c')
pep = model.metabolites.get_by_id('pep_c')
pyr = model.metabolites.get_by_id('pyr_c')

index = [r for r in pyr.reactions if pyr in r.reactants]
df = pd.DataFrame(index = index, columns = ['UP_WT', 'DOWN_WT', 'UP_PG', 
                                            'DOWN_PG', 'UP_Hexa', 'DOWN_Hexa'])
                                            
optimize_minimal_flux(model)
print model.solution.f

for r in index:
    
    if r.id in model.solution.x_dict:

        if pep in r.products:
            df['UP_WT'][r] = model.solution.x_dict[r.id]
        if accoa in r.products:
            df['DOWN_WT'][r] = model.solution.x_dict[r.id]
            
draw_flux('res/pyruvate_biforcation_WT.html')
#
ko = 'G6PDH2r,PGL,PFK,ICL,PGM,PFL'
ki = 'PRK+RBC'
knockout_reactions(model, ko)
knockin_reactions(model, ki)
optimize_minimal_flux(model)
print model.solution.f

for r in index:
    if r.id in model.solution.x_dict:
        if pep in r.products:
            df['UP_PG'][r] = model.solution.x_dict[r.id]
        if accoa in r.products:
            df['DOWN_PG'][r] = model.solution.x_dict[r.id]
            
draw_flux('res/pyruvate_biforcation_PG.html')
#
add_metabolite_exchange(model, 'xu5p_D')
eat_xylose = model.reactions.get_by_id('EX_xu5p_D_e')
eat_xylose.lower_bound = -10

optimize_minimal_flux(model)
print model.solution.f

for r in index:
    if r.id in model.solution.x_dict:
        if pep in r.products:
            df['UP_Hexa'][r] = model.solution.x_dict[r.id]
        if accoa in r.products:
            df['DOWN_Hexa'][r] = model.solution.x_dict[r.id]
            
draw_flux('res/pyruvate_biforcation_Hexa.html')


df.dropna(how='all', inplace=True)


import numpy as np
import matplotlib.pyplot as plt

up_cmap = plt.get_cmap('Blues')
down_cmap = plt.get_cmap('Reds')

UP = df[['UP_WT', 'UP_Hexa', 'UP_PG']].replace(0, np.NaN).dropna(how='all').values.flatten()
DOWN = df[['DOWN_WT', 'DOWN_Hexa', 'DOWN_PG']].replace(0, np.NaN).dropna(how='all').values.flatten()
biomass = 10-(UP+DOWN)


N = 3
ind = np.arange(N)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence
fontsize=15
fig = plt.figure(figsize=(6,4))
ax = plt.axes()
p = plt.bar(ind+width/2, DOWN/UP, width, color='#FFD97F', zorder=3)

#p0 = plt.bar(ind+width, DOWN,   width, color='#FFD97F', label='DOWN')
#p1 = plt.bar(ind+width, UP,   width, color='#66CD6D', bottom=DOWN, label='UP')
#p2 = plt.bar(ind+width, biomass, width, color='#8F5EB2', bottom=DOWN+UP, label='BIOMASS')

plt.ylabel('PDH / PPS', size=fontsize)
plt.xticks(ind+width, ('WT', 'Hexa', 'PG'))
plt.legend(loc=4)
plt.grid(zorder=0)
[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
ax.tick_params(axis='both', which='both', top='off', right='off')

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., height, '%.2f'%height,
                ha='center', va='bottom', size=fontsize)
autolabel(p)               
ax.set_yticks(np.arange(0,13,3))
plt.tight_layout()
plt.savefig('res/pyruvate biforcation ratios.pdf')
