import scipy.sparse
from itertools import chain
import os, sys, pickle
from copy import deepcopy
import matplotlib.pyplot as plt
from cobra.flux_analysis.single_deletion import single_reaction_deletion_fba as sd 
from analysis_toolbox import model_summary, plot_multi_PPP
from models import init_wt_model, knockout_reactions, knockin_reactions
from optknock import OptKnock
from draw_flux import DrawFlux
from html_writer import HtmlWriter

def main():

    main_html = HtmlWriter('res/fba.html')
    main_html.write('<h1>Flux Balance Analysis</h1>\n')

    model = init_wt_model('full', {'pyr' : -10}) 
    ko_reactions = 'PGM,ICL,G6PDH2r,PGL,MALS,PFK,'
    ki_reactions = 'PRK+RBC'

    ''' lower to upper leaks'''
    ko = 'GLYCK,' # glyoxylate metabolism - Lior's gene 
    ko += 'DXPS'
    ko_reactions += ko
       
    models = {'WT' : model}

    if ko_reactions:
        for k in models.keys():
            m = deepcopy(models[k])
            knockout_reactions(m, ko_reactions)
            models[k + ' -%s' % ko_reactions] = m

    if ki_reactions:
        for k in models.keys():
            m = deepcopy(models[k])
            knockin_reactions(m, ki_reactions)
            models[k + ' +%s' % ki_reactions] = m

    # Run the optimization for the objective reaction and medium composition
    # set in the file.
    main_html.write('<table border="1">\n')
    main_html.write('<tr><td><b>Model Name</b></td><td><b>Growth Yield</b></td></tr>\n')
    growths = {}
    for name, model in sorted(models.iteritems()):
    
        print "solving %s model\n" % name,
        ok = OptKnock(model)
        ok.prepare_FBA_primal()
        ok.solve()
        growths[name] = ok.get_objective_value()
    
        if growths[name] is None:
            main_html.write('<tr><td>%s</td><td>infeasible</td></tr>\n' % name)
            print ": No solution \n"
        else:
            print ': f = %.3g \n' % growths[name]
            main_html.write('<tr><td>')
            html = main_html.branch(name)
            main_html.write('</td><td>%.3g</td></tr>\n' % growths[name])
            html.write('<h1>Model name: %s</h1>\n' % name)
            html.write('<h2>Growth Yield: %.3g</h2>\n' % growths[name])
            ok.draw_svg(html)
            ok.model_summary(html)
    main_html.write('</table>\n')

if __name__ == "__main__":
    main()
