from obiwan import priors
import os

def test_ELG():
    elg= priors.EmptyClass()
    
    d= priors.Data()
    outdir='.'
    #outdir='/home/kaylan/mydata/priors_data'
    #d.fetch(outdir)
    #elg.data= d.load_elg(DR=3)

    #elg.model= KDE_Model('elg',elg.data,outdir)
    #elg.model.kde= elg.model.get_kde()

    p= priors.Plot(outdir)
    #p.elg(elg.data, elg.model.df_wcut,elg.model.df_for_kde,
    #      elg.model.kde)
    assert(p.outdir == os.path.join(outdir,'plots'))
