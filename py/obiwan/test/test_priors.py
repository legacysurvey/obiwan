from obiwan import priors

def test_ELG():
    kwargs= dict(savekde=True, loadkde= False,
                 savefig= True,
                 alpha= 0.25, DR=3, outdir='.',
                 rlimit=23.4+1)
    elg= priors.ELG(**kwargs)
    assert(elg.rlimit == 24.4)
