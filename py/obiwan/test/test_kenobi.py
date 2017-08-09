from obiwan import kenobi

def test_output_dir():
    kwargs=dict(decals_sim_dir= 'hello',
                objtype='elg',
                brickname='1238p245',
                rowst=1)
    ans= 'hello/elg/123/1238p245/rs1'
    assert(kenobi.get_savedir(**kwargs) == ans)
