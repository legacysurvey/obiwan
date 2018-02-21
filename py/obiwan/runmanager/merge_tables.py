
class HeatmapTable(bricklist,nproc,outfn):
    bricks = fits_table(os.path.join(os.environ['LEGACY_SURVEY_DIR'],
                                     'survey-bricks.fits.gz'))
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricklist= np.array_split(bricklist, comm.size)[comm.rank]
    else:
        class MyComm(object):
            def __init__(self):
                self.rank=0
        comm= MyComm()

    one_row_per_brick= []
    one_row_per_brick.append( HeatmapTable().run(bricklist) )
    one_row_per_brick = comm.gather(one_row_per_brick, root=0)
    if comm.rank == 0:
        allRows = merge_tables(one_row_per_brick)
        allRows.writeto(outfn)
        print('Wrote', outfn)


