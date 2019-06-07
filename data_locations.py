from starter2 import *
enzo_directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
snapshot_location ='/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe'
if 'machine' in os.environ:
    if os.environ['machine'] in ['shipsterns']:
        enzo_directory = '/home/dcollins/scratch/u05-r4-l4-128'
        snapshot_location = '/home/dcollins/scratch/Paper19/track_index_fix'
