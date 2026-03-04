# Usage: 
# python3 generate_neighborlists.py <summary_table>
# or (uncomment lines):
# python3 -s <instruction_file> -main <encompass_main_folder>

import sys
from encompass.struct_comparisons.analysis_parallel import neighborlists

####### standard init
# from encompass.pipeline.initialize_repository import initialize_repository
# options, locations = initialize_repository()
# summarytable_fn = locations['SYSFILES']['summarytable']
#######


####### ...or custom initialisation for toni
# first create the destination folders, then run
t_locations = {'FSYSPATH' : {}}
t_locations['FSYSPATH']['seqneighs'] = 'test_seq/' # dest folder for seq neighbors
t_locations['FSYSPATH']['strneighs'] = 'test_str/' # dest folder for str neighbors
t_locations['FSYSPATH']['totneighs'] = 'test_tot/' # dest folder for tot neighbors

t_options = {'PARAMETERS' : {}}
t_options['PARAMETERS'][('seqid_thr', 'sequence_identity_threshold')] = 0.85
t_options['PARAMETERS'][('tmscore_thr', 'tmscore_threshold')] = 0.6

summarytable_fn = sys.argv[1]
#######

# command
pdbi_chs = neighborlists(t_options, summarytable_fn, t_locations, no_print=False)

#print all the pdbi_ch entries considered
print("\n".join(sorted(pdbi_chs)))
