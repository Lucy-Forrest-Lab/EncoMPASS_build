# Name: run_encompass.py
# Language: python3
# Description: Main wrapper function to generate the EncoMPASS structure set.
#              Can be resumed at different stages.
# Author: Lucy Forrest
# Date: 2025-02-23

from encompass.pipeline.supporting_functions import *
from encompass.pipeline.initialize_repository import *
from encompass.sources.combine_sources import *
from encompass.sources.complete_information import *
from encompass.struct_comparisons.structure_alignment import *


def setup_stage(locations, stage_num, file_label):
  """Create messages, check for and read checkpoints for previous stage
  :param locations:
  :param stage_num:  stage of build, integer
  :param file_label: prefix of checkpoint dictionary filename before stage number
  :return: empty data structure if pickle doesn't exist (assuming first stage), otherwise read checkpoint
  """
  t = time.ctime()
  print(f"\n--- Starting Build Stage {stage_num} ---\nStart Time: {t}")

  if stage_num == 1:
    print(f"WARNING: checkpoint not used for stage {stage_num}, so initializing empty data structure")
    str_data = []
  else:
    pkl_number = stage_num-1
    read_pickle_filepath = locations['FSYSPATH']['cache'] + f"{file_label}{pkl_number}.pkl"
    if not os.path.exists(read_pickle_filepath):
      print(f"FATAL ERROR: No checkpoint {read_pickle_filepath} for previous stage {pkl_number}")
      exit()
    else:
      print("Reading", read_pickle_filepath)
      str_data = read_checkpoint(read_pickle_filepath, locations['SYSFILES']['data_structure_template'])

  return str_data


def finalize_stage(locations, str_data, finished_stage_num, file_label):
    """Write checkpoint for a finished stage and return the next stage number."""
    write_pickle_filepath = locations['FSYSPATH']['cache'] + f"{file_label}{finished_stage_num}.pkl"
    write_checkpoint(str_data, write_pickle_filepath)

    t_end = time.ctime()
    print(
        f"Completed Build Stage {finished_stage_num} and stored data structure in "
        f"{write_pickle_filepath}\nEnd Time: {t_end}"
    )
    following_stage = finished_stage_num + 1
    return following_stage


def run_encompass_pipeline(options=None, locations=None, checkpoint_label: str = "str_data_stage") -> None:
    """
    Run the multi-stage EncoMPASS build pipeline.

    If options and locations are not provided, they are obtained from initialize_repository().
    This is the function that the CLI subcommand will call.
    """
    if options is None or locations is None:
        options, locations = initialize_repository()

    # pass stage with instruction file
    next_stage: int = options['ALL']['start_stage']
    end_stage: int = options['ALL']['end_stage']

    # Execute each stage including and following the provided flag for next_stage
    while next_stage != end_stage:
        if next_stage == 1:
            _ = setup_stage(locations, next_stage, checkpoint_label)  # return will be empty, don't expect pickle for stage 1
            id_list, secd = mainlist_from_OPM(options, locations)
            str_data = scan_OPM(options, locations, id_list)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 2:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            str_data = scan_PDB(options, locations, str_data=str_data)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 3:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            str_data = combine_sources(options, locations, str_data)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 4:
            str_data = setup_stage(locations, next_stage, checkpoint_label)  # read pickle
            # The actual set of structures is all entries that have not been already eliminated; check types
            str_set = {pdbi for pdbi in str_data if 'eliminated' not in str_data[pdbi]['status']}
            print("MUSCLE RUN", str_set, [str_data[x]['status'] for x in str_data])
            str_data = run_MUSCLE(options, locations, str_data, str_set, only_analysis=False)  # False in production
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 5:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            str_data, run_in_PPM = parallel_generate_full_structure(options, locations, str_data)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 6:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            # TODO consider putting code from run_PPM_wrapper in this section to launch run_PPM directly
            str_data = run_PPM_wrapper(options, locations, str_data, only_analysis=False)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 7:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            str_data = structure_sorter(options, locations, str_data)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 8:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            str_data = post_insertion_analysis(options, locations, str_data)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 9:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            str_data, unique_chains = find_redundant_chains(locations, str_data)
            # combined into one stage because of unique chain list
            str_data = assign_uniprot_acc(locations, str_data, unique_chains)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 10:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            str_data = convert_coords_to_enc(locations, str_data, checkpoint_label)
            make_structurewise_table(str_data, locations['SYSFILES']['structurewise_table'])
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 11:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            exelist, pdbich_nte_list = comparison_lists(options, locations, str_data)
            pickle.dump(exelist, open(locations['FSYSPATH']['cache'] + "exelist.pkl", "wb"))
            pickle.dump(pdbich_nte_list, open(locations['FSYSPATH']['cache'] + "pdbich_nte_list.pkl", "wb"))
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 12:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            exelist = pickle.load(open(locations['FSYSPATH']['cache'] + "exelist.pkl", "rb"))
            pdbich_nte_list = pickle.load(open(locations['FSYSPATH']['cache'] + "pdbich_nte_list.pkl", "rb"))
            structure_alignment(options, locations, str_data, exelist, pdbich_nte_list, only_table=False, only_gather=False)
            structure_alignment(options, locations, str_data, exelist, pdbich_nte_list,
                                only_table=False, only_gather=False)
            next_stage = finalize_stage(locations, str_data, next_stage, checkpoint_label)

        elif next_stage == 13:
            str_data = setup_stage(locations, next_stage, checkpoint_label)
            complete_straln(
                exelist_name=locations['SYSFILES']['H_scheduledalns'],
                summary_name=locations['SYSFILES']['summarytable'],
            )
            finalize_stage(locations, str_data, next_stage, checkpoint_label)
            next_stage = -1  # End

        else:
            print(
                "ERROR: there are only 13 stages; change the argument for start_stage, "
                "currently", options['ALL']['start_stage']
            )
            exit()

    final_checkpoint_name = locations['FSYSPATH']['cache'] + "/str_data_completegen.pkl"
    # rename str_data_$tot_stages.pkl to str_data_completegen.pkl
    shutil.copyfile(
        f"{locations['FSYSPATH']['cache']}/{checkpoint_label}10.pkl",
        final_checkpoint_name,
    )

    # Completed all stages
    print(
        "All build stages completed. Saved",
        final_checkpoint_name,
        locations['SYSFILES']['structurewise_table'],
    )


def main(argv=None) -> None:
    """
    CLI-compatible entrypoint.

    The `argv` parameter is accepted for compatibility with console scripts,
    but currently not used, because initialize_repository() handles its own
    configuration / argument parsing.
    """
    # We deliberately ignore argv here; initialize_repository() is the config front-end.
    run_encompass_pipeline()


if __name__ == "__main__":
    main()
