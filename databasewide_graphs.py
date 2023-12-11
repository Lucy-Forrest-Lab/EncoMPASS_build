# Create plots describing the entire Encompass database
# written by Edoardo Sarti 
# edits by Lucy Forrest
import pandas as pd

from supporting_functions import *
import argparse
import matplotlib
import pickle as pkl


def plot_databasewide_data(table, pdbi_list, plot_path):
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams['font.size'] = 14
    matplotlib.rcParams['font.family'] = 'Arial'

    # thresholds for when stacked data are combined
    MAXNONTMCH = 4
    MAXCH = 5

    # Number of TM chains per complex - old vs stacked
    tmchain_count_per_non_tmchains = define_tmchains_for_stacked_plots(MAXNONTMCH, pdbi_list, table)
    plot_ntm_per_complex_stacked(MAXNONTMCH, table, tmchain_count_per_non_tmchains, plot_path)

    # TM chains vs chains in each complex - scatter plot or heatmap
    plot_ntmch_vs_nch_heatmap(table, plot_path)

    # TM regions per TM chain - old vs stacked
    tmregion_count_per_tmchains = define_ntmregions_for_stacked_plots(MAXCH, table)
    plot_tmreg_per_tmch_stacked(MAXCH, MAXNONTMCH, tmregion_count_per_tmchains, table, plot_path)

    # TM regions per complex - hist and heatmap vs TM chains
    num_tm_chains_in_complex, num_tm_regions_in_complex = define_num_tmreg_and_tmch_per_complex(pdbi_list, table)
    plot_tmreg_per_complex(num_tm_regions_in_complex, plot_path)
    plot_tmreg_per_complex_zoom(num_tm_regions_in_complex, plot_path)
    plot_tmres_per_complex_vs_tmchs_heatmap(num_tm_chains_in_complex, num_tm_regions_in_complex, plot_path)
    plot_tmregs_per_complex_vs_tmchs_heatmap_zoom(num_tm_chains_in_complex, num_tm_regions_in_complex, plot_path)

    # TM regions vs number of TM chains per chain - scatter and heatmap
    plot_tmreg_vs_tmch_heatmap(table, plot_path)

    # TM region coverage and SS plots
    plot_tmcoverage_histogram(table, plot_path)
    plot_tmcoverage_vs_chlength(table, plot_path)
    plot_sscoverage_vs_tmspans(table, plot_path)


def read_and_parse_chain_table(chain_structural_properties_file, output_path, exclusions):
    if os.path.isfile(output_path + "chain_structural_properties_filtered.pkl"):
        print("Reading checkpoint file from previous run", output_path + "chain_structural_properties_filtered.pkl")
        table = pkl.load(open(output_path + "chain_structural_properties_filtered.pkl", "rb"))
    else:
        table = read_and_prune_chain_list(chain_structural_properties_file, output_path, exclusions)
        print("Parsed", chain_structural_properties_file)
    return table


def read_and_parse_pairs_table(structure_alignment_pairs_file, output_path):
    if os.path.isfile(output_path + "pairwise_list_filtered.pkl"):
        print("Reading checkpoint file from previous run", output_path + "pairwise_list_filtered.pkl")
        pair_data = pkl.load(open(output_path + "pairwise_list_filtered.pkl", "rb"))
    else:
        pair_data = read_and_prune_pair_list(structure_alignment_pairs_file, output_path)
        print("Parsed", structure_alignment_pairs_file)

    return pair_data


def parse_nonredundant_list(non_redundant_entry_file, table, exclude_list):
    non_redundant_entry_list = load_non_redundant_entries(non_redundant_entry_file)

    nr_list = np.unique(non_redundant_entry_list['Representatives'])
    print("Using non-redundant list of entries from", non_redundant_entry_file, "containing", len(nr_list), "entries")

    for exclude_i in exclude_list:
        exclude_index = np.argwhere(nr_list == exclude_i)
        nr_list = np.delete(nr_list, exclude_index)

    print("Returning pdbi_list containing", len(nr_list), "entries, after excluding", len(exclude_list), "entries")
    table = reduce_chain_table_redundancy(nr_list, table)
    print("\nFiltered dataset now contains nr_list:", len(nr_list), "entries, in:", table)
    return nr_list, table


def load_non_redundant_entries(redundancy_list):
    print("Reading", redundancy_list)
    column_names = ["Representatives", "Neighbors"]
    df = pd.read_csv(redundancy_list, header=0, sep='\t', names=column_names)
    non_redundant_entries = df.drop("Neighbors", axis='columns')
    print("Number of non-redundant entries in ", redundancy_list, "is", len(non_redundant_entries))
    return non_redundant_entries


def reduce_chain_table_redundancy(nr_list, table):
    print("Original table contains this many entries:", len(table))
    nr_table = pd.DataFrame()
    for rep in nr_list:
        nr_table = nr_table.append(table[table['PDB_CODE'] == rep])
    table = nr_table
    print("Table filtered with nr_list now contains non-redundant entries:", len(table))
    return table


def plot_tmreg_per_complex(num_tm_regions_in_complex, plot_path):
    plot_name = "tmregs_in_complex"
    plt.clf()

    matplotlib.rcParams['font.size'] = 30
    fig, ax = plt.subplots(figsize=(24, 8))
    bins = np.arange(0.5, max(num_tm_regions_in_complex) + 1.5, 1)

    # - ticks (more ticks than labels)
    ax.set_xticks([i for i in range(0, 251, 25)])
    print(plot_name, "ticks are:", len([i for i in range(0, 250, 10)]), len([str(x) for x in range(0, 251, 50)]))
    ax.set_xticklabels([str(x) for x in range(0, 251, 25)])

    ax.hist(num_tm_regions_in_complex, bins=bins)
    ax.set_xlabel('Total num. of transmembrane regions per complex')
    ax.set_ylabel('Num. of complexes')

    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(num_tm_regions_in_complex, bins, datafilename)


def plot_tmreg_per_complex_zoom(num_tm_regions_in_complex, plot_path):
    plot_name = "tmregs_in_complex_zoom"
    plt.clf()

    matplotlib.rcParams['font.size'] = 30
    fig, ax = plt.subplots(figsize=(24, 12))
    #    fig, ax = plt.subplots()
    bins = np.arange(0.5, max(num_tm_regions_in_complex) + 1.5, 1)

    ax.set_xlim([99, 210])
    ax.set_ylim([0, 5])
    ax.hist(num_tm_regions_in_complex, bins=bins)
    ax.set_xlabel('Total num. of transmembrane regions per complex')
    ax.set_ylabel('Num. of complexes')

    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(num_tm_regions_in_complex, bins, datafilename)


def plot_tmregs_per_complex_vs_tmchs_heatmap_zoom(num_tm_chains_in_complex, num_tm_regions_in_complex, plot_path):
    plot_name = "tmregs_in_complex_vs_tmchs_heatmap_zoom"
    maxx, minx, maxy, miny = max(num_tm_regions_in_complex), min(num_tm_regions_in_complex), \
                             max(num_tm_chains_in_complex), min(num_tm_chains_in_complex)
    maxx = 30  # zoom in on region up to 30 regions or chains
    maxy = 30
    y2x_ratio = (maxy - miny) / (maxx - minx)
    rect_histx, rect_histy, rect_scatter = define_heatmap_xyhist_dimensions()

    plt.clf()
    fig = plt.figure(figsize=(8, 8 * y2x_ratio))
    ax, axHistx, axHisty = create_heatmap_xyhist_axes(rect_histx, rect_histy, rect_scatter)
    print(plot_name, ", NTMS_in_COMPLEX", minx, maxx, "TMCH", miny, maxy)
    mtx = np.ones((maxx - minx + 1, maxy - miny + 1))
    zoom_ntms_in_complex, zoom_tmch_in_complex = [], []
    for i in range(len(num_tm_regions_in_complex)):
        if num_tm_regions_in_complex[i] <= maxx and num_tm_chains_in_complex[i] <= maxy:
            mtx[num_tm_regions_in_complex[i] - minx, num_tm_chains_in_complex[i] - miny] += 1
            if num_tm_regions_in_complex[i] <= maxx and num_tm_chains_in_complex[i] <= maxy:
                zoom_ntms_in_complex.append(num_tm_regions_in_complex[i])
                zoom_tmch_in_complex.append(num_tm_chains_in_complex[i])
    ax.imshow(np.log(mtx.T), cmap='Blues', origin='lower', interpolation='none')
    ax.set_xlabel('Num. of transmembrane regions per complex')
    ax.set_ylabel('Num. of TM chains per complex')
    tick_gap = 5
    set_gapped_ticks(ax, tick_gap, mtx)
    create_xy_histograms(ax, axHistx, axHisty, maxx, maxy, minx, miny, zoom_ntms_in_complex, zoom_tmch_in_complex)
    outfilename = plot_path + plot_name
    write_plots(outfilename)


def plot_tmres_per_complex_vs_tmchs_heatmap(num_tm_chains_in_complex, num_tm_regions_in_complex, plot_path):
    plot_name = "tmregs_in_complex_vs_tmchs_heatmap"
    # - main sizes
    maxx, minx, maxy, miny = max(num_tm_regions_in_complex), min(num_tm_regions_in_complex), \
                             max(num_tm_chains_in_complex), min(num_tm_chains_in_complex)
    y2x_ratio = (maxy - miny) / (maxx - minx)
    matplotlib.rcParams['font.size'] = 14
    rect_histx, rect_histy, rect_scatter = define_plot_size_y2x(y2x_ratio)
    plt.clf()
    fig = plt.figure(figsize=(30, 30 * y2x_ratio))
    ax, axHistx, axHisty = create_heatmap_xyhist_axes(rect_histx, rect_histy, rect_scatter)

    # - main plot
    maxx, minx, maxy, miny = max(num_tm_regions_in_complex), min(num_tm_regions_in_complex), max(
        num_tm_chains_in_complex), min(num_tm_chains_in_complex)
    row_labels = [i for i in range(minx, maxx + 1)]
    col_labels = [i for i in range(miny, maxy + 1)]
    print(plot_name, ": tm_regions_in_complex", minx, maxx, "tm_chains_in_complex", miny, maxy)
    mtx = np.ones((maxx - minx + 1, maxy - miny + 1))
    for i in range(len(num_tm_regions_in_complex)):
        mtx[num_tm_regions_in_complex[i] - minx, num_tm_chains_in_complex[i] - miny] += 1
    ax.imshow(np.log(mtx.T), cmap='Blues', origin='lower', interpolation='none')
    ax.set_xlabel('Num. of transmembrane regions per complex')
    ax.set_ylabel('Num. of TM chains per complex')

    tick_gap = 25
    set_gapped_ticks(ax, tick_gap, mtx)

    create_xy_histograms(ax, axHistx, axHisty, maxx, maxy, minx, miny,
                         num_tm_regions_in_complex, num_tm_chains_in_complex)

    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_data_matrix(mtx - 1, row_labels, col_labels, datafilename)


def plot_sscoverage_vs_tmspans(table, plot_path):
    plot_name = "sscovintmregs_hist"
    plt.clf()
    fig, ax = plt.subplots()
    i, j = 'SS_COVERAGE_IN_TMs', 'TM_COVERAGE'
    bins = np.arange(-0.01, 1.01, 0.05)
    ax.hist(table[i] / table[j], bins=bins)  # data here
    ax.set_xlabel('Fraction of secondary structure in the membrane per chain')
    ax.set_ylabel('Num. of membrane-spanning chains')
    # ax.set_title('Secondary structure coverage in the membrane vs TM coverage')
    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(table[i] / table[j], bins, datafilename)


def plot_num_struct_align_per_tmregion(MAXDIFF, ncomp_stacked, pair_data, plot_path):
    plot_name = "naln_stacked"
    i = 'NTM1'
    plt.clf()
    fig, ax = plt.subplots()
    bins = np.arange(0.5, max(pair_data[i]) + 1.5, 1)
    colors = define_9color_scale(MAXDIFF)
    legend = define_9color_legend(MAXDIFF)

    ax.hist([ncomp_stacked[i] for i in range(-MAXDIFF, MAXDIFF + 1)], bins=bins, stacked=True, label=legend,
            color=colors, edgecolor='black', linewidth=0.5)
    add_color_legend(ax)
    ax.set_xlabel('Num. of transmembrane regions per chain')
    ax.set_ylabel('Num. of structure-based alignments')
    # ax.set_title('Number of structure alignments')
    outfilename = plot_path + plot_name
    plt.savefig(outfilename + '.png', bbox_inches="tight")
    plt.savefig(outfilename + '.pdf', bbox_inches="tight")
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(ncomp_stacked, bins, datafilename, stacked=True)


def plot_tmcoverage_vs_chlength(table, plot_path):
    plot_name = "tmcov_vs_chlen"
    plt.clf()
    fig, ax = plt.subplots()
    i, j = 'CHAIN_LENGTH', 'TM_COVERAGE'
    ax.scatter(table[i].to_list(), (table[j] / table[i]).to_list())

    ax.set_xlabel('Chain length [residues]')
    ax.set_ylabel('Fraction of chain in the membrane')
    outfilename = plot_path + plot_name
    write_plots(outfilename)

    # heatmap
    plot_name = "tmcov_vs_chlen_heatmap"
    chain_len, tmcov = table['CHAIN_LENGTH'].to_numpy(), \
                       table['TM_COVERAGE'].to_numpy() / table['CHAIN_LENGTH'].to_numpy()
    maxx, minx, maxy, miny = max(chain_len), min(chain_len), max(tmcov), min(tmcov)
    bw_row = 10
    bw_col = 0.1
    print("Max Min x, y", maxx, minx, maxy, miny)
    y2x_ratio = ((maxy - miny) * bw_row) / ((maxx - minx) * bw_col)
    print("RATIO", y2x_ratio, 4 / y2x_ratio, 4)
    rect_histx, rect_histy, rect_scatter = define_plot_size_y2x(y2x_ratio)

    plt.clf()
    fig = plt.figure(figsize=(8 / y2x_ratio, 8))
    ax, axHistx, axHisty = create_heatmap_xyhist_axes(rect_histx, rect_histy, rect_scatter)

    # - main plot
    bw_row = 10
    row_labels = [i for i in range(minx, maxx + bw_row, bw_row)]
    bw_col = 0.1
    col_labels = [i for i in np.arange(miny, maxy + bw_col, bw_col)]
    print(plot_name, " chain length:", minx, maxx, "coverage:", miny, maxy)
    mtx = np.ones((int((maxx - minx + bw_row) / bw_row), int((maxy - miny + bw_col) / bw_col)))
    for i in range(len(chain_len)):
        mtx[int((chain_len[i] - minx) / bw_row), int((tmcov[i] - miny) / bw_col)] += 1
    ax.imshow(np.log(mtx.T), cmap='Blues', origin='lower', interpolation='none')
    ax.set_xlabel('Chain length [residues]')
    ax.set_ylabel('Fraction of chain in the membrane')
    create_xy_histograms(ax, axHistx, axHisty, maxx, maxy, minx, miny, chain_len - minx, tmcov - miny)

    # - x-axis histogram - modifications to default
    bins = np.arange(minx - 0.5 * bw_row, maxx + 1.5 * bw_row, bw_row)
    axHistx.set_xlim((minx, maxx))

    # - y-axis histogram - modifications to default
    bins = np.arange(miny - 0.5 * bw_col, maxy + 1.5 * bw_col, bw_col)
    axHisty.set_ylim((miny, maxy))

    # fig.suptitle('Transmembrane coverage vs chain length')
    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_data_matrix(mtx - 1, row_labels, col_labels, datafilename)


def plot_tmcoverage_histogram(table, plot_path):
    plot_name = "tmcov_hist"
    plt.clf()
    fig, ax = plt.subplots()
    i, j = 'TM_COVERAGE', 'CHAIN_LENGTH'
    bins = np.arange(-0.05, 1.05, 0.05)
    table['PDBCH'] = list(zip(table['PDB_CODE'], table['CHAIN_ID']))
    ax.hist(table[i] / table[j], bins=bins)
    ax.set_xlabel('Fraction of chain in the membrane')
    ax.set_ylabel('Num. of membrane-spanning chains')
    # ax.set_title('Transmembrane coverage of TM chains')
    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(table[i] / table[j], bins, datafilename)


def plot_tmreg_vs_tmch_heatmap(table, plot_path):
    plot_name = "tmregs_vs_tmchains_heatmap"
    tmch, ntm = table['N_TM_CHAINS'].to_numpy(), table['N_TMs'].to_numpy()
    maxx, minx, maxy, miny = max(tmch), min(tmch), max(ntm), min(ntm)
    y2x_ratio = (maxy - miny) / (maxx - minx)

    rect_histx, rect_histy, rect_scatter = define_heatmap_xyhist_dimensions()
    plt.clf()
    fig = plt.figure(figsize=(8, 8 * y2x_ratio))
    ax, axHistx, axHisty = create_heatmap_xyhist_axes(rect_histx, rect_histy, rect_scatter)

    #    tmch, ntm = table['N_TM_CHAINS'].to_numpy(), table['N_TMs'].to_numpy()
    #    maxx, minx, maxy, miny = max(tmch), min(tmch), max(ntm), min(ntm)
    row_labels = [i for i in range(minx, maxx + 1)]
    col_labels = [i for i in range(miny, maxy + 1)]
    print("tmregs_vs_tmchains_heatmap.png Num_TM_chains", minx, maxx, "Num_TM_regions", miny, maxy)
    miny = 1  # We don't track NTM==0
    mtx = np.ones((maxx - minx + 1, maxy - miny + 1))
    for i in range(len(tmch)):
        if ntm[i] > 0:
            mtx[tmch[i] - minx, ntm[i] - miny] += 1
    ax.imshow(np.log(mtx.T), cmap='Blues', origin='lower', interpolation='none')
    ax.set_xlabel('Num. membrane-spanning chains per complex')
    ax.set_ylabel('Num. transmembrane regions per chain')

    tick_gap = 10
    set_gapped_ticks(ax, tick_gap, mtx)

    create_xy_histograms(ax, axHistx, axHisty, maxx, maxy, minx, miny, tmch, ntm)
    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_data_matrix(mtx - 1, row_labels, col_labels, datafilename)


def plot_tmreg_vs_tmch_scatter(table, plot_path):
    plot_name = "tmregs_vs_tmchains"
    plt.clf()
    fig, ax = plt.subplots()
    num_tm_chains, num_tm_regions = table['N_TM_CHAINS'].to_numpy(), table['N_TMs'].to_numpy()
    ax.scatter(table[num_tm_chains], table[num_tm_regions])
    ax.set_xlabel('#TM chains in a complex')
    ax.set_ylabel('#TM regions in a TM chain')
    ax.set_title('Number of TM regions in a chain vs. number of chains in a complex')
    outfilename = plot_path + plot_name
    write_plots(outfilename)


def plot_ntmch_vs_nch_scatter(table, plot_path):
    plot_name = "tmchains_vs_chains"
    plt.clf()
    fig, ax = plt.subplots()
    num_chains, num_tmchains = table['N_CHAINS'].to_numpy(), table['N_TM_CHAINS'].to_numpy()
    ax.scatter(num_chains, num_tmchains)
    ax.set_xlabel('Num. of chains per complex')
    ax.set_ylabel('Num. of membrane-spanning chains per complex')
    outfilename = plot_path + plot_name
    write_plots(outfilename)


def plot_tmreg_per_tmch_stacked(MAXCH, MAXNONTMCH, ntms_stacked, table, plot_path):
    plot_name = "tmregs_hist_stacked"
    plt.clf()
    fig, ax = plt.subplots()
    i = 'N_TMs'
    bins = np.arange(0.5, max(table[i]) + 0.5, 1)
    colors = define_4color_scale(MAXNONTMCH)  # previously set to MAXCH - 1, which here is same as MAXNONTMCH
    legend = define_4color_stack_legend(MAXCH)

    # - main plot
    ax.hist([ntms_stacked[i] for i in range(1, MAXCH + 1)], bins=bins, stacked=True, label=legend, color=colors,
            edgecolor='black', linewidth=0.5)
    add_color_legend(ax)

    ax.set_xlabel('Num. of transmembrane regions per chain')
    ax.set_ylabel('Num. of membrane-spanning chains')
    outfilename = plot_path + plot_name
    plt.savefig(outfilename + '.png', bbox_inches="tight")
    plt.savefig(outfilename + '.pdf', bbox_inches="tight")
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(ntms_stacked, bins, datafilename, stacked=True)


def plot_tmreg_per_tmch_hist(table, plot_path):
    plot_name = "tmregs_hist"
    plt.clf()
    fig, ax = plt.subplots()
    i = 'N_TMs'
    bins = np.arange(0.5, max(table[i]) + 0.5, 1)
    ax.hist(table[i], bins=bins)
    ax.set_xlabel('#TM regions in a chain')
    ax.set_ylabel('#TM chains in EncoMPASS')
    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(table[i].to_list(), bins, datafilename)


def plot_ntmch_vs_nch_heatmap(table, plot_path):
    plot_name = "tmchains_vs_chains_heatmap"

    # - main plot: heatmap
    plt.clf()

    num_chains_per_pdbi = table['N_CHAINS'].to_numpy()
    tmchains_per_pdbi = table['N_TM_CHAINS'].to_numpy()
    max_chains = max(num_chains_per_pdbi)
    max_tmchains = max(tmchains_per_pdbi)
    print("max_chains", max_chains, "max_tmchains", max_tmchains)

    maxx, minx, maxy, miny = max(max_chains, max_tmchains), min(tmchains_per_pdbi), \
                             max(max_chains, max_tmchains), min(num_chains_per_pdbi)
    y2x_ratio = (maxy - miny) / (maxx - minx)
    fig = plt.figure(figsize=(8, 8 * y2x_ratio))
    rect_histx, rect_histy, rect_scatter = define_heatmap_xyhist_dimensions()
    ax, axHistx, axHisty = create_heatmap_xyhist_axes(rect_histx, rect_histy, rect_scatter)
    row_labels = [i for i in range(minx, maxx + 1)]  # TODO ask why i for i
    col_labels = [i for i in range(miny, maxy + 1)]
    print(plot_name, " axis range to use: NCH", minx, maxx, "TMCH", miny, maxy)

    mtx = np.ones((maxx - minx + 1, maxy - miny + 1))
    for i in range(len(num_chains_per_pdbi)):
        mtx[tmchains_per_pdbi[i] - minx, num_chains_per_pdbi[i] - miny] += 1
    ax.imshow(np.log(mtx.T), cmap='Blues', origin='lower', interpolation='none')
    ax.set_ylabel('Num. of chains per complex')
    ax.set_xlabel('Num. of membrane-spanning chains')

    tick_gap = 10
    set_gapped_ticks(ax, tick_gap, mtx)
    create_xy_histograms(ax, axHistx, axHisty, maxx, maxy, minx, miny, tmchains_per_pdbi, num_chains_per_pdbi)

    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_data_matrix(mtx - 1, row_labels, col_labels, datafilename)


def plot_ntm_per_complex(table, plot_path):
    plot_name = "tmchains_hist"
    plt.clf()
    i = 'N_TMs'
    tmchains = table[i]
    fig, ax = plt.subplots()
    bins = np.arange(0.5, max(table['N_TM_CHAINS']) + 1.5, 1)
    ax.hist(tmchains, bins=bins)  # TODO check this works
    ax.set_xlabel('#TM chains')
    ax.set_ylabel('#Complexes in EncoMPASS')
    ax.set_title('Number of TM chains per complex')
    outfilename = plot_path + plot_name
    write_plots(outfilename)
    datafilename = outfilename + '.txt'
    save_discrete_histogram_data(tmchains, bins, datafilename)


def plot_ntm_per_complex_stacked(max_non_tm_chains, table, tmchains_stacked, plot_path):
    plot_name = "tmchains_stacked_hist"
    plt.clf()
    fig, ax = plt.subplots()
    bins = np.arange(0.5, max(table['N_TM_CHAINS']) + 1.5, 1)
    colors = define_4color_scale(max_non_tm_chains)
    legend = define_4color_stack_legend(max_non_tm_chains)

    ax.hist([tmchains_stacked[i] for i in range(max_non_tm_chains + 1)],
            bins=bins, stacked=True, label=legend, color=colors, edgecolor='black', linewidth=0.5)
    add_color_legend(ax)
    ax.set_xlabel('Num. of membrane-spanning chains')
    ax.set_ylabel('Num. of complexes')

    outfilename = plot_path + plot_name
    plt.savefig(outfilename + '.png', bbox_inches="tight")  # TODO - figure out how to pass to write_plots
    plt.savefig(outfilename + '.pdf', bbox_inches="tight")
    datafilename = outfilename + ".txt"
    save_discrete_histogram_data(tmchains_stacked, bins, datafilename, stacked=True)


def define_tmchains_per_pdbi(pdbi_list, table):  # TODO - check obsolete
    tmchains_per_pdbi = []  # TM chains per pdb code
    for pdbi in pdbi_list:
        num_tm_chains = table[table['PDB_CODE'] == pdbi]['N_TM_CHAINS'].iloc[0]
        tmchains_per_pdbi.append(num_tm_chains)
    return tmchains_per_pdbi


def define_num_comparisons_for_stacked_plots(max_diff_in_tm_regions, pair_data):
    ncomp_stacked = {}  # number of comparisons per TM region, divided by difference in TM regions
    for i in range(-max_diff_in_tm_regions, max_diff_in_tm_regions + 1):
        h1 = pair_data[pair_data['NTM1'] == pair_data['NTM2'] - i]['NTM1']
        ncomp_stacked[i] = h1
    return ncomp_stacked


def define_num_tmreg_and_tmch_per_complex(pdbi_list, table):
    num_tm_chains_in_complex = []
    num_tm_regions_in_complex = []
    for pdbi in pdbi_list:
        sumtms = np.sum(table[table['PDB_CODE'] == pdbi]['N_TMs'])  # Sum of TMs in every row with a certain PDB_CODE
        num_tm_regions_in_complex.append(sumtms)
        num_tm_chains_in_complex.append(table[table['PDB_CODE'] == pdbi]['N_TM_CHAINS'].iloc[0])
    #    print("Total TM regions in complex", pdbi, sumtms)

    num_tm_chains_in_complex = np.array(num_tm_chains_in_complex)
    num_tm_regions_in_complex = np.array(num_tm_regions_in_complex)
    return num_tm_chains_in_complex, num_tm_regions_in_complex


def define_ntmregions_for_stacked_plots(max_num_chains, table):
    tmregion_count_per_num_tmchains = {}
    for i in range(1, max_num_chains):
        tmregion_count_per_num_tmchains[i] = table[table['N_CHAINS'] == i]['N_TMs'].to_list()
    tmregion_count_per_num_tmchains[max_num_chains] = table[table['N_CHAINS'] >= max_num_chains]['N_TMs'].to_list()
    #    print("tmregion_count_per_num_tmchains=", tmregion_count_per_num_tmchains)

    return tmregion_count_per_num_tmchains


def define_tmchains_for_stacked_plots(MAXNONTMCH, pdbi_list, table):
    tmchain_count_per_non_tmchains = {}  # ? dictionary of TM chains for specific number of non-TM chains
    num_chains_per_pdbi = []
    for pdbi in pdbi_list:

        num_chs_in_current_pdbi = table[table['PDB_CODE'] == pdbi]['N_CHAINS'].iloc[0]
        num_chains_per_pdbi.append(num_chs_in_current_pdbi)
        num_tm_chains = table[table['PDB_CODE'] == pdbi]['N_TM_CHAINS'].iloc[0]

        num_nontm_chains = min(num_chs_in_current_pdbi - num_tm_chains, MAXNONTMCH)
        if num_nontm_chains not in tmchain_count_per_non_tmchains:  # initializes non-tm rows of structure
            tmchain_count_per_non_tmchains[num_nontm_chains] = []
        tmchain_count_per_non_tmchains[num_nontm_chains].append(num_tm_chains)

    return tmchain_count_per_non_tmchains


def save_discrete_histogram_data(mainlist, bins, output_data_file, stacked=False):
    if not stacked:
        clist = {'': mainlist}
        w = ''
    else:
        clist = mainlist
        w = 'STACK'
        # print("STACK THIS", clist)
    with open(output_data_file, 'w') as f:
        f.write('#BIN\tNPOINTS\t{0}\n'.format(w))
        for ist, l in clist.items():
            dig = np.digitize(l, bins=bins)
            for i, b in enumerate(bins):
                f.write('{0}\t{1}\t{2}\n'.format(bins[i] - 0.5, len(dig[dig == i]), ist))


def save_data_matrix(mtx, row_labels, col_labels, output_matrix_file):
    with open(output_matrix_file, 'w') as f:
        for i in range(mtx.shape[0]):
            for j in range(mtx.shape[1]):
                f.write('{0}\t{1}\t{2}\n'.format(row_labels[i], col_labels[j], mtx[i, j]))


def read_and_prune_chain_list(chain_structural_properties_file, output_folder, exclusions):

    # NOTE: assume rows of pandas-generated table are NOT numbered
    structural_properties = ['PDB_CODE', 'CHAIN_ID', 'UNIPROT_ID', 'N_CHAINS', 'N_TM_CHAINS', 'CHAIN_LENGTH',
                             'N_TEs', 'N_TMs', 'TM_COVERAGE', 'EXTERNAL_NTERM_LENGTH', 'EXTERNAL_MIDPARTS_LENGTH',
                             'EXTERNAL_CTERM_LENGTH', 'SS_COVERAGE', 'SS_COVERAGE_IN_TMs']
    print("Reading and pruning membrane-only chains from", chain_structural_properties_file)

    df = pd.read_csv(chain_structural_properties_file, header=0, sep='\t', names=structural_properties)
    df.info()
    table = df.drop(['UNIPROT_ID'], axis=1)
    table = table[table['N_TM_CHAINS'] > 0]
    table = table[table['N_CHAINS'] > 0]

    if len(exclusions) > 0:
        print("Need to exclude: \n", exclusions)
        print("Table before excluding: ")
        print(table)
        for drop_i in exclusions:
            # print("dropping entry: ", drop_i)
            table = table.loc[~(table['PDB_CODE'] == drop_i)]
        print("Table after excluding: ")
        print(table)

    pkl.dump(table, open(output_folder + "chain_structural_properties_filtered.pkl", "wb"))
    print("Created pruned checkpoint file", output_folder + "chain_structural_properties_filtered.pkl")

    return table


def read_and_prune_pair_list(structure_alignment_pairs_file, output_folder):
    # NOTE: assume rows of pandas-generated table are NOT numbered
    structural_properties = ['PDB_CODE1', 'CHAIN_ID1', 'PDB_CODE2', 'CHAIN_ID2', 'IS_BETA1', 'IS_BETA2',
                             'NTM1', 'NTM2', 'DUPL1', 'DUPL2', 'TMSET1', 'TMSET2']
    print("Reading and pruning IS_BETA2, DUPL* and TMSET* columns from", structure_alignment_pairs_file)
    df = pd.read_csv(structure_alignment_pairs_file, header=1, sep='\t', names=structural_properties)

    table = df.drop(['IS_BETA2', 'DUPL1', 'DUPL2', 'TMSET1', 'TMSET2'], axis=1)
    pkl.dump(table, open(output_folder + "pairwise_list_filtered.pkl", "wb"))
    print("Created pruned checkpoint file", output_folder + "pairwise_list_filtered.pkl")
    return table


def write_plots(root_filename):
    plt.savefig(root_filename + '.png')
    plt.savefig(root_filename + '.pdf')


def add_color_legend(ax):
    handles, labels = ax.get_legend_handles_labels()
    lgd = plt.legend(handles[::-1], labels[::-1], loc=1, ncol=1)
    ax.add_artist(lgd)


def create_heatmap_xyhist_axes(rect_histx, rect_histy, rect_scatter):
    ax = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    return ax, axHistx, axHisty


def set_gapped_ticks(ax, tick_gap, mtx):
    xM, yM = mtx.shape
    ax.set_xticks([i for i in range(0, xM, tick_gap)])
    ax.set_xticklabels([str(x) for x in range(1, xM + 1, tick_gap)])
    ax.set_yticks([i for i in range(0, yM, tick_gap)])
    ax.set_yticklabels([str(x) for x in range(1, yM + 1, tick_gap)])


def define_plot_size_y2x(y2x_ratio):
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = bottom + height + 0.02
    left_h = left + width + 0.02 * y2x_ratio
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2 * y2x_ratio, height]
    return rect_histx, rect_histy, rect_scatter


def create_xy_histograms(ax, axHistx, axHisty, maxx, maxy, minx, miny, xdata, ydata):
    # - x-axis histogram
    plt.locator_params(nbins=4)
    bins = np.arange(minx - 0.5, maxx + 1.5, 1)
    axHistx.hist(xdata - minx + 1, bins=bins, color='blue')
    axHistx.set_xlim((ax.get_xlim()[0] + 1, ax.get_xlim()[1] + 1))
    axHistx.set_yscale('log')
    axHistx.set_xticks([])
    axHistx.set_ylim(bottom=1)
    axHistx.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0, subs=(1.0,), numticks=100))
    axHistx.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=100))
    axHistx.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    # - y-axis histogram
    plt.locator_params(nbins=4)
    bins = np.arange(miny - 0.5, maxy + 1.5, 1)
    axHisty.hist(ydata - miny + 1, bins=bins, orientation='horizontal', color='blue')
    axHisty.set_ylim((ax.get_ylim()[0] + 1, ax.get_ylim()[1] + 1))
    axHisty.set_xscale('log')
    axHisty.set_yticks([])
    axHisty.set_xlim(left=1)
    axHisty.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0, subs=(1.0,), numticks=100))
    axHisty.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=100))
    axHisty.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def define_9color_legend(max_difference):
    legend = []
    for ic in range(-max_difference, max_difference + 1):
        legend.append('{0:+d} TM regions'.format(ic))
    legend[0] = '<{0:+d} TM regions'.format(1 - max_difference)
    legend[-1] = '>{0:+d} TM regions'.format(max_difference - 1)
    return legend


def define_4color_scale(max_num_colors):
    colhues_half1 = int(max_num_colors / 2)
    colhues_half2 = max_num_colors - colhues_half1 + 1
    colors = []
    for i in reversed(range(1, colhues_half1 + 1)):
        colors.append(
            (colhues_half1 / (colhues_half1 + 1), colhues_half1 / (colhues_half1 + 1), i / (colhues_half1 + 1)))
    for i in reversed(range(colhues_half2)):
        colors.append((colhues_half2 / (colhues_half2 + 1), i / (colhues_half2 + 1), 1 / (colhues_half2 + 1)))
    #    print("COLORS", colors)
    return colors


def define_9color_scale(max_difference_in_TM_region_count):
    colhues_half1 = int(max_difference_in_TM_region_count / 2)
    colhues_half2 = max_difference_in_TM_region_count - colhues_half1 + 1
    colors = []
    for i in range(colhues_half1 + 1):
        colors.append((1 / (colhues_half1 + 1), i / (colhues_half1 + 1), colhues_half1 / (colhues_half1 + 1)))
    for i in range(2, colhues_half2 + 1):
        colors.append(
            (i / (colhues_half2 + 1), colhues_half2 / (colhues_half2 + 1), colhues_half2 / (colhues_half2 + 1)))
    colhues_half1 = int(max_difference_in_TM_region_count / 2) + 1
    colhues_half2 = max_difference_in_TM_region_count - colhues_half1 + 1
    for i in reversed(range(1, colhues_half1)):
        colors.append(
            (colhues_half1 / (colhues_half1 + 1), colhues_half1 / (colhues_half1 + 1), i / (colhues_half1 + 1)))
    for i in reversed(range(colhues_half2)):
        colors.append((colhues_half2 / (colhues_half2 + 1), i / (colhues_half2 + 1), 1 / (colhues_half2 + 1)))
    #    print("COLORS", colors)
    return colors


def define_4color_stack_legend(max_num_subunits_to_analyze):
    legend = ['all membrane subunits']
    for ic in range(1, max_num_subunits_to_analyze):
        legend.append('{0:d} non-membrane subunits'.format(ic))
    legend.append('>{0:d} non-membrane subunits'.format(max_num_subunits_to_analyze - 1))
    return legend


def define_heatmap_xyhist_dimensions():
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    rect_scatter = [left, bottom, width, height]  # central heatmap plot
    rect_histx = [left, bottom_h, width, 0.2]  # x-axis histogram
    rect_histy = [left_h, bottom, 0.2, height]  # y-axis histogram
    return rect_histx, rect_histy, rect_scatter


def read_exclusions(delete_list):
    exclude_pdbs = set({})
    if delete_list:
        with open(delete_list, 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    exclude_pdbs.add(line[0:4])
                    #print("Excluding:", line[0:4])
    return exclude_pdbs


def parse_args():
    parser = argparse.ArgumentParser(description="\nCreate plots describing entire Encompass database\n")
    parser.add_argument("-c", "--chain_file", nargs="?", required=True, type=str,
                        help="Table containing properties of each chain, e.g. structurewise_table.txt")
    parser.add_argument("-p", "--pair_file", nargs="?", required=True, type=str,
                        help="Pair-wise list reflecting structural alignments, e.g. scheduled_alignments.dat")
    parser.add_argument("-o", "--output_folder", nargs="?", required=True, type=str,
                        help="Folder in which plots will be generated")

    # Optional
    parser.add_argument("-nr", "--nonredundant_chain_list", type=str,
                        help="optional file containing list of non-redundant pdb codes")
    parser.add_argument("-nr_label", "--nr_threshold_label", type=str,
                        help="label for plot files, e.g. nr40 or nr90, otherwise tag is just nr. Use with -nr")
    parser.add_argument("-d", "--delete_list", nargs="?", default=None,
                        help="path to a specific list of pdb entries to exclude (e.g. 1okc); one per line") #TODO replace with define_locations
    parser.add_argument("-nobeta", "--exclude_beta_barrels", action='store_true', default=False,
                        help="exclude all beta-barrel proteins from plots. Defined by BETA1 column in pair-file.")

    args = parser.parse_args()

    return args


def main():
    inputs = parse_args()
    print("Input file for chain data is", inputs.chain_file)
    print("Input file for beta-barrel assignment and aligned pairs is", inputs.pair_file)
    print("Writing plots and output data to: ", inputs.output_folder)
    if inputs.nonredundant_chain_list:
        print("Input file containing subset for reducing redundancy is:", inputs.nonredundant_chain_list)
    else:
        print("Creating plots for entire database")

    exclusions = []
    if inputs.delete_list:
        print("Excluding all PDB entries in the provided delete list:", inputs.delete_list)
        exclusions = read_exclusions(inputs.delete_list)
        print(len(exclusions), "deleted chains", exclusions)

    pair_data = read_and_parse_pairs_table(inputs.pair_file, inputs.output_folder)

    tag = ""
    if inputs.exclude_beta_barrels:
        beta_table = pair_data[pair_data['IS_BETA1']]
        beta_list = np.unique(beta_table['PDB_CODE1'])
        exclusions.update(beta_list)
        tag += "alpha_"
        print("Found the following", len(beta_list), "beta barrels to exclude:\n", beta_list)

    table = read_and_parse_chain_table(inputs.chain_file, inputs.output_folder, exclusions)

    if not inputs.nonredundant_chain_list:
        pdbi_list = np.unique(table['PDB_CODE'])  # PDB codes in complete dataset
        print("The full, redundant dataset contains ", len(pdbi_list), "entries")
    else:
        pdbi_list, table = parse_nonredundant_list(inputs.nonredundant_chain_list, table, exclusions)
        if inputs.nr_threshold_label:
            tag += inputs.nr_threshold_label
            tag += "_"
        else:
            tag += "nr_"
    full_name = inputs.output_folder + "/" + tag
    print("Output filenames will be prepended with: ", full_name)

    plot_databasewide_data(table, pdbi_list, full_name)

    # this plot requires pairwise data, which does not apply to non-redundant sets
    MAXDIFF = 4
    if not inputs.nonredundant_chain_list:
        ncomp_stacked = define_num_comparisons_for_stacked_plots(MAXDIFF, pair_data)
        plot_num_struct_align_per_tmregion(MAXDIFF, ncomp_stacked, pair_data, full_name)


if __name__ == '__main__':
    main()
