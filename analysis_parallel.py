from supporting_functions import *
from sklearn.cluster import AgglomerativeClustering
import run_locusts_on_analysis
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42


def cart2pol(x, y):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def con(x, r=[]):
    res = []
    for i in range(2, len(x), 2):
        res.append(r[int(i / 2)] ** 2 - x[i] ** 2 - x[i + 1] ** 2)
    return res


def equilibrium_ENM(posMX, relaxdMX, stiffMX, radialc):
    con1 = {'type': 'eq', 'fun': con, 'args': [radialc]}
    bounds = [[radialc[0], radialc[0]], [0, 0]] + [[None, None]] * (posMX.shape[0] * 2 - 2)

    newposMX = scipy.optimize.minimize(energy_ENM, posMX.flatten(),
                                       args=(relaxdMX, stiffMX),
                                       method='trust-constr',
                                       constraints=[con1]).x.reshape((-1, 2))
    return newposMX


def energy_ENM(posMX, relaxdMX, stiffMX):
    # The argument P is a vector (flattened matrix).
    # We convert it to a matrix here.
    P = posMX.reshape((-1, 2))
    L = relaxdMX  # .reshape(())
    K = stiffMX  # .reshape((-1, 1))
    # We compute the distance matrix.
    D = scipy.spatial.distance_matrix(P, P)
    # The potential energy is the sum of the
    # gravitational and elastic potential energies.
    return .5 * (K * (D - L) ** 2).sum()  # element-wise (Hadamard) multiplication...


def colorbar_index(ncolors, cmap, ntm, cax=None, horiz=False):
    if horiz:
        cmap = cmap_discretize(cmap, ncolors)
        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(-0.5, ncolors + 0.5)
        colorbar = plt.colorbar(mappable, cax=cax, ticks=[], orientation='horizontal', fraction=0.046, pad=0.04)
        colorbar.outline.set_visible(False)
        colorbar.ax.get_yaxis().set_ticks([])
        for j, lab in [(0, '#TM -5'), (5.5, '#TM ({0})'.format(ntm)), (11, '#TM +5')]:
            colorbar.ax.text(j, 5, lab, ha='center', va='center')
    else:
        cmap = cmap_discretize(cmap, ncolors)
        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(-0.5, ncolors + 0.5)
        colorbar = plt.colorbar(mappable, fraction=0.046, pad=0.04)
        colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
        colorbar.set_ticklabels(['#TM -5', '', '', '', '', '#TM ({0})'.format(ntm), '', '', '', '', '#TM +5'])
    return


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

		cmap: colormap instance, eg. cm.jet. 
		N: number of colors.

	Example
		x = resize(arange(100), (5,100))
		djet = cmap_discretize(cm.jet, 5)
		imshow(x, cmap=djet)
	"""

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki]) for i in range(N + 1)]

    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)


def neighborlists(options, sorted_table_fn, an_locations,
                  no_print=False): 
    ### WARNING: ASSUMES summary_table IS SORTED!
    SEQTHR = float(options['PARAMETERS'][('seqid_thr', 'sequence_identity_threshold')])
    STRTHR = float(options['PARAMETERS'][('tmscore_thr', 'tmscore_threshold')])

    pdbi_chs = []  # PDBI_CHs for which a neighbor list has been compiled
    pdbi_chs2 = set()

    numtotneighs = {}

    wref = 'XXX'
    open_files = {}
    with open(sorted_table_fn) as stf:
        for line in stf:
            pdbi1, ch1, pdbi2, ch2, alnlen, seq_seqid, str_seqid, rmsd, tmscore = line.split()
            pdbi_ch1 = pdbi1 + "_" + ch1
            pdbi_ch2 = pdbi2 + "_" + ch2
            pdbi_chs2.add(pdbi_ch2)
            if not no_print and pdbi1 != wref:
                wref = pdbi1
                for pc in open_files:
                    for f in open_files[pc]:
                        f.close()
                open_files = {}
            if pdbi_ch1 not in pdbi_chs:
                pdbi_chs.append(pdbi_ch1)
                numtotneighs[pdbi_ch1] = 0
                if not no_print:
                    seqneigh_fn = an_locations['FSYSPATH']['seqneighs'] + 'seqneigh_' + pdbi_ch1 + '.txt'
                    strneigh_fn = an_locations['FSYSPATH']['strneighs'] + 'strneigh_' + pdbi_ch1 + '.txt'
                    totneigh_fn = an_locations['FSYSPATH']['totneighs'] + 'totneigh_' + pdbi_ch1 + '.txt'
                    segment_fn = an_locations['FSYSPATH']['totneighs'] + 'segment_' + pdbi_ch1 + '.txt'
                    seqneigh_f = open(seqneigh_fn, 'w')
                    strneigh_f = open(strneigh_fn, 'w')
                    totneigh_f = open(totneigh_fn, 'w')
                    segment_f = open(segment_fn, 'w')
                    open_files[pdbi_ch1] = [seqneigh_f, strneigh_f, totneigh_f, segment_f]

            if float(seq_seqid) > SEQTHR or float(tmscore) > STRTHR:
                numtotneighs[pdbi_ch1] += 1
                if not no_print:
                    open_files[pdbi_ch1][2].write("{0}\t{1}\t{2}\n".format(pdbi_ch2, seq_seqid, tmscore))
            if no_print:
                continue
            if float(seq_seqid) > SEQTHR:
                open_files[pdbi_ch1][0].write("{0}\t{1}\n".format(pdbi_ch2, seq_seqid))
            if float(tmscore) > STRTHR:
                open_files[pdbi_ch1][1].write("{0}\t{1}\n".format(pdbi_ch2, tmscore))
            open_files[pdbi_ch1][3].write("{0}\t{1}\t{2}\n".format(pdbi_ch2, seq_seqid, tmscore))

    print("ENTRIES ONLY IN COL2", pdbi_chs2 - set(pdbi_chs))
    print("ENTRIES ONLY IN COL1", set(pdbi_chs) - pdbi_chs2)
    pdbi_chs = sorted(list(set(pdbi_chs)))
    return pdbi_chs, numtotneighs


def cdist(c1, c2):
    return ((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2) ** 0.5


def find_disposition(tm, ccoord, pcoord, placed, j):
    n = len(placed)

    sampled_pcoord = np.zeros((360, 2))
    sampled_ccoord = np.zeros((360, 2))
    for a in range(sampled_pcoord.shape[0]):
        theta = 2 * math.pi * a / sampled_pcoord.shape[0]
        sampled_pcoord[a] = (tm[0][j], theta)
        sampled_ccoord[a] = (tm[0][j] * math.cos(theta), tm[0][j] * math.sin(theta))

    # Sum the potential
    V = np.zeros(360)
    for a in range(V.shape[0]):
        for i in placed:
            #			print(a, i, ccoord[i], sampled_ccoord[a])
            d = cdist(ccoord[i], sampled_ccoord[a])
            d0 = tm[i][j]
            V[a] += (d - d0) ** 2

    am = V.argmin()
    return sampled_pcoord[am], sampled_ccoord[am]


def gif_instr_creator(options, locations, neighbors, names):
    this_name = 'gif_creator'
    tic = time.time()
    script_folder = os.getcwd() + '/'
    fnull = open(os.devnull, 'w')

    pymolpath = options['pymol_path']

    for pdbname in sorted(list(names.keys())):
        pdb_filename = locations['']
        pml_filename = locations['FSYS']['images_pml_whole'] + pdbname + '.pml'
        fig_filename = locations['FSYS']['images_figs_whole'] + pdbname + '.png'
        pml_file = open(pml_filename, 'w')
        pml_file.write(create_whole_png_text(pdb_filename, fig_filename))
        pml_file.close()
        os.system("cd {0}; {1} -ucq {2}; cd -".format(locations['analysis'], pymolpath, pml_filename))

    for struct in sorted(list(neighbors.keys())):
        pml_filename = locations['FSYS']['images_pml_chain'] + struct + '.pml'
        fig_filename = locations['FSYS']['images_figs_chain'] + struct + '.png'
        pdbname = struct[:4]
        chain = struct[5]
        pml_file = open(pml_filename, 'w')
        pml_file.write(create_chain_png_text(pdb_filename, chain, fig_filename))
        pml_file.close()
        os.system("cd {0}; {1} -ucq {2}; cd -".format(locations['analysis'], pymolpath, pml_filename))

    toc = time.time()
    print_log(this_name, "Time spent: {0}".format(toc - tic))


def gen_cartesian_randomphi(distV):
    # Returns a vector of length len(distV)+1 od 2D polar coords, with the first one being (0,0)

    cres = np.zeros((len(distV), 2))
    pres = np.zeros((len(distV), 2))
    for i in range(len(distV)):
        randphi = 2 * np.pi * random.random()
        if distV[i] == 0:
            randphi = 0
        pres[i][:] = (distV[i], randphi)
        cres[i][:] = pol2cart(distV[i], randphi)
    return cres, pres


def polar(options, locations, entry_list, ntms, only_write=False, do_only=None, do_eps=False, do_pdf=False):
    this_name = 'polar'  ## CHANGE

    print("! Creating polar plot", do_pdf)
    print("NTMS", ntms)

    # Figure parameters
    figsize = (8, 8)
    points_per_inch = 72.272  # Conversion factor
    clickmap = {}  # Dict containing for each entry a list of 3 same-length lists with names, xcoord and ycoord of the pixels

    N = len(entry_list)
    # Extract TMscore and sequence identity AS METRICS
    tm = []  # Includes self as first
    tmlabels = []  # Includes self as first

    # Record all the TMscores from all entries
    for pdbi_ch1 in entry_list:
        tm.append([1])
        tmlabels.append([pdbi_ch1])
        strneigh_fn = locations['FSYSPATH']['strneighs'] + 'strneigh_' + pdbi_ch1 + '.txt'
        with open(strneigh_fn) as strneigh_f:
            for line in strneigh_f:
                pdbi_ch2, tmscore = line.split()
                tm[-1].append(float(tmscore))
                tmlabels[-1].append(pdbi_ch2)

    # Main loop
    if do_only:
        en = [(b, c) for b, c in enumerate(entry_list) if c in do_only]
        print("DO ONLY", do_only)
        print("EN", en)
    else:
        en = list(enumerate(entry_list))
    for i1, pdbi_ch1 in en:
        # Build the TM matrix
        subN = len(tm[i1])
        subtm = np.ones((subN, subN))  # indexes relative to pdbi_ch1
        print("ITSELF", pdbi_ch1)
        print("NEIGHBORS", tmlabels[i1])
        print("NUMBER OF NEIGHBORS", subN)
        for j1, tm1 in enumerate(tm[i1]):  # [(ii+1, x) for ii, x in enumerate(tm[i1])]:
            # j1 in [1, len(tm[i1])+1)
            # subtm[0, j1] = subtm[j1, 0] = tm1
            j1_ind = entry_list.index(tmlabels[i1][j1])
            for j2, tmlabel2 in [(h, c) for h, c in enumerate(tmlabels[i1]) if
                                 h > j1]:  # nonj2 is an index relative to a neighbor of pdbi_ch1, must be converted
                if tmlabel2 in tmlabels[j1_ind]:
                    j2_ind = tmlabels[j1_ind].index(tmlabel2)
                    subtm[j1, j2] = subtm[j2, j1] = tm[j1_ind][j2_ind]
        subtm = 1. - subtm

        # cpos0, ppos0 = gen_cartesian_randomphi(subtm[0])  # pos0.shape == subtm.shape
        ccoords, pcoords = -np.ones((subtm.shape[0], 2)), -np.ones((subtm.shape[0], 2))

        build_matrix = 10 * np.ones(subtm.shape)
        build_matrix[0] = np.copy(subtm[0])
        placed = set()
        placed0 = set()  # set([0])
        while len(placed) + len(placed0) < subtm.shape[0]:
            # Flatten the matrix and decide an order for putting points
            linear_form = []
            for nn in range(build_matrix.shape[0]):
                for mm in range(build_matrix.shape[1]):
                    linear_form.append((nn, mm, build_matrix[nn][mm]))
            linear_form = sorted(linear_form,
                                 key=lambda x: (x[2], x[0]))  # Order: from the closest to the farthest from the center
            for tr in linear_form:
                if tr[0] == 0 and tr[2] == 0 and tr[1] not in placed0:
                    placed0.add(tr[1])
                    pcoords[tr[1]][0] = 0.0
                    pcoords[tr[1]][1] = 0.0
                    ccoords[tr[1]][0] = 0.0
                    ccoords[tr[1]][1] = 0.0
                    continue
                elif tr[1] not in placed and tr[1] not in placed0:
                    pcoords[tr[1]], ccoords[tr[1]] = find_disposition(build_matrix, ccoords, pcoords, placed, tr[1])
                    build_matrix[tr[1]] = np.copy(subtm[tr[1]])
                    placed.add(tr[1])
                    break

        ntm_diff = []  # Difference in number of TM regions
        coord_fn = locations['FSYSPATH']['polar_data'] + '{0}_neigh.txt'.format(pdbi_ch1)
        pdbi1, ch1 = pdbi_ch1.split('_')
        with open(coord_fn, 'w') as coord_f:
            coord_f.write('#PDBid\t{1:8}\t{2:8}\t{3}\n'.format(pdbi_ch1, 'rho(=TM)', 'theta', '#MemReg'))
            for i in range(ccoords.shape[0]):
                pdbi2, ch2 = tmlabels[i1][i].split('_')
                coord_f.write('{0}\t{1:8.5f}\t{2:8.5f}\t{3}\n'.format(tmlabels[i1][i], pcoords[i][0], pcoords[i][1],
                                                                      ntms[pdbi2][ch2]))
                ntm_diff.append(0.5 - 0.1 * (ntms[pdbi1][ch1] - ntms[pdbi2][ch2]))

        rhos, thetas = pcoords[:, 0], pcoords[:, 1]

        # Plot
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
        ax.axes.get_xaxis().set_ticks([])
        plt.yticks([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4], ['', '0.9', '', '0.8', '', '0.7', '', '0.6'])
        ax.set_axisbelow(True)
        im = ax.scatter(thetas, rhos, c=ntm_diff, cmap='coolwarm_r', edgecolors='k')
        im.set_clim(0.0, 1.0)
        cmap = plt.get_cmap('coolwarm_r')
        colorbar_index(ncolors=11, cmap=cmap, ntm=ntms[pdbi1][ch1], horiz=True)
        ax.set_rmin(0.0)
        ax.set_rmax(0.4)
        fig.canvas.draw()
        figname = locations['FSYSPATH']['polar_figs'] + 'p_' + pdbi_ch1 + ".png"
        plt.savefig(figname, dpi=fig.dpi)
        if do_eps:
            plt.savefig(figname[:-4] + '.eps', dpi=fig.dpi)
        if do_pdf:
            plt.savefig(figname[:-4] + '.pdf')
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
        ax.axes.get_xaxis().set_ticks([])
        plt.yticks([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4], ['', '0.9', '', '0.8', '', '0.7', '', '0.6'])
        ax.plot(thetas, rhos, 'bo')
        im = ax.scatter(thetas, rhos, c=ntm_diff, cmap='coolwarm_r', edgecolors='k')
        im.set_clim(0.0, 1.0)

        cmap = plt.get_cmap('coolwarm_r')
        colorbar_index(ncolors=11, cmap=cmap, ntm=ntms[pdbi1][ch1], horiz=True)

        ax.set_rmin(0.0)
        ax.set_rmax(0.4)
        fig.canvas.draw()

        pix_xy = ax.transData.transform(np.vstack([thetas, rhos]).T)  # Take pixels
        pix_x, pix_y = pix_xy.T
        width, height = fig.canvas.get_width_height()
        pix_y = height - pix_y
        clickmap = list(zip(tmlabels[i1], pix_x, pix_y))
        clickmapfilename = locations['FSYSPATH']['polar_maps'] + '{0}_cmap.txt'.format(pdbi_ch1)
        clickmapfile = open(clickmapfilename, 'w')

        print(clickmap)

        groups = set()
        for sname, px, py in clickmap:
            added_to_group = False
            for g in groups:
                not_in_this_group = False
                for elem in g:
                    d = ((elem[1] - px) ** 2 + (elem[2] - py) ** 2) ** 0.5
                    if (d > 5):
                        not_in_this_group = True
                        break
                if not not_in_this_group:
                    gg = set()
                    gg.add(g)
                    groups = groups - gg
                    g = set(g)
                    g.add((sname, px, py))
                    groups.add(frozenset(g))
                    added_to_group = True
                    break
            if not added_to_group:
                fs = set()
                fs.add((sname, px, py))
                groups.add(frozenset(fs))

        clickmapfilename = locations['FSYSPATH']['polar_maps'] + '{0}_cmap.txt'.format(pdbi_ch1)
        with open(clickmapfilename, 'w') as cmf:
            cmf.write("[\n")
            for ng, g in enumerate(groups):
                cx = int(sum([n[1] for n in list(g)]) / len(list(g)))
                cy = int(sum([n[2] for n in list(g)]) / len(list(g)))
                ns = [n[0] for n in list(g)]
                if ng == len(groups) - 1:
                    cmf.write("{{\"X\":{0},\"Y\":{1},\"radius\":3,\"chains\":\"{2}\"}}\n".format(cx, cy, "".join(
                        [x + "," for x in ns[:-1]]) + ns[-1]))
                else:
                    cmf.write("{{\"X\":{0},\"Y\":{1},\"radius\":3,\"chains\":\"{2}\"}},\n".format(cx, cy, "".join(
                        [x + "," for x in ns[:-1]]) + ns[-1]))
            cmf.write("]\n")


def scatterplots_parallel(an_locations, entry_list, ntms, pkl_fn, do_eps=False, do_pdf=False):
    print("!About to open dictionary:", pkl_fn)
    print("!Creating scatter plots:", do_pdf)
    with open(pkl_fn, 'rb') as f:
        packing = pickle.load(f)
    (fig, (rect_scatter, rect_histx, rect_histy, binwidth), (si_tot, tm_tot, labels_tot)) = packing
    labelorder = [x[0] for x in labels_tot]
    fig.clf()
    print("ENTRY LIST", entry_list, len(entry_list))
    for i1, pdbi_ch1 in enumerate(entry_list):
        fig.clf()
        print("INIT CYCLE", i1, pdbi_ch1, time.time())

        labid = labelorder.index(pdbi_ch1)
        XX = si_tot[labid]
        YY = tm_tot[labid]
        labels = labels_tot[labid]

        pdbi1, ch1 = pdbi_ch1.split('_')
        layer_fn = an_locations['FSYSPATH']['densityscatter_figs'] + 'ds_' + pdbi_ch1 + '_layer.png'
        figfilename = an_locations['FSYSPATH']['densityscatter_figs'] + 'ds_' + pdbi_ch1 + '.png'

        ax = fig.add_axes(rect_scatter)
        axHistx = fig.add_axes(rect_histx)
        axHisty = fig.add_axes(rect_histy)

        levels = np.array(
            [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
             0.95, 0.99])
        cmap = copy.copy(cm.get_cmap('Blues'))

        print("SCATTER", time.time())

        plt.locator_params(nbins=6)
        tt = ntms[pdbi1][ch1]  # table[table["PDB1"] == pdbi]["NTM1"].iloc[0].astype(int).values[0]
        ndif = []
        for i2, pdbi_ch2 in enumerate(labels):
            pdbi2, ch2 = pdbi_ch2.split('_')
            ndif.append(ntms[pdbi2][ch2] - tt)
        ax.set_axisbelow(True)
        cmap.set_under("magenta")
        cmap.set_over("yellow")
        cmap.set_bad("red")
        im = ax.scatter(XX, YY, c=(np.array(ndif) + 5) / 10, cmap='coolwarm_r', edgecolors='k', zorder=1)
        im.set_clim(0.0, 1.0)
        cbaxes = fig.add_axes([0.1, 0.01, 0.8, 0.03])
        cmap = plt.get_cmap('coolwarm_r')
        colorbar_index(ncolors=11, cmap=cmap, ntm=ntms[pdbi1][ch1], cax=cbaxes, horiz=True)

        print("CANVAS DRAW", time.time())

        fig.canvas.draw()
        ax.set_xlabel("Sequence Identity")
        ax.set_ylabel("Structure Similarity")
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])

        print("XHIST", time.time())

        bins = np.arange(0, 1 + binwidth, binwidth)
        plt.locator_params(nbins=4)
        axHistx.hist(XX, bins=bins, color='pink')
        axHistx.set_xlim(ax.get_xlim())
        axHistx.set_xticks([])

        print("YHIST", time.time())

        plt.locator_params(nbins=4)
        axHisty.hist(YY, bins=bins, orientation='horizontal', color='pink')
        axHisty.set_ylim(ax.get_ylim())
        axHisty.set_yticks([])

        print("TRANSFORMS", time.time())

        pix_xy = ax.transData.transform(np.vstack([XX, YY]).T)
        pix_x, pix_y = pix_xy.T
        width, height = fig.canvas.get_width_height()
        pix_y = height - pix_y

        snames = labels  # [s for s in sorted(list(table[struct].keys()))]
        clickmap = list(zip(snames, pix_x, pix_y))
        clickmapfilename = an_locations['FSYSPATH']['densityscatter_maps'] + '{0}_cmap.txt'.format(pdbi_ch1)
        clickmapfile = open(clickmapfilename, 'w')

        print("GROUPS", time.time())

        groups = set()
        for sname, px, py in clickmap:
            added_to_group = False
            for g in groups:
                not_in_this_group = False
                for elem in g:
                    d = ((elem[1] - px) ** 2 + (elem[2] - py) ** 2) ** 0.5
                    if (d > 5):
                        not_in_this_group = True
                        break
                if not not_in_this_group:
                    gg = set()
                    gg.add(g)
                    groups = groups - gg
                    g = set(g)
                    g.add((sname, px, py))
                    groups.add(frozenset(g))
                    added_to_group = True
                    break
            if not added_to_group:
                fs = set()
                fs.add((sname, px, py))
                groups.add(frozenset(fs))

        print("WRITE CLICKMAP", time.time())

        clickmapfile.write("[\n")
        for ng, g in enumerate(groups):
            cx = int(sum([n[1] for n in list(g)]) / len(list(g)))
            cy = int(sum([n[2] for n in list(g)]) / len(list(g)))
            ns = [n[0] for n in list(g)]
            if ng == len(groups) - 1:
                clickmapfile.write("{{\"X\":{0},\"Y\":{1:5d},\"radius\":3,\"chains\":\"{2}\"}}\n".format(cx, cy,
                                                                                                         "".join(
                                                                                                             [x + ","
                                                                                                              for x in
                                                                                                              ns[
                                                                                                              :-1]]) +
                                                                                                         ns[-1]))
            else:
                clickmapfile.write("{{\"X\":{0},\"Y\":{1:5d},\"radius\":3,\"chains\":\"{2}\"}},\n".format(cx, cy,
                                                                                                          "".join(
                                                                                                              [x + ","
                                                                                                               for x in
                                                                                                               ns[
                                                                                                               :-1]]) +
                                                                                                          ns[-1]))
        clickmapfile.write("]\n")
        clickmapfile.close()

        fig.savefig(layer_fn, dpi=fig.dpi, transparent=True)
        print("!Creating png:", layer_fn)
        if do_eps:
            fig.savefig(layer_fn[:-4] + '.eps', dpi=fig.dpi, transparent=True)
            print("!Creating eps:", layer_fn[:-4] + '.eps')
        if do_pdf:
            fig.savefig(layer_fn[:-4] + '.pdf', dpi=fig.dpi, transparent=True)
            print("!Creating pdf layer:", layer_fn[:-4] + '.pdf')
            print("!, will need to be combined with background:", an_locations['SYSFILES']['background'])
        plt.close(fig)
        ax.remove()

        print(" ".join(["composite", layer_fn, an_locations['SYSFILES']['background'], figfilename]))
        msg = subprocess.run(["composite", layer_fn, an_locations['SYSFILES']['background'], figfilename],
                             stderr=subprocess.PIPE).stderr.decode('utf8')
        print(msg)
        print("PNG", figfilename)
        os.remove(layer_fn)


def scatterplots_background(an_locations, entry_list, pkl_fn="sp.pkl", do_pdf=True):
    N = len(entry_list)
    # Extract TMscore and sequence identity AS METRICS
    tm = []  # Includes self as first
    si = []
    labels = []  # Includes self as first

    # Record all the TMscores from all entries
    for pdbi_ch1 in entry_list:
        tm.append([1])
        si.append([1])
        labels.append([pdbi_ch1])
        totneigh_fn = an_locations['FSYSPATH']['totneighs'] + 'segment_' + pdbi_ch1 + '.txt'
        with open(totneigh_fn) as totneigh_f:
            for line in totneigh_f:
                pdbi_ch2, seqid, tmscore = line.split()
                tm[-1].append(float(tmscore))
                si[-1].append(float(seqid))
                labels[-1].append(pdbi_ch2)

    xmin, xmax = 0, 1
    ymin, ymax = 0, 1
    xedges = np.mgrid[0:1:200j]
    yedges = np.mgrid[0:1:200j]
    H, xedges, yedges = np.histogram2d([item for sublist in tm for item in sublist],
                                       [item for sublist in si for item in sublist], bins=(xedges, yedges))

    # Take log of histogram, adjust -inf to 0
    H[H == 0] = 0.00000001
    H = np.log(H)
    H[H > 1] = 1

    # Build KDE on log-histogram, by creating a fake set of data
    # on a 4-fold graph where the real graph is the lower left quadrant
    DDx = []
    DDy = []
    for x in range(len(H)):
        for y in range(len(H[x])):
            for i in range(math.ceil(H[x][y])):
                DDx.append(x / len(H))
                DDy.append(y / len(H[x]))
                DDx.append(2 - (x / len(H)))
                DDy.append(y / len(H[x]))
                DDx.append(x / len(H))
                DDy.append(2 - (y / len(H[x])))
                DDx.append(2 - (x / len(H)))
                DDy.append(2 - (y / len(H[x])))
    DDx = np.array(DDx)
    DDy = np.array(DDy)
    xmin = DDx.min()
    xmax = DDx.max()
    ymin = DDy.min()
    ymax = DDy.max()
    X, Y = np.mgrid[0:2:400j, 0:2:400j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([DDx, DDy])
    kernel = scipy.stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions), X.shape)
    Z = Z.T

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    binwidth = 0.05

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    fig_not = plt.figure(1, figsize=(8, 8))
    ax_not = plt.axes(rect_scatter)
    ax_not.set_xlim([0, 1])
    ax_not.set_ylim([0, 1])
    ZZ = np.array(Z)
    levels = np.array(
        [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
         0.99])
    cfset = ax_not.contourf(X, Y, ZZ, levels=levels, cmap="Blues")
    cset = ax_not.contour(X, Y, ZZ, levels=levels, colors='k')
    plt.locator_params(nbins=6)
    ax_not.imshow(np.rot90(ZZ), cmap="Blues", extent=[0, 2, 0, 2], aspect=1)
    ax_not.set_xlabel("Sequence Identity")
    ax_not.set_ylabel("Structure Similarity")
    fig_not.savefig(an_locations['SYSFILES']['background'], dpi=fig_not.dpi)
    fig_not.savefig(an_locations['SYSFILES']['background'][:-4] + '.eps', dpi=fig_not.dpi)
    fig_not.savefig(an_locations['SYSFILES']['background'][:-4] + '.pdf', dpi=fig_not.dpi)

    packing = (fig_not, (rect_scatter, rect_histx, rect_histy, binwidth), (si, tm, labels))
    with open(pkl_fn, 'wb') as f:
        pickle.dump(packing, f)
    plt.close(fig_not)

    return


def scatterplots(options, locations, entry_list, ntms, an_locations, do_only=None, do_eps=False, do_pdf=False):
    this_name = scatterplots.__name__

    # Figure parameters
    figsize = (8, 8)
    points_per_inch = 72.272  # Conversion factor
    clickmap = {}  # Dict containing for each entry a list of 3 same-length lists with names, xcoord and ycoord of the pixels

    N = len(entry_list)
    # Extract TMscore and sequence identity AS METRICS
    tm = []  # Includes self as first
    si = []
    labels = []  # Includes self as first

    # Record all the TMscores from all entries
    for pdbi_ch1 in entry_list:
        tm.append([1])
        si.append([1])
        labels.append([pdbi_ch1])
        totneigh_fn = an_locations['FSYSPATH']['totneighs'] + 'segment_' + pdbi_ch1 + '.txt'
        with open(totneigh_fn) as totneigh_f:
            for line in totneigh_f:
                pdbi_ch2, seqid, tmscore = line.split()
                tm[-1].append(float(tmscore))
                si[-1].append(float(seqid))
                labels[-1].append(pdbi_ch2)

    xmin, xmax = 0, 1
    ymin, ymax = 0, 1
    xedges = np.mgrid[0:1:200j]
    yedges = np.mgrid[0:1:200j]
    H, xedges, yedges = np.histogram2d([item for sublist in tm for item in sublist],
                                       [item for sublist in si for item in sublist], bins=(xedges, yedges))

    # Take log of histogram, adjust -inf to 0
    H[H == 0] = 0.00000001
    H = np.log(H)
    H[H > 1] = 1

    # Build KDE on log-histogram, by creating a fake set of data
    # on a 4-fold graph where the real graph is the lower left quadrant
    DDx = []
    DDy = []
    for x in range(len(H)):
        for y in range(len(H[x])):
            for i in range(math.ceil(H[x][y])):
                DDx.append(x / len(H))
                DDy.append(y / len(H[x]))
                DDx.append(2 - (x / len(H)))
                DDy.append(y / len(H[x]))
                DDx.append(x / len(H))
                DDy.append(2 - (y / len(H[x])))
                DDx.append(2 - (x / len(H)))
                DDy.append(2 - (y / len(H[x])))
    DDx = np.array(DDx)
    DDy = np.array(DDy)
    X, Y = np.mgrid[0:2:400j, 0:2:400j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([DDx, DDy])
    kernel = scipy.stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions), X.shape)
    Z = Z.T

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    binwidth = 0.05

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    fig_not = plt.figure(1, figsize=(8, 8))
    ax_not = plt.axes(rect_scatter)
    ax_not.set_xlim([0, 1])
    ax_not.set_ylim([0, 1])
    ZZ = np.array(Z)
    levels = np.array(
        [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
         0.99])
    cfset = ax_not.contourf(X, Y, ZZ, levels=levels, cmap="Blues")
    cset = ax_not.contour(X, Y, ZZ, levels=levels, colors='k')
    plt.locator_params(nbins=6)
    ax_not.imshow(np.rot90(ZZ), cmap="Blues", extent=[0, 2, 0, 2], aspect=1)
    ax_not.set_xlabel("Sequence Identity")
    ax_not.set_ylabel("Structure Similarity")
    fig_not.savefig(an_locations['FSYSPATH']['databasewide'] + 'densitylines.png', dpi=fig_not.dpi)
    fig_not.savefig(an_locations['FSYSPATH']['databasewide'] + 'densitylines.eps', dpi=fig_not.dpi)

    buf = io.BytesIO()
    pickle.dump(fig_not, buf)

    plt.close(fig_not)

    for i1, pdbi_ch1 in enumerate(entry_list):
        pdbi1, ch1 = pdbi_ch1.split('_')
        figfilename = an_locations['FSYSPATH']['densityscatter_figs'] + 'ds_' + pdbi_ch1 + '.png'
        ZZ = np.array(Z)
        buf.seek(0)
        fig = pickle.load(buf)
        ax = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)

        levels = np.array(
            [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
             0.95, 0.99])
        cmap = copy.copy(cm.get_cmap('Blues'))

        plt.locator_params(nbins=6)
        tt = ntms[pdbi1][ch1]  # table[table["PDB1"] == pdbi]["NTM1"].iloc[0].astype(int).values[0]
        XX = si[i1]
        YY = tm[i1]
        ndif = []
        for i2, pdbi_ch2 in enumerate(labels[i1]):
            pdbi2, ch2 = pdbi_ch2.split('_')
            ndif.append(ntms[pdbi2][ch2] - tt)
        ax.set_axisbelow(True)
        cmap.set_under("magenta")
        cmap.set_over("yellow")
        cmap.set_bad("red")
        im = ax.scatter(XX, YY, c=(np.array(ndif) + 5) / 10, cmap='coolwarm_r', edgecolors='k', zorder=1)
        im.set_clim(0.0, 1.0)
        cbaxes = fig.add_axes([0.1, 0.01, 0.8, 0.03])
        cmap = plt.get_cmap('coolwarm_r')
        colorbar_index(ncolors=11, cmap=cmap, ntm=ntms[pdbi1][ch1], cax=cbaxes, horiz=True)

        fig.canvas.draw()
        ax.set_xlabel("Sequence Identity")
        ax.set_ylabel("Structure Similarity")

        bins = np.arange(0, 1 + binwidth, binwidth)
        plt.locator_params(nbins=4)
        axHistx.hist(XX, bins=bins, color='pink')
        axHistx.set_xlim(ax.get_xlim())
        axHistx.set_xticks([])

        plt.locator_params(nbins=4)
        axHisty.hist(YY, bins=bins, orientation='horizontal', color='pink')
        axHisty.set_ylim(ax.get_ylim())
        axHisty.set_yticks([])

        pix_xy = ax.transData.transform(np.vstack([XX, YY]).T)
        pix_x, pix_y = pix_xy.T
        width, height = fig.canvas.get_width_height()
        pix_y = height - pix_y

        snames = labels[i1]  # [s for s in sorted(list(table[struct].keys()))]
        clickmap = list(zip(snames, pix_x, pix_y))
        clickmapfilename = an_locations['FSYSPATH']['densityscatter_maps'] + '{0}_cmap.txt'.format(pdbi_ch1)
        clickmapfile = open(clickmapfilename, 'w')

        groups = set()
        for sname, px, py in clickmap:
            added_to_group = False
            for g in groups:
                not_in_this_group = False
                for elem in g:
                    d = ((elem[1] - px) ** 2 + (elem[2] - py) ** 2) ** 0.5
                    if (d > 5):
                        not_in_this_group = True
                        break
                if not not_in_this_group:
                    gg = set()
                    gg.add(g)
                    groups = groups - gg
                    g = set(g)
                    g.add((sname, px, py))
                    groups.add(frozenset(g))
                    added_to_group = True
                    break
            if not added_to_group:
                #                print("New group created!")
                fs = set()
                fs.add((sname, px, py))
                groups.add(frozenset(fs))

        clickmapfile.write("[\n")
        for ng, g in enumerate(groups):
            cx = int(sum([n[1] for n in list(g)]) / len(list(g)))
            cy = int(sum([n[2] for n in list(g)]) / len(list(g)))
            ns = [n[0] for n in list(g)]
            if ng == len(groups) - 1:
                clickmapfile.write("{{\"X\":{0},\"Y\":{1:5d},\"radius\":3,\"chains\":\"{2}\"}}\n".format(cx, cy,
                                                                                                         "".join(
                                                                                                             [x + ","
                                                                                                              for x in
                                                                                                              ns[
                                                                                                              :-1]]) +
                                                                                                         ns[-1]))
            else:
                clickmapfile.write("{{\"X\":{0},\"Y\":{1:5d},\"radius\":3,\"chains\":\"{2}\"}},\n".format(cx, cy,
                                                                                                          "".join(
                                                                                                              [x + ","
                                                                                                               for x in
                                                                                                               ns[
                                                                                                               :-1]]) +
                                                                                                          ns[-1]))
        clickmapfile.write("]\n")

        clickmapfile.close()

        fig.savefig(figfilename, dpi=fig.dpi)
        if do_eps:
            fig.savefig(figfilename[:-4] + '.eps', dpi=fig.dpi)
        if do_pdf:
            fig.savefig(figfilename[:-4] + '.pdf', dpi=fig.dpi)
        plt.close(fig)


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    window_size = np.abs(np.int(window_size))
    order = np.abs(np.int(order))
    if window_size % 2 != 1 or window_size < 1:
        window_size += 1
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


def smoothen(l, k):
    k1 = int(k / 2)
    k2 = k - k1
    newl = [l[0]] * k1 + l + [l[-1]] * k2

    res = []
    for i in range(len(l)):
        res.append(sum(newl[i:i + k]) / k)
    return res


def straln_profiles(locations, entry_list, str_data, an_locations, do_eps=False, do_pdf=False):
    import matplotlib

    cmap = matplotlib.cm.get_cmap('hsv')
    colors = colorblind_palette()
    print("COLORS", colors)
    print(entry_list)
    dgts = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
    for pdbi_ch1 in entry_list:
        if os.path.exists(an_locations['FSYSPATH']['distributions_figs'] + 'distr_' + pdbi_ch1 + '.png'):
            continue
        pdbi1, ch1 = pdbi_ch1[:4], pdbi_ch1[5]
        ref_ch1 = ch1
        if ch1 in str_data[pdbi1]['ENCOMPASS']['structure']['redundant_chains']:
            ref_ch1 = str_data[pdbi1]['ENCOMPASS']['structure']['redundant_chains'][ch1]
            straln_fn = an_locations['FSYSPATH']['distributions_data'] + '/distr_' + pdbi1 + "_" + ref_ch1 + '.txt'
        else:
            straln_fn = an_locations['FSYSPATH']['distributions_data'] + '/distr_' + pdbi_ch1 + '.txt'
        anc = {}
        print(straln_fn)
        if not os.path.exists(straln_fn):
            print(pdbi_ch1, "NOT THERE", straln_fn)
            ref_ch1 = ch1
            straln_fn = an_locations['FSYSPATH']['distributions_data'] + '/distr_' + pdbi_ch1 + '.txt'
        if not os.path.exists(straln_fn):
            print(pdbi_ch1, "NOT THERE AGAIN", straln_fn)
            continue

        with open(straln_fn) as straln_f:
            for line in straln_f:
                if not line.strip():
                    continue
                fields = line.split()
                if 'BEGIN' in line:
                    pdbi_ch2 = fields[4]
                    if pdbi_ch2 not in anc:
                        anc[pdbi_ch2] = {}
                elif 'END' in line:
                    continue
                else:
                    compl = ''
                    if fields[0][-1] not in dgts:
                        nri = int(fields[0][:-1])
                        compl = fields[0][-1]
                    else:
                        compl = ''
                        nri = int(fields[0])
                    anc[pdbi_ch2][(nri, compl)] = float(fields[1])

        print(pdbi_ch1, "LEN ANC", len(anc))

        polarc_fn = an_locations['FSYSPATH']['polar_data'] + '{0}_neigh.txt'.format(pdbi_ch1)
        coords = {}
        pdbis = set()
        with open(polarc_fn) as polarc_f:
            for line in polarc_f:
                if not line.strip() or line.strip().startswith('#'):
                    continue
                pdbi_ch2, rho, theta, _ = line.split()
                if float(rho) > 0.001 and pdbi_ch2 in anc:
                    coords[pdbi_ch2] = pol2cart(float(rho), float(theta))
                    pdbis.add(pdbi_ch2)
        if not pdbis:
            print(pdbi_ch1, "ALL ENTRIES ARE IDENTICAL", polarc_fn)
            continue
        print(pdbi_ch1, "POLARC", polarc_fn)
        print(pdbi_ch1, "PDBIS", pdbis)
        Nall = len(pdbis)
        groups = []
        CLUSTER_RADIUS = 0.3
        X = []
        X_labels = []
        for pdbi_ch2 in sorted(list(pdbis)):
            X_labels.append(pdbi_ch2)
            X.append(coords[pdbi_ch2])
        if len(X_labels) > 1:
            clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=CLUSTER_RADIUS,
                                                 compute_full_tree=True).fit(X)
            for ix in range(max(clustering.labels_) + 1):
                group = []
                for iix, il in enumerate(clustering.labels_):
                    if il == ix:
                        group.append(X_labels[iix])
                groups.append(group)
        elif len(X_labels) == 1:
            group = X_labels[:]
            groups.append(group)

        lines = []
        ich1 = str_data[pdbi1]['ENCOMPASS']['structure']['kchains'].index(ref_ch1)
        iress = []
        for x in str_data[pdbi1]['ENCOMPASS']['structure']['chains'][ich1]['residues']:
            iress.append((x[0][0], x[0][1]))
        for group in groups:
            # print("GROUP", group)
            ganc = {}
            norm = {}
            for pdbi_ch2 in group:
                for ires in anc[pdbi_ch2]:
                    if ires not in ganc:
                        ganc[ires] = 0
                        norm[ires] = 0
                    ganc[ires] += anc[pdbi_ch2][ires]
                    norm[ires] += 1
            for ires in ganc:
                ganc[ires] /= norm[ires]
            segy = []
            segx = []
            segments = []
            for ires in sorted([x for x in ganc], key=lambda x: (x[0], x[1])):
                # print("PROBE", ires)
                if ires not in iress:
                    print("WARNING: ires not found", pdbi_ch1, ires)
                    continue
                idx = iress.index(ires)
                if segx:
                    if idx - previdx > 1:  # MIGLIORA QUI, CI POSSONO COMUNQUE ESSERE DEI BUCHI
                        segments.append((segx, segy))
                        segx = []
                        segy = []
                segx.append(ires)
                previdx = idx
                segy.append(ganc[ires])
            if segx:
                segments.append((segx, segy))
            lines.append(segments)

        if not lines:
            print(pdbi_ch1, "NO LINES! WHY??")
            print(pdbi_ch1, "ITS GROUPS:", groups)
            print(pdbi_ch1, "ITS FRTMALIGN RES", sorted([x for x in ganc], key=lambda x: (x[0], x[1])))
            print(pdbi_ch1, "ITS RECORDED RES", iress)
            continue
        try:
            xmin = min([int(x[0][0][0][0]) for x in lines])
            xmax = max([int(x[-1][0][-1][0]) for x in lines])
        except:
            print(pdbi_ch1, "ERROR HERE")
            print(pdbi_ch1, lines)
            print(pdbi_ch1, "ITS FRTMALIGN RES", sorted([x for x in ganc], key=lambda x: (x[0], x[1])))
            print(pdbi_ch1, "ITS RECORDED RES", iress)
            continue
        fig, ax = plt.subplots(figsize=(max(10, 10 * (xmax - xmin) / 150), 8))
        ax.set_xlim([xmin, xmax])
        ncol = len(lines)
        if ncol > 12:
            print("WARNING: there are {0} lines and only 12 colors!".format(ncol))
        colors = colorblind_palette(n=ncol)
        for il, segments in enumerate(lines):
            for segx, segy in segments:
                thickness = max(10 * len(groups[il]) / Nall, 1)
                # print(segx, segy)
                X = [x[0] for x in segx]
                Y = smoothen(segy, 10)
                ax.plot(X, Y, c=colors[il % len(colors)], linewidth=thickness)
        ax.set_xlim([xmin, xmax])
        # ymax = max(Y)
        ax.set_ylim([0, 10])
        ax.set_xlabel('Residue Number')
        ax.set_ylabel(r'Distance ($\AA$)')
        ax.set_title(r'$C_{\alpha}$ distances of all other topologically related chains')

        pdbi1, ch1 = pdbi_ch1.split('_')
        if pdbi1 in str_data and ch1 in str_data[pdbi1]['ENCOMPASS']['structure']['ktmchains']:
            ich1 = str_data[pdbi1]['ENCOMPASS']['structure']['ktmchains'].index(ch1)
            for ini, fin in str_data[pdbi1]['ENCOMPASS']['structure']['chains'][ich1]['TM_regions'][
                'TM_regions_extrema']:
                plt.axvspan(ini[0][0], fin[0][0], facecolor='0.1', alpha=0.2)

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        inset_ax = inset_axes(ax,
                              width="20%",  # width = 30% of parent_bbox
                              height="20%",  # height : 1 inch
                              axes_class=matplotlib.projections.get_projection_class('polar')
                              )
        for igroup, group in enumerate(groups):
            xs, ys = [], []
            for pdbi_ch2 in group:
                r, t = cart2pol(coords[pdbi_ch2][0], coords[pdbi_ch2][1])
                xs.append(t)
                ys.append(r)
            inset_ax.scatter(xs, ys, c=colors[igroup % len(colors)], s=1)
        inset_ax.set_rmax(0.4)
        inset_ax.set_xticklabels([])
        inset_ax.set_yticklabels([])

        figfilename = an_locations['FSYSPATH']['distributions_figs'] + 'distr_' + pdbi_ch1 + '.png'
        plt.savefig(figfilename)
        if do_eps:
            plt.savefig(figfilename[:-4] + '.eps')
        if do_pdf:
            plt.savefig(figfilename[:-4] + '.pdf')
        plt.close(fig)


def straln_profiles_preprocess(options, locations):
    # ASCII codes of letters and numbers
    ords = [i for i in range(ord('a'), ord('z') + 1)] + [i for i in range(ord('0'), ord('9') + 1)]
    for i in ords:
        # Collect filenames with the wanted pattern and .txt.gz extension
        atomzip_fns = [y.decode('utf8').split()[0] for y in subprocess.Popen(
            "ls {0}/stralns_?{1}??_?.ATOM.txt.gz".format(locations['FSYSPATH']['stralns'], chr(i)), shell=True,
            stdout=subprocess.PIPE).stdout]

        s1 = s2 = "XXX"
        for fn in atomzip_fns:
            # gunzip file
            p = subprocess.run(["gunzip", fn]) 

            tmpfn = fn[:-3]
            firstline = True
            badpair = False
            with open(tmpfn) as f:
                anc1, anc2 = [], []
                iline = 0
                for line in f:
                    iline += 1
                    if not line.strip():
                        continue
                    if 'INIT' in line:
                        line = line[line.index("INIT"):]
                        fields = line.split()
                        if firstline:
                            distrfilename = locations['FSYSPATH']['distributions_data'] + 'distr_' + fields[1] + '.txt'
                            distrfile = open(distrfilename, 'w')
                            print("START WRITING", distrfilename)
                            firstline = False
                        else:
                            l = len(anc1)
                            distrfile.write('BEGIN CHAIN1: {0} CHAIN2: {1}\n'.format(s1, s2))
                            for j in range(l):
                                if anc1[j][0] != False and anc2[j][0] != False:
                                    dist = ((anc1[j][1] - anc2[j][1]) ** 2 + (anc1[j][2] - anc2[j][2]) ** 2 + (
                                                anc1[j][3] - anc2[j][3]) ** 2) ** 0.5
                                    distrfile.write("{0:10d} {1:20.6f}\n".format(anc1[j][0], dist))
                        anc1, anc2 = [], []
                        anc = anc1
                        s1 = fields[1]
                        s2 = fields[2]
                    elif line.strip().startswith('TER'):
                        if len(anc2) == 0:
                            anc = anc2
                        continue
                    elif line.strip().startswith('ATOM'):
                        fields = line.split()
                        nres = int(fields[4])
                        try:
                            cx, cy, cz = float(line[30:38].strip()), float(line[38:46].strip()), float(
                                line[46:54].strip())
                            anc.append((nres, cx, cy, cz))
                        except:
                            print('ERROR: the pair', s1, s2, 'CONTAINS AN ILL-FORMATTED LINE:\n', line,
                                  'THIS PAIR WILL NOT BE TAKEN INTO CONSIDERATION')
                            anc.append((False, False, False, False))

                distrfile.write('BEGIN CHAIN1: {0} CHAIN2: {1}\n'.format(s1, s2))
                l = len(anc1)
                for j in range(l):
                    if anc1[j][0] != False and anc2[j][0] != False:
                        dist = ((anc1[j][1] - anc2[j][1]) ** 2 + (anc1[j][2] - anc2[j][2]) ** 2 + (
                                    anc1[j][3] - anc2[j][3]) ** 2) ** 0.5
                        distrfile.write("{0:10d} {1:20.6f}\n".format(anc1[j][0], dist))
                distrfile.write('\nEND\n\n')
                distrfile.close()
            p = subprocess.run(["gzip", tmpfn])


def make_DSSP(dssp_path, struct_path, dssp_out_filename, fnull, opm_data, struct, opm_check=True,
              leave_original_dssp=False):
    if opm_check:
        pdbname = struct[:4]
        chain = struct[5]
        segments = opm_data[pdbname][chain]['segments']

    # Run DSSP
    print(dssp_path, '-i', struct_path, '-o', dssp_out_filename)
    p = subprocess.Popen([dssp_path, '-i', struct_path, '-o', dssp_out_filename], stdout=fnull, stderr=fnull)
    p.wait()

    # If no output file is created, return empty
    if not os.path.exists(dssp_out_filename):
        print("Struct", dssp_out_filename, "jumped")
        return []

    # Parse DSSP file and create a handy dictionary
    # Create a DSSP_segment list with elements such as ('A', 10, 34), where 'A' = helix and 'B' = sheet. DSSP: B, E = beta, H, G, I = alpha.
    dssp_out_file = open(dssp_out_filename, 'r')
    text = dssp_out_file.read().split('\n')
    dssp_out_file.close()

    dssp_to_edo = {'H': 'A', 'G': 'A', 'I': 'A', 'B': 'B', 'E': 'B'}
    residues = []
    dssp_segments = []
    ss_set = set()
    ss_dict = {}
    start_next = False
    for line in text:
        if not line:
            continue
        fields = line.split()
        if start_next and fields[1] != '!':
            if line[16] in list(dssp_to_edo.keys()):
                ss_type = dssp_to_edo[line[16]]
                ss_set.add(line[11] + '_' + line[5:10].strip())
                if leave_original_dssp:
                    ss_dict[line[11] + '_' + line[5:10].strip()] = line[16]
                else:
                    ss_dict[line[11] + '_' + line[5:10].strip()] = ss_type
        if fields[0].strip() == '#':
            start_next = True

    # This is the only (optional) modification from the original DSSP file
    if opm_check and not leave_original_dssp:
        print("Merge with OPM")
        for seg in opm_data[struct[:4]][struct[5]]['segments']:
            for nres in range(int(seg[0]), int(seg[1]) + 1):
                cres = struct[5] + '_' + str(nres)
                ss_set.add(cres)
                if cres not in list(ss_dict.keys()):
                    ss_dict[cres] = 'A'
                elif ss_dict[cres] != 'A':
                    print("Inconsistency!", cres, ss_dict[cres], 'A')

    new_dssp_segments = []
    dssp_seg = []
    for nres in opm_data[struct[:4]]['FROM_PDB'][struct[5]]['RESIDS']:
        nres = int(nres)
        cres = struct[5] + '_' + str(nres)
        if cres in ss_set:
            if not dssp_seg:
                dssp_seg.append(ss_dict[cres])
                dssp_seg.append(nres)
            elif dssp_seg[0] == ss_dict[cres]:
                dssp_seg.append(nres)
            elif dssp_seg[0] != ss_dict[cres]:
                new_dssp_segments.append(dssp_seg)
                dssp_seg = []
                dssp_seg.append(ss_dict[cres])
                dssp_seg.append(nres)
        else:
            if dssp_seg:
                new_dssp_segments.append(dssp_seg)
                dssp_seg = []
    return new_dssp_segments


def draw_topology(options, locations, neighbors, names, opm_data, an_locations):
    def create_cylinder(init_p, end_p, radius, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        text += "draw cylinder {{{0} {1} {2}}} {{{3} {4} {5}}} radius {6} resolution 100\n".format(init_p[0], init_p[1],
                                                                                                   init_p[2], end_p[0],
                                                                                                   end_p[1], end_p[2],
                                                                                                   radius)
        return text

    def create_slab(init_p, end_p, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        text += "draw cylinder {{{0} {1} {2}}} {{{3} {4} {5}}} radius {6}\n".format(init_p[0], init_p[1], init_p[2],
                                                                                    end_p[0], end_p[1], end_p[2], 1)
        return text

    def create_sphere(center_p, radius, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        text += "draw sphere {{{0} {1} {2}}} radius {3}\n".format(center_p[0], center_p[1], center_p[2], radius)
        return text

    def create_Bezier(BezierPoints, nsegments, radius, color='', gappy=False):
        print("BezierPoints", BezierPoints)
        centers = []
        text = ''
        N = len(BezierPoints) - 1
        for t in np.arange(0, 1 + 1 / nsegments, 1 / nsegments):
            BezierFormula = np.zeros_like(BezierPoints[0])
            for n in range(N + 1):
                BezierFormula += ((1 - t) ** (N - n)) * (t ** n) * ncr(N, n) * BezierPoints[n]
            centers.append(BezierFormula)
        if color:
            text = "draw color {0}\n".format(color)
        for nc in range(len(centers) - 1):
            text += create_sphere(centers[nc], radius, color='')
            if not gappy:
                text += create_cylinder(centers[nc], centers[nc + 1], radius, color='')
        text += create_sphere(centers[len(centers) - 1], radius, color='')
        return text

    def create_triangle(p1, p2, p3, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        text += "draw triangle {{{0} {1} {2}}} {{{3} {4} {5}}} {{{6} {7} {8}}}\n".format(p1[0], p1[1], p1[2], p2[0],
                                                                                         p2[1], p2[2], p3[0], p3[1],
                                                                                         p3[2])
        return text

    def create_rectangle(p1, p2, p3, p4, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        text += create_triangle(p1, p2, p3, color='')
        text += create_triangle(p3, p4, p1, color='')
        return text

    def create_convex_polygon(points, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        if len(points) < 3:
            raise NameError("For creating a polygon you have to pass at least 3 points")
        text += create_triangle(points[0], points[1], points[2], color='')
        for i in range(3, len(points)):
            text += create_triangle(points[i - 1], points[i], points[0], color='')
        return text

    def create_prism(points, hvector, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        text += create_convex_polygon(points, color='')
        text += create_convex_polygon(points + hvector, color='')
        for i in range(len(points) - 1):
            rectpoints = np.array([points[i], points[i + 1], points[i + 1] + hvector, points[i] + hvector])
            text += create_convex_polygon(rectpoints, color='')
        rectpoints = np.array([points[-1], points[0], points[0] + hvector, points[-1] + hvector])
        text += create_convex_polygon(rectpoints, color='')
        return text

    def create_arrow(init_p, end_p, normal, color=''):
        text = ''
        if color:
            text = "draw color {0}\n".format(color)
        vec = end_p - init_p
        normvec = vec / np.linalg.norm(vec)
        normal = np.array([0, 1, 0])
        widthvec = np.cross(normvec, normal)
        p1 = init_p + widthvec
        p2 = end_p - 4 * normvec + widthvec
        p3 = end_p - 4 * normvec - widthvec
        p4 = init_p - widthvec
        points = np.array([p1, p2, p3, p4])
        text += create_prism(points - normal / 4, normal / 2, color='')
        p1 = end_p - 4 * normvec + widthvec * 2
        p2 = end_p
        p3 = end_p - 4 * normvec - widthvec * 2
        points = np.array([p1, p2, p3])
        text += create_prism(points - normal / 4, normal / 2, color='')
        return text

    dssp_path = options['dssp_path']
    fnull = open(os.devnull, 'w')
    for struct in sorted(list(neighbors.keys())):
        print('struct', struct)
        struct_path = locations['FSYSPATH']['chains'] + struct + '.pdb'
        dssp_out_filename = an_locations['FSYSPATH']['DSSP'] + 'dssp_' + struct + '.txt'
        top_figname = an_locations['FSYSPATH']['topologies_figs'] + 'top_' + struct + '.jpeg'

        dssp_segments = make_DSSP(dssp_path, struct_path, dssp_out_filename, fnull, opm_data, struct)

        print(dssp_segments)

        for n_dssp_seg in range(len(dssp_segments)):
            if dssp_segments[n_dssp_seg][0][0] == 'A' and len(
                    dssp_segments[n_dssp_seg]) < 7:  # No helices with less than 6 CA
                dssp_segments[n_dssp_seg][0][0] == 'X'
            if dssp_segments[n_dssp_seg][0][0] == 'B' and len(
                    dssp_segments[n_dssp_seg]) < 4:  # No sheets with less than 3 CA in a row
                dssp_segments[n_dssp_seg][0][0] == 'X'

        # Read atom positions
        struct_file = open(struct_path, 'r')
        text = struct_file.read().split('\n')
        struct_file.close()
        bb_at = ['N', 'CA', 'C']
        residue = {}
        residue_list = []
        residue_dict = {}
        for line in text:
            if not line:
                continue
            fields = line.split()
            if not (fields[0] == 'ATOM' and line[13:16].strip() in bb_at):
                continue
            if not residue:
                residue['index'] = int(line[22:26].strip())
                residue_dict[int(line[22:26].strip())] = {}
                residue_dict[int(line[22:26].strip())][line[13:16].strip()] = np.array(
                    [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
                prev_index = int(line[22:26].strip())
            elif residue and int(line[22:26].strip()) != prev_index:
                residue_list.append(residue)
                residue_dict[int(line[22:26].strip())] = {}
                residue = {'index': int(line[22:26].strip())}
                prev_index = int(line[22:26].strip())
            residue[line[13:16].strip()] = np.array(
                [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
            residue_dict[int(line[22:26].strip())][line[13:16].strip()] = np.array(
                [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])

        exclude = set()
        for resid in list(residue_dict.keys()):
            if len(list(residue_dict[resid].keys())) != 3:
                exclude.add(resid)
        for resid in exclude:
            del residue_dict[resid]

        # Retrieve limits
        pdbname_path = locations['FSYSPATH']['whole'] + pdbname + '_opm.pdb'
        pdbname_file = open(pdbname_path, 'r')
        text = pdbname_file.read().split('\n')
        pdbname_file.close()
        lim_sup = ''
        lim_inf = ''
        for line in text:
            if not line or ('ATOM' not in line and 'HETATM' not in line):
                continue
            fields = line.split()
            if fields[3] == 'DUM' and fields[2] == 'O' and not lim_sup:
                lim_sup = float(line[46:54].strip())
            elif fields[3] == 'DUM' and fields[2] == 'N' and not lim_inf:
                lim_inf = float(line[46:54].strip())
            if lim_sup and lim_inf:
                break
        print("lim_sup", lim_sup, "lim_inf", lim_inf)

        new_dssp_segments = []
        for dssp_seg in dssp_segments:
            if not (dssp_seg[0][0] == 'A' or dssp_seg[0][0] == 'B'):
                continue
            if len(dssp_seg) < 3:  # If there only is one residue...
                continue
            init_point = ''
            end_point = ''
            init_mem = ''
            end_mem = ''
            if dssp_seg[0][0] == 'A':
                tail = 3
            elif dssp_seg[0][0] == 'B':
                tail = 1
            if len(dssp_seg) - tail < 3:
                continue
            for i in range(1, len(dssp_seg) - tail):
                if dssp_seg[0][0] == 'A':
                    gcenter = (residue_dict[dssp_seg[i]]['CA'] + residue_dict[dssp_seg[i + 1]]['CA'] +
                               residue_dict[dssp_seg[i + 2]]['CA'] + residue_dict[dssp_seg[i + 3]]['CA']) / 4
                elif dssp_seg[0][0] == 'B':
                    gcenter = (residue_dict[dssp_seg[i]]['N'] + residue_dict[dssp_seg[i]]['CA'] +
                               residue_dict[dssp_seg[i]]['C'] + residue_dict[dssp_seg[i + 1]]['N'] +
                               residue_dict[dssp_seg[i + 1]]['CA'] + residue_dict[dssp_seg[i + 1]]['C']) / 6

                if type(init_point) == str and init_point == '':
                    init_point = np.copy(gcenter)
                    continue

            end_point = np.copy(gcenter)
            print("init, end", dssp_seg[1], init_point, dssp_seg[-1], end_point)
            dist_vec = end_point - init_point
            norm_dist_vec = dist_vec / np.linalg.norm(dist_vec)
            if dssp_seg[0][0] == 'A':
                add_dist_init = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[1]]['CA'] - init_point)))
                add_dist_end = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[-1]]['CA'] - end_point)))
            elif dssp_seg[0][0] == 'B':
                add_dist_init = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[1]]['N'] - init_point)))
                add_dist_end = abs(np.dot(norm_dist_vec, (residue_dict[dssp_seg[-1]]['C'] - end_point)))
            init_point = init_point - add_dist_init * norm_dist_vec
            end_point = end_point + add_dist_end * norm_dist_vec

            if not ((init_point[2] > lim_sup and end_point[2] > lim_sup) or (
                    init_point[2] < lim_inf and end_point[2] < lim_inf)):
                if np.dot(end_point - init_point, np.array([0, 0, 1])) == 0:
                    init_mem = init_point
                    end_mem = end_point
                else:
                    line_param_sup = (lim_sup - init_point[2]) / norm_dist_vec[2]
                    point_sup = init_point + norm_dist_vec * line_param_sup
                    line_param_inf = (lim_inf - init_point[2]) / norm_dist_vec[2]
                    point_inf = init_point + norm_dist_vec * line_param_inf
                    if np.dot(end_point - init_point, point_inf - point_sup) > 0:
                        init_mem = point_sup
                        end_mem = point_inf
                    else:
                        init_mem = point_inf
                        end_mem = point_sup
                    i_im = init_mem - init_point
                    im_em = end_mem - init_mem
                    em_e = end_point - end_mem
                    if np.dot(i_im, im_em) < 0:
                        init_mem = init_point
                    if np.dot(em_e, im_em) < 0:
                        end_mem = end_point
            else:
                init_mem = ''
                end_mem = ''
            dssp_seg[0] = [dssp_seg[0][0], init_point, end_point, init_mem, end_mem]
            new_dssp_segments.append(dssp_seg)

        dssp_segments[:] = new_dssp_segments[:]

        # Perform rotaions! Mind to only rotate atoms C, CA, and N. Place first pivot at (0, 0, z). Translate all the complex. Then place second pivot at (x, 0, z) rotating all the complex. From the third on, place the pivot at (x, 0, z) rotating only residues after the last residue of the preceding loop / the corresponding ss structure (depending if pivot is initial or final). Stop at the (n-1)-th pivot.
        pivots = []
        for dssp_seg in dssp_segments:
            print(dssp_seg)
            if type(dssp_seg[0][3]) != str and type(dssp_seg[0][4]) != str:
                if not pivots:
                    pivots.append([dssp_seg[0][1], sorted(list(residue_dict.keys()))])
                    pivots.append([dssp_seg[0][2], [x for x in sorted(list(residue_dict.keys())) if x > dssp_seg[-1]]])
                else:
                    pivots.append([dssp_seg[0][1], [x for x in sorted(list(residue_dict.keys())) if x >= dssp_seg[1]]])
                    pivots.append([dssp_seg[0][2], [x for x in sorted(list(residue_dict.keys())) if x > dssp_seg[-1]]])
        if not pivots:
            pivots.append(
                [residue_dict[[x for x in sorted(list(residue_dict.keys()))][0]], sorted(list(residue_dict.keys()))])
            pivots.append([residue_dict[[x for x in sorted(list(residue_dict.keys()))][-1]], []])
        print("pivots", pivots)

        for resid in sorted(list(residue_dict.keys())):
            print(resid, list(residue_dict[resid].keys()))
            for atname in bb_at:
                residue_dict[resid][atname] = residue_dict[resid][atname] - pivots[0][0]
        trasl = pivots[0][0]
        for npiv in range(len(pivots)):
            pivots[npiv][0] = pivots[npiv][0] - trasl
        for dssp_seg in dssp_segments:
            for i in range(1, 5):
                if type(dssp_seg[0][i]) != str:
                    dssp_seg[0][i] = dssp_seg[0][i] - trasl

        for npiv in range(len(pivots) - 1):
            vec = pivots[npiv + 1][0] - pivots[npiv][0]
            vec[2] = 0.0
            ax = np.array([1, 0, 0])
            vec = vec / (np.linalg.norm(vec))
            print("norm vec:", np.linalg.norm(vec))
            R = np.zeros((3, 3))
            R[0][0] = vec[0] * ax[0] + vec[1] * ax[1]
            R[0][1] = -(vec[0] * ax[1] - ax[0] * vec[1])
            R[1][0] = -R[0][1]
            R[1][1] = R[0][0]
            R[2][2] = 1.0
            print("axis:", ax, "vector:", vec, "matrix:", R)
            for resid in pivots[npiv][1]:
                for atname in bb_at:
                    residue_dict[resid][atname] = np.dot(R, (residue_dict[resid][atname] - np.array(
                        [pivots[npiv][0][0], pivots[npiv][0][1], 0]))) + np.array(
                        [pivots[npiv][0][0], pivots[npiv][0][1], 0])
            # Rotate all pivots after this
            for nnpiv in range(npiv + 1, len(pivots)):
                pivots[nnpiv][0] = np.dot(R, (
                            pivots[nnpiv][0] - np.array([pivots[npiv][0][0], pivots[npiv][0][1], 0]))) + np.array(
                    [pivots[npiv][0][0], pivots[npiv][0][1], 0])
            # Rotate all endpoints after this
            for dssp_seg in dssp_segments:
                if dssp_seg[1] >= pivots[npiv][1][0]:
                    for i in range(1, 5):
                        if type(dssp_seg[0][i]) != str:
                            dssp_seg[0][i] = np.dot(R, (dssp_seg[0][i] - np.array(
                                [pivots[npiv][0][0], pivots[npiv][0][1], 0]))) + np.array(
                                [pivots[npiv][0][0], pivots[npiv][0][1], 0])

        # Compressing representation
        if neighbors[struct][0] == 'alpha':
            for npiv in range(2, len(pivots) - 1, 2):
                vec = pivots[npiv + 1][0] - pivots[npiv][0]
                # The compression is just a backtranslation of all pivots after the one considered. The backtranslation is along the x axis and amounts to the double sine of the tilt angle of that segment
                backtrans = np.array([-1, 0, 0]) * 2 * abs(vec[0])
                print("backtrans", backtrans, np.linalg.norm(backtrans))
                # We have to check that no pivot falls in the semi-plain before the one defined by the previous segment and the y-axis
                vecprev = pivots[npiv - 1][0] - pivots[npiv - 2][0]
                plainnorm = np.cross(vecprev, np.array([0, 1, 0])) / np.linalg.norm(
                    np.cross(vecprev, np.array([0, 1, 0])))
                do_compression = True
                # Check if moved helix does not clash
                checkv_n = np.linalg.norm(pivots[npiv + 1][0] + backtrans - pivots[npiv - 1][0])
                checkv = (pivots[npiv + 1][0] + backtrans - pivots[npiv - 1][0]) / checkv_n
                if not (np.dot(checkv, plainnorm) * np.dot(plainnorm, np.array([1, 0, 0])) > 0 and checkv_n * np.dot(
                        checkv, plainnorm) > 5):
                    print("No compression for", npiv, "due to", npiv + 1, checkv_n, checkv, plainnorm,
                          np.dot(checkv, [1, 0, 0]), np.dot(checkv, plainnorm), checkv_n * np.dot(checkv, plainnorm))
                    do_compression = False
                # Check if other helices do not clash with moved helix
                veccompr = vec + backtrans
                plainnorm2 = np.cross(vecprev, np.array([0, 1, 0])) / np.linalg.norm(
                    np.cross(vecprev, np.array([0, 1, 0])))
                for nnpiv in range(npiv + 2, min(len(pivots), npiv + 4)):
                    checkv2_n = np.linalg.norm(pivots[nnpiv][0] + backtrans - pivots[npiv][0])
                    checkv2 = (pivots[nnpiv][0] + backtrans - pivots[npiv][0]) / checkv2_n
                    if not (np.dot(checkv2, plainnorm2) * np.dot(plainnorm2,
                                                                 np.array([1, 0, 0])) > 0 and checkv_n * np.dot(checkv2,
                                                                                                                plainnorm2) > 6):
                        print("No compression for", npiv, "due to", nnpiv, checkv2_n, checkv2, plainnorm2,
                              np.dot(checkv2, [1, 0, 0]), np.dot(checkv2, plainnorm2),
                              checkv2_n * np.dot(checkv2, plainnorm2))
                        do_compression = False
                if do_compression:
                    print("Compression!", npiv)
                    for resid in pivots[npiv + 1][1]:
                        for atname in bb_at:
                            residue_dict[resid][atname] = residue_dict[resid][atname] + backtrans
                    for nnpiv in range(npiv + 1, len(pivots)):
                        pivots[nnpiv][0] = pivots[nnpiv][0] + backtrans
                    for dssp_seg in dssp_segments:
                        if dssp_seg[1] == pivots[npiv][1][0]:
                            for i in [2, 4]:
                                if type(dssp_seg[0][i]) != str:
                                    dssp_seg[0][i] = dssp_seg[0][i] + backtrans
                        elif dssp_seg[1] > pivots[npiv][1][0]:
                            for i in range(1, 5):
                                if type(dssp_seg[0][i]) != str:
                                    dssp_seg[0][i] = dssp_seg[0][i] + backtrans

        # Further compression of the x axis
        if neighbors[struct][0] == 'alpha':
            cfactor = 2
        elif neighbors[struct][0] == 'beta':
            cfactor = 2
        for resid in list(residue_dict.keys()):
            for atname in bb_at:
                residue_dict[resid][atname][0] = residue_dict[resid][atname][0] / cfactor
        for nnpiv in range(npiv + 1, len(pivots)):
            pivots[nnpiv][0][0] = pivots[nnpiv][0][0] / cfactor
        for dssp_seg in dssp_segments:
            for i in range(1, 5):
                if type(dssp_seg[0][i]) != str:
                    dssp_seg[0][i][0] = dssp_seg[0][i][0] / cfactor

        # For each DSSP_segment, write the command to create a cylinder (A) or a slab (B) with those endpoints.
        text = ''
        text += 'color Display Background white\n'
        text += 'rotate x by -90\n'
        big = 1.5
        small = 0.2
        BezierCA = {}
        for dssp_seg in dssp_segments:
            if dssp_seg[0][0] == 'A':
                text += create_cylinder(dssp_seg[0][1], dssp_seg[0][2], big, 'orange')
            elif dssp_seg[0][0] == 'B':
                normal = np.cross(residue_dict[dssp_seg[2]]['N'] - residue_dict[dssp_seg[1]]['CA'],
                                  residue_dict[dssp_seg[2]]['CA'] - residue_dict[dssp_seg[1]]['N'])
                normal = normal / np.linalg.norm(normal)
                text += create_arrow(dssp_seg[0][1], dssp_seg[0][2], normal, 'blue2')
            vec_dist = dssp_seg[0][2] - dssp_seg[0][1]
            norm_vec_dist = vec_dist / (np.linalg.norm(vec_dist))
            init_handle = dssp_seg[0][1] - norm_vec_dist
            end_handle = dssp_seg[0][2] + norm_vec_dist
            BezierCA[dssp_seg[1]] = {-1: init_handle, 0: dssp_seg[0][1], 1: ''}
            BezierCA[dssp_seg[-1]] = {-1: '', 0: dssp_seg[0][2], 1: end_handle}

        # for each DSSP_segment loop (not counting first and last) take the endpoint (if any) of the preceding SS piece. This is the first Bezier point. Then, take the normal vector of the preceding SS piece. The point at 1 A from the first point in that direction is the first Bezier handle. Second and third handles are the C and N atoms. Second point is the CA. Fourth handle is the point 1 A away of the CA in the opposite direction of the vector from N to C of that residue. Then, build the Bezier curve. Take parameter at (0, 0.05, 0.1, ..., 1.0) and in each site put a sphere, and between each two sites put a cylinder.

        init = True
        inside_a_ss = False
        for resid in sorted(list(residue_dict.keys())):
            if init:
                if resid in list(BezierCA.keys()):
                    inside_a_ss = True
                    continue
                else:
                    BezierCA[resid] = {-1: '', 0: residue_dict[resid]['CA'], 1: end_handle}
                init = False
            if not inside_a_ss:
                if resid in list(BezierCA.keys()):
                    if BezierCA[resid][1] == '':
                        inside_a_ss = True
                        continue
                res_N = residue_dict[resid]['N']
                res_C = residue_dict[resid]['C']
                init_handle = residue_dict[resid]['CA'] - (res_C - res_N) / (np.linalg.norm(res_C - res_N))
                end_handle = residue_dict[resid]['CA'] + (res_C - res_N) / (np.linalg.norm(res_C - res_N))
                BezierCA[resid] = {-1: init_handle, 0: residue_dict[resid]['CA'], 1: end_handle}
            if inside_a_ss and resid in list(BezierCA.keys()) and (BezierCA[resid][-1] == ''):
                inside_a_ss = False

        resids_in_BezierCA = sorted(list(BezierCA.keys()))
        init = False
        end = False
        gap = False
        count = 0
        BezierPoints = []
        for nresid in range(len(resids_in_BezierCA)):
            resid = resids_in_BezierCA[nresid]
            if BezierCA[resid][-1] == '':
                BezierPoints = []
                init = True
            if BezierCA[resid][1] == '' or count == 3 or (
                    nresid + 1 < len(resids_in_BezierCA) and type(BezierCA[resid][-1]) != str and np.linalg.norm(
                    BezierCA[resids_in_BezierCA[nresid + 1]][0] - BezierCA[resid][0]) > 4.5):
                end = True
            if type(BezierCA[resid][-1]) != str and nresid - 1 > 0 and np.linalg.norm(
                    BezierCA[resid][0] - BezierCA[resids_in_BezierCA[nresid - 1]][0]) > 4.5:
                end = True
                gap = True
            if end:
                BezierPoints.append(BezierCA[resid][-1])
            BezierPoints.append(BezierCA[resid][0])
            count += 1
            if init:
                BezierPoints.append(BezierCA[resid][1])
                init = False
            if end:
                print("BezierPoints", BezierPoints)
                if gap:
                    text += create_Bezier(BezierPoints, 10, small, 'gray', gappy=True)
                else:
                    text += create_Bezier(BezierPoints, 100, small, 'gray')
                init = False
                end = False
                gap = False
                count = 1
                BezierPoints = []
                BezierPoints.append(BezierCA[resid][0])
                BezierPoints.append(BezierCA[resid][1])

        limright = 0
        limleft = 0
        for resid in sorted(list(residue_dict.keys())):
            if limright < residue_dict[resid]['CA'][0]:
                limright = residue_dict[resid]['CA'][0]
            if limleft > residue_dict[resid]['CA'][0]:
                limleft = residue_dict[resid]['CA'][0]

        limup = lim_sup - trasl[2]
        limdown = lim_inf - trasl[2]
        if int(opm_data[struct[:4]][struct[5]]['ntm']) < 13:
            text = 'display resize {0} {1}\n'.format(8000, 6000) + text
        else:
            text = 'display resize {0} {1}\n'.format(16000, 6000) + text
        text += 'display projection Orthographic\n'
        text += 'display depthcue off\n'

        p1 = np.array([limright, 50, limup])
        p2 = np.array([limright, 50, limdown])
        p3 = np.array([limleft, 50, limdown])
        p4 = np.array([limleft, 50, limup])
        text += create_rectangle(p1, p2, p3, p4, color='silver')

        scale = 0.01
        text += 'scale to {0}\n'.format(scale)
        text += 'translate to {0} {1} {2}\n'.format(-(limright - limleft) / 2 * scale, 0,
                                                    -(limup + limdown) / 2 * scale)
        text += 'axes location off\n'

        tgafig = an_locations['FSYSPATH']['topologies_data'] + 'tmp_' + struct + '.tga'
        text += 'render TachyonInternal "{0}"\n'.format(tgafig)
        text += 'quit\n'

        vmdscript_filename = an_locations['FSYSPATH']['topologies_data'] + struct + '.tcl'
        vmdscript_file = open(vmdscript_filename, 'w')
        vmdscript_file.write(text)
        vmdscript_file.close()

        p = subprocess.Popen([options['vmd_path'], "-dispdev", "text", "-e", vmdscript_filename], stdout=fnull)
        p.wait()

        im = Image.open("{0}".format(tgafig))
        pix = np.asarray(im)

        pix = pix[:, :, 0:3]  # Drop the alpha channel
        idx = np.where(pix - 255)[0:2]  # Drop the color when finding edges
        box = list(map(min, idx))[::-1] + list(map(max, idx))[::-1]

        region = im.crop(box)
        region_pix = np.asarray(region)

        fig, ax = plt.subplots()
        img = ax.imshow(region_pix)
        ax.axis('off')
        fig.savefig(top_figname, bbox_inches='tight', dpi=1000)
        fig.savefig(top_figname[:-4] + '.eps', bbox_inches='tight', dpi=1000)
        os.remove(tgafig)


def analysis_archive(locations, analysis_path, generate=True, analysis_name='analysis/'):
    """
	Creates a completely separate and self-contained file system to be shown to the users
	"""

    if (not analysis_path) and (analysis_path != ''):
        raise NameError("ERROR: analysis directory not specified. Add -an <path> or --analysisdir <path> flag.")

    an_locations = {'MAIN': analysis_path, 'FSYS': collections.OrderedDict(), 'FSYSPATH': collections.OrderedDict(),
                    'SYSFILES': collections.OrderedDict(), 'OPT': collections.OrderedDict(),
                    'FIXED': collections.OrderedDict()}
    an_locations['FSYS']['analysis'] = analysis_name  # analysis/
    an_locations['FSYS']['ancache'] = 'analysis_cache/'
    an_locations['FSYS']['cache'] = an_locations['FSYS']['ancache'] + 'cache/'
    an_locations['FSYS']['images'] = an_locations['FSYS']['analysis'] + 'images/'  # analysis/images/
    an_locations['FSYS']['images_figs'] = an_locations['FSYS'][
                                              'images'] + 'images_figs/'  # analysis/images/images_figs/
    an_locations['FSYS']['images_pml'] = an_locations['FSYS'][
                                             'images'] + 'images_pymolscripts/'  # analysis/images/images_pymolscripts/
    an_locations['FSYS']['images_figs_whole'] = an_locations['FSYS'][
                                                    'images_figs'] + 'whole_structs/'  # analysis/images/images_figs/whole_structs/
    an_locations['FSYS']['images_figs_chain'] = an_locations['FSYS'][
                                                    'images_figs'] + 'chains/'  # analysis/images/images_figs/chains/
    an_locations['FSYS']['images_pml_whole'] = an_locations['FSYS'][
                                                   'images_pml'] + 'whole_structs/'  # analysis/images/images_pymolscripts/whole_structs/
    an_locations['FSYS']['images_pml_chain'] = an_locations['FSYS'][
                                                   'images_pml'] + 'chains/'  # analysis/images/images_pymolscripts/chains/
    an_locations['FSYS']['estimators'] = an_locations['FSYS']['analysis'] + 'estimators/'  # analysis/estimators/
    an_locations['FSYS']['nlists'] = an_locations['FSYS']['analysis'] + 'neighbor_lists/'  # analysis/neighbor_lists/
    an_locations['FSYS']['seqneighs'] = an_locations['FSYS'][
                                            'nlists'] + 'seq_neighbors/'  # analysis/neighbor_lists/seq_neighbors/
    an_locations['FSYS']['strneighs'] = an_locations['FSYS'][
                                            'nlists'] + 'str_neighbors/'  # analysis/neighbor_lists/str_neighbors/
    an_locations['FSYS']['totneighs'] = an_locations['FSYS'][
                                            'nlists'] + 'tot_neighbors/'  # analysis/neighbor_lists/tot_neighbors/
    an_locations['FSYS']['distributions'] = an_locations['FSYS'][
                                                'analysis'] + 'distributions/'  # analysis/distributions/
    an_locations['FSYS']['distributions_figs'] = an_locations['FSYS'][
                                                     'distributions'] + 'distributions_figs/'  # analysis/distributions/distributions_figs/
    an_locations['FSYS']['distributions_data'] = an_locations['FSYS'][
                                                     'distributions'] + 'distributions_data/'  # analysis/distributions/distributions_data/
    an_locations['FSYS']['densityscatter'] = an_locations['FSYS'][
                                                 'analysis'] + 'densityscatter/'  # analysis/densityscatter/
    an_locations['FSYS']['densityscatter_figs'] = an_locations['FSYS'][
                                                      'densityscatter'] + 'densityscatter_figs/'  # analysis/densityscatter/densityscatter_figs/
    an_locations['FSYS']['densityscatter_maps'] = an_locations['FSYS'][
                                                      'densityscatter'] + 'densityscatter_maps/'  # analysis/densityscatter/densityscatter_maps/
    an_locations['FSYS']['polar'] = an_locations['FSYS']['analysis'] + 'polarplots/'  # analysis/polarplots/
    an_locations['FSYS']['polar_figs'] = an_locations['FSYS'][
                                             'polar'] + 'polarplots_figs/'  # analysis/polarplots/polarplots_figs/
    an_locations['FSYS']['polar_data'] = an_locations['FSYS'][
                                             'polar'] + 'polarplots_data/'  # analysis/polarplots/polarplots_data/
    an_locations['FSYS']['polar_maps'] = an_locations['FSYS'][
                                             'polar'] + 'polarplots_maps/'  # analysis/polarplots/polarplots_maps/
    an_locations['FSYS']['DSSP'] = an_locations['FSYS']['analysis'] + 'DSSP/'  # analysis/DSSP/
    an_locations['FSYS']['topologies'] = an_locations['FSYS']['analysis'] + 'topologies/'  # analysis/topologies/
    an_locations['FSYS']['topologies_figs'] = an_locations['FSYS'][
                                                  'topologies'] + 'topologies_figs/'  # analysis/topologies/topologies_figs/
    an_locations['FSYS']['topologies_data'] = an_locations['FSYS'][
                                                  'topologies'] + 'topologies_data/'  # analysis/topologies/topologies_data/
    an_locations['FSYS']['databasewide'] = an_locations['FSYS'][
                                               'analysis'] + 'databasewide_graphs/'  # analysis/databasewide_graphs/
    # FSYSPATH
    for name, val in an_locations['FSYS'].items():
        an_locations['FSYSPATH'][name] = analysis_path + val
    # SYSFILES
    an_locations['SYSFILES']['neighborstable'] = an_locations['FSYSPATH']['analysis'] + 'neighbors_table.txt'
    an_locations['SYSFILES']['neighborstlist'] = an_locations['FSYSPATH']['analysis'] + 'neighbors_tot_list.txt'
    an_locations['SYSFILES']['neighborsselist'] = an_locations['FSYSPATH']['analysis'] + 'neighbors_seq_list.txt'
    an_locations['SYSFILES']['neighborsstlist'] = an_locations['FSYSPATH']['analysis'] + 'neighbors_str_list.txt'
    an_locations['SYSFILES']['strdata'] = an_locations['FSYSPATH']['cache'] + 'str_data_straln.pkl'
    an_locations['SYSFILES']['strdatatmpl'] = an_locations['FSYSPATH']['cache'] + 'str_data_template.json'
    an_locations['SYSFILES']['background'] = an_locations['FSYSPATH']['databasewide'] + 'densitylines.png'

    # Generate filesystem
    if generate:
        if not os.path.exists(analysis_path):
            os.mkdir(analysis_path)
        for index, duple in enumerate(an_locations['FSYSPATH'].items()):
            if not os.path.exists(duple[1]):
                os.mkdir(duple[1])
        print('Analysis filesystem was generated')

        shutil.copyfile(locations['SYSFILES']['data_structure_template'], an_locations['SYSFILES']['strdatatmpl'])
        shutil.copyfile(locations['FSYSPATH']['cache'] + os.path.basename(an_locations['SYSFILES']['strdata']),
                        an_locations['SYSFILES']['strdata'])

    return an_locations


def graphical_analysis(options, locations, str_data, analysisdir='None'):
    if analysisdir == 'None':
        analysisdir = locations['FSYSPATH']['analysis']
    print("INIT ARCHIVE")
    an_locations = analysis_archive(locations, analysisdir, generate=True)  # True for Amar
    print("NEIGHBOR LISTS")
    pdbi_chs = neighborlists(options, locations['SYSFILES']['summarytable'], an_locations, no_print=False)
    do_only = None  # pdbi_chs[:20]
    print("FINAL_N_OF_CHAINS", len(pdbi_chs))
    ntms = extract_ntms(locations['SYSFILES']['structurewise_table'])
    scatterplots(options, locations, pdbi_chs, ntms, an_locations, do_eps=True)
    # return
    print("POLAR")
    polar(options, locations, pdbi_chs, ntms, an_locations,
          do_only=do_only)  # , only_write=False, do_only=rlist, do_eps=True)	# True for Amar
    straln_profiles(locations, do_only, str_data, an_locations, do_eps=True)  # , do_only=rlist, do_eps=True)
    return
    draw_topology(options, locations, neighbors, names, opm_data, an_locations)


def parallel_graphical_analysis(options, locations, str_data, do_only, ntms, sp_pkl_fn, analysisdir='None'):
    if analysisdir == 'None':
        analysisdir = locations['FSYSPATH']['analysis']

    print("INIT ARCHIVE", analysisdir)
    an_locations = analysis_archive(locations, analysisdir, generate=False)  # True for Amar
    print(an_locations)
    scatterplots_parallel(an_locations, do_only, ntms, sp_pkl_fn, do_eps=False)
    all_entries = []
    range_list = [x for x in range(ord('a'), ord('z') + 1)]
    range_list += [x for x in range(ord('0'), ord('9') + 1)]
    for ic1 in range_list:
        for ic2 in range_list:
            all_entries += [os.path.basename(x[0][-10:-4]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
                "ls {0}/strneigh_?{1}{2}*".format(an_locations['FSYSPATH']['strneighs'], chr(ic1), chr(ic2)),
                shell=True, stdout=subprocess.PIPE).stdout] if x]
    print(all_entries)
    polar(options, an_locations, all_entries, ntms, an_locations, do_only=do_only)
    straln_profiles(an_locations, do_only, str_data, an_locations, do_eps=True)


def parallel_densityscatter(options, locations, str_data, do_only, ntms, sp_pkl_fn, analysisdir='None'):
    if analysisdir == 'None':
        analysisdir = locations['FSYSPATH']['analysis']

    print("INIT ARCHIVE", analysisdir)
    an_locations = analysis_archive(locations, analysisdir, generate=False)  # True for Amar
    print(an_locations)
    scatterplots_parallel(an_locations, do_only, ntms, sp_pkl_fn, do_eps=False)


def parallel_polar(options, locations, str_data, do_only, ntms, sp_pkl_fn, analysisdir='None'):
    if analysisdir == 'None':
        analysisdir = locations['FSYSPATH']['analysis']

    print("INIT ARCHIVE", analysisdir)
    an_locations = analysis_archive(locations, analysisdir, generate=False)  # True for Amar
    print(an_locations)
    all_entries = []
    range_list = [x for x in range(ord('a'), ord('z') + 1)]
    range_list += [x for x in range(ord('0'), ord('9') + 1)]
    for ic1 in range_list:
        for ic2 in range_list:
            all_entries += [os.path.basename(x[0][-10:-4]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
                "ls {0}/strneigh_?{1}{2}*".format(an_locations['FSYSPATH']['strneighs'], chr(ic1), chr(ic2)),
                shell=True, stdout=subprocess.PIPE).stdout] if x]
    print(all_entries)
    polar(options, locations, all_entries, ntms, an_locations, do_only=do_only)


def parallel_profile(options, locations, str_data, do_only, ntms, sp_pkl_fn, analysisdir='None'):
    if analysisdir == 'None':
        analysisdir = locations['FSYSPATH']['analysis']

    print("INIT ARCHIVE", analysisdir)
    an_locations = analysis_archive(locations, analysisdir, generate=False)  # True for Amar
    print(an_locations)
    all_entries = []
    range_list = [x for x in range(ord('a'), ord('z') + 1)]
    range_list += [x for x in range(ord('0'), ord('9') + 1)]
    for ic1 in range_list:
        for ic2 in range_list:
            all_entries += [os.path.basename(x[0][-10:-4]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
                "ls {0}/strneigh_?{1}{2}*".format(an_locations['FSYSPATH']['strneighs'], chr(ic1), chr(ic2)),
                shell=True, stdout=subprocess.PIPE).stdout] if x]
    do_only = all_entries
    straln_profiles(locations, do_only, str_data, an_locations, do_eps=True)


def parallel_polprof(options, locations, str_data, do_only, ntms, sp_pkl_fn, analysisdir='None'):
    if analysisdir == 'None':
        analysisdir = locations['FSYSPATH']['analysis']

    print("INIT ARCHIVE", analysisdir)
    an_locations = analysis_archive(locations, analysisdir, generate=False)  # True for Amar
    all_entries = []
    range_list = [x for x in range(ord('a'), ord('z') + 1)]
    range_list += [x for x in range(ord('0'), ord('9') + 1)]
    for ic1 in range_list:
        for ic2 in range_list:
            all_entries += [os.path.basename(x[0][-10:-4]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
                "ls {0}/strneigh_?{1}{2}*".format(an_locations['FSYSPATH']['strneighs'], chr(ic1), chr(ic2)),
                shell=True, stdout=subprocess.PIPE).stdout] if x]
    print(all_entries)
    polar(options, locations, all_entries, ntms, an_locations, do_only=do_only)
    straln_profiles(locations, do_only, str_data, an_locations, do_eps=True)


def dsc_test(analysisdir):
    print("INIT")
    options, locations = initialize_repository()
    print("ANLOC")
    an_locations = analysis_archive(locations, analysisdir, generate=False)
    print("NEIGHLIST")
    pdbi_chs = neighborlists(options, locations['SYSFILES']['summarytable'], an_locations, no_print=False)
    print("NTMS")
    ntms = extract_ntms(locations['SYSFILES']['structurewise_table'])
    print("SCPL BKG")
    pkl_fn = an_locations['FSYS']['analysis'] + "sp.pkl"
    scatterplots_background(an_locations, pdbi_chs, pkl_fn=an_locations['MAIN'] + pkl_fn)
    print("###### SCPL ######")

    do_onlys = [pdbi_chs[x:x + 10] for x in range(0, len(pdbi_chs), 100)]
    do_onlys = do_onlys[:5]
    do_only = sum(do_onlys, [])
    # do_only = pdbi_chs[:10]
    scatterplots_parallel(an_locations, do_only, ntms, an_locations['MAIN'] + pkl_fn, do_eps=False)
    # exit(1)


def extract_ntms(fn):
    xinf = {}
    first = True
    # print(fn)
    with open(fn) as f:
        for line in f:
            if first:
                first = False
                continue
            fields = line.split('\t')
            # print(fields)
            if fields[0] not in xinf:
                xinf[fields[0]] = {}
            if (len(fields) == 14):  # contains uniprot code
                xinf[fields[0]][fields[1]] = int(fields[7])
                # xinf[fields[1]][fields[2]] = int(fields[7])
            elif (len(fields) == 13):  # no uniprot code column
                print(fields[0], "has", len(fields), "data columns. ntms = ", fields[6])
                xinf[fields[0]][fields[1]] = int(fields[6])
    return xinf


def test_graphical_analysis(options, locations, str_data, analysisdir):
    if analysisdir == 'None':
        analysisdir = locations['FSYSPATH']['analysis']

    print("INIT ARCHIVE")
    an_locations = analysis_archive(locations, analysisdir,
                                    generate=True)  # True for Amar
    print("NEIGHBOR LISTS")
    pdbi_chs = neighborlists(options, locations['SYSFILES']['summarytable'], an_locations, no_print=False)
    print("FINAL_N_OF_CHAINS", len(pdbi_chs))
    ntms = extract_ntms(locations['SYSFILES']['structurewise_table'])

    pkl_fn = an_locations['FSYS']['analysis'] + "sp.pkl"

    scatterplots_background(an_locations, pdbi_chs, pkl_fn=an_locations['MAIN'] + pkl_fn)
    straln_profiles_preprocess(options, locations, an_locations)

    do_only = pdbi_chs[:10]

    scatterplots_parallel(an_locations, do_only, ntms, sp_pkl_fn, do_eps=False)
    all_entries = []
    for ic1 in range(ord('a'), ord('z') + 1):
        for ic2 in range(ord('a'), ord('z') + 1):
            all_entries += [os.path.basename(x[0][-10:-4]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
                "ls {0}/strneigh_?{1}{2}*".format(an_locations['FSYSPATH']['strneighs'], chr(ic1), chr(ic2)),
                shell=True, stdout=subprocess.PIPE).stdout] if x]
    print(all_entries)
    polar(options, locations, all_entries, ntms, an_locations, do_only=do_only)
    straln_profiles(locations, do_only, str_data, an_locations, do_eps=True)

def update_graphs(options, filters, locations, opm_data, table, more_list, less_list, analysisdir='None'):
    if analysisdir == 'None':
        analysisdir = options['analysisdir']
    an_locations = analysis_archive(locations, analysisdir, generate=False)
    neighborlists_creator(options, locations, opm_data, table, an_locations, only_dicts=False)  # True for Amar
    modified_list = list(set(more_list) | set(less_list))
    scatterplots(options, locations, table, an_locations)
    return

    polar(options, locations, table, neighbors, an_locations, only_write=False, do_only=redo_graphs)  # True for Amar
    straln_profiles(locations, neighbors, opm_data, table, an_locations, do_only=redo_graphs)


### Checks
def sc_bkg_UT(localmain):
    options, locations = initialize_repository()
    an_locations = analysis_archive(locations, localmain, generate=False)
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')
    straln_profiles_preprocess(options, locations, an_locations, batch_job_code="check_" + timestamp + "_FrTMAlign000")
    exit(1)


def ds_UT(localmain):
    print("\n\n\n#### initialize_repository\n\n")
    options, locations = initialize_repository()
    print("\n\n\n#### analysis_archive\n\n")
    an_locations = analysis_archive(locations, localmain, generate=False)
    print("\n\n\n#### pdbi_chs\n\n")
    pdbi_chs = neighborlists(options, locations['SYSFILES']['summarytable'], an_locations, no_print=True)
    pickle.dump(pdbi_chs, open(locations['FSYSPATH']['cache'] + "pdbi_chs.pkl", 'wb'))
    print("\n\n\n#### background\n\n")
    # scatterplots_background(an_locations, pdbi_chs, pkl_fn=an_locations['MAIN']+pkl_fn)
    # Background already present
    print("\n\n\n#### do_only\n\n")
    do_onlys = [pdbi_chs[x:x + 10] for x in range(0, len(pdbi_chs), 100)]
    do_onlys = do_onlys[:5]
    do_only = sum(do_onlys, [])
    # do_only = pdbi_chs[:10]
    print("\n\n\n#### ntms\n\n")
    ntms = pickle.load(open(an_locations['FSYSPATH']['analysis'] + "ntms.pkl", 'rb'))
    print("\n\n\n#### pkl\n\n")
    pkl_fn = an_locations['FSYSPATH']['analysis'] + "sp.pkl"
    print("\n\n\n#### scatterplots_parallel\n\n")
    print(do_only)
    scatterplots_parallel(an_locations, do_only, ntms, pkl_fn, do_eps=False)
    print("\n\n\n#### DONE")


def profiles_UT(localmain):
    print("\n\n\n#### initialize_repository\n\n")
    options, locations = initialize_repository()
    print("\n\n\n#### analysis_archive\n\n")
    an_locations = analysis_archive(locations, localmain, generate=False)
    print("\n\n\n#### pdbi_chs\n\n")
    pdbi_chs = pickle.load(open(locations['FSYSPATH']['cache'] + "pdbi_chs.pkl", 'rb'))
    print("\n\n\n#### do_only\n\n")
    do_onlys = [pdbi_chs[x:x + 10] for x in range(0, len(pdbi_chs), 100)]
    do_onlys = do_onlys[:5]
    do_only = sum(do_onlys, [])
    # do_only = pdbi_chs[:10]
    str_data_pkl_fn = locations['FSYSPATH']['cache'] + 'str_data_completegen.pkl'
    str_data_templ_fn = locations['FSYSPATH']['cache'] + 'str_data_template.json'
    str_data = read_checkpoint(str_data_pkl_fn, str_data_templ_fn)
    print("\n\n\n#### straln_profiles\n\n")
    straln_profiles(locations, do_only, str_data, an_locations, do_eps=True)
    print("\n\n\n#### DONE")


def missing_entries_gen():
    options, locations = initialize_repository()
    pdbi_chs = pickle.load(open(locations['FSYSPATH']['analyses'] + "pdbi_chs.pkl", 'rb'))
    do_only = []
    with open(options['PATHS'][('me', 'missing_entries')]) as f:
        for line in f:
            if line.strip():
                do_only.append(line.strip())
    str_data_pkl_fn = locations['FSYSPATH']['cache'] + 'str_data_straln.pkl'
    str_data_templ_fn = locations['FSYSPATH']['cache'] + 'str_data_template.json'
    str_data = read_checkpoint(str_data_pkl_fn, str_data_templ_fn)

    ntms = extract_ntms(locations['SYSFILES']['structurewise_table'])
    pkl_fn = locations['FSYSPATH']['analyses'] + "sp.pkl"

    scatterplots_background(locations, pdbi_chs, pkl_fn=pkl_fn)
    scatterplots_parallel(locations, do_only, ntms, pkl_fn, do_pdf=True)
    all_entries = []
    ords = list(range(ord('a'), ord('z') + 1)) + list(range(ord('0'), ord('9') + 1))
    for ic1 in ords:
        for ic2 in ords:
            all_entries += [os.path.basename(x[0][-10:-4]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
                "ls {0}/strneigh_?{1}{2}*".format(locations['FSYSPATH']['strneighs'], chr(ic1), chr(ic2)), shell=True,
                stdout=subprocess.PIPE).stdout] if x]
    all_entries = sorted(all_entries)
    polar(options, locations, all_entries, ntms, do_only=do_only, do_pdf=True)
    straln_profiles(locations, do_only, str_data, locations, do_pdf=True)
    print("Finished creating plots for: ", do_only)

if __name__ == "__main__":
    from supporting_functions import *
    from initialize_repository import *
    from combine_sources import *

    kw = ['PSCATTER', 'PPOLAR', 'PPROFILE', 'PARALLEL', 'PPOLPROF']
    if sys.argv[-1] in kw:  # PARALLEL MUST BE PUT AS LAST ARGUMENT
        print("Initiate parallel")
        str_data_pkl_fn = sys.argv[1]
        str_data_templ_fn = sys.argv[2]
        analysisdir = sys.argv[3]
        sp_pkl_fn = sys.argv[4]
        ntms_fn = sys.argv[5]
        do_only_fn = sys.argv[6]

        str_data = read_checkpoint(str_data_pkl_fn, str_data_templ_fn)
        do_only = []
        with open(do_only_fn) as df:
            for line in df:
                do_only.append(line.strip())

        with open(ntms_fn, 'rb') as f:
            ntms = pickle.load(f)
        print("Call subroutine")
        if sys.argv[-1] == 'PARALLEL':
            parallel_graphical_analysis({}, {}, str_data, do_only, ntms, sp_pkl_fn, analysisdir=analysisdir)
        elif sys.argv[-1] == 'PSCATTER':
            parallel_densityscatter({}, {}, str_data, do_only, ntms, sp_pkl_fn, analysisdir=analysisdir)
        elif sys.argv[-1] == 'PPOLAR':
            parallel_polar({}, {}, str_data, do_only, ntms, sp_pkl_fn, analysisdir=analysisdir)
        elif sys.argv[-1] == 'PPROFILE':
            parallel_profile({}, {}, str_data, do_only, ntms, sp_pkl_fn, analysisdir=analysisdir)
        elif sys.argv[-1] == 'PPOLPROF':
            parallel_polprof({}, {}, str_data, do_only, ntms, sp_pkl_fn, analysisdir=analysisdir)
        exit(1)
    else:
        timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')
        print("########### MAIN - 00 - TIMESTAMP", timestamp)
        locusts_tmp = 'analysis_parallel_locusts_cache_{0}/'.format(timestamp)

        print("\n\n########### MAIN - 01 - INITIALIZE")
        options, locations = initialize_repository()
        str_data_pkl_fn = locations['FSYSPATH']['cache'] + 'str_data_straln.pkl'
        str_data_templ_fn = locations['FSYSPATH']['cache'] + 'str_data_template.json'
        str_data = read_checkpoint(str_data_pkl_fn, str_data_templ_fn)

        print("\n\n########### MAIN - 02 - NEIGHBOR LISTS")
        if options['ALL']['no_neighbor_lists_analysis']:
            print("MAIN - 02a - LOAD NEIGHBORS")
            pdbi_chs = pickle.load(open(locations['FSYSPATH']['analyses'] + "pdbi_chs.pkl", "rb"))
            numtotneighs = pickle.load(open(locations['FSYSPATH']['analyses'] + "totn.pkl", "rb"))
        else:
            print("MAIN - 02b - FIND NEIGHBORS")
            pdbi_chs, numtotneighs = neighborlists(options, locations['SYSFILES']['summarytable'], locations, no_print=False)
            pickle.dump(pdbi_chs, open(locations['FSYSPATH']['analyses'] + "pdbi_chs.pkl", "wb"))
            pickle.dump(numtotneighs, open(locations['FSYSPATH']['analyses'] + "totn.pkl", "wb"))

        print("MAIN - 02 - PDBI CHS", pdbi_chs)
        print("MAIN - 02 - NUMTOTNEIGHS", numtotneighs)

        print("\n\n########### MAIN - 03 - BATCHES")
        BATCH_SIZE = options['ALL']['batch_size_analysis']
        print("MAIN - 03 - BATCH SIZE", BATCH_SIZE)

        if options['ALL']['filter_entries_analysis']:
            print("MAIN - 03a - FILTER ENTRIES")

            kkww = [('densityscatter_figs', 'ds_XXX.png'), ('polar_figs', 'p_XXX.png'), ('distributions_figs', 'distr_XXX.png')]

            chk = []
            if sys.argv[-1][1:] == 'PARALLEL':
                chk == kkww
            else:
                chk = [kkww[kw.index(sys.argv[-1][1:])]]

            proto_do_onlys = set()
            for pdbi_ch in pdbi_chs:
                for d, fn in chk:
                    chk_fn = locations['FSYSPATH'][d] + '/' + fn.replace('XXX', pdbi_ch)
                    if not os.path.exists(chk_fn):
                        print('MISSING', chk_fn)
                        proto_do_onlys.add(pdbi_ch)

            pdo = sorted(list(proto_do_onlys))
            print(pdo)
            print(len(pdo))

            print("MAIN - 03a - BATCHES")
            do_onlys = [pdo[x:x+BATCH_SIZE] for x in range(0, len(pdo), BATCH_SIZE)]
        else:
            print("MAIN - 03b - ALL ENTRIES")
            # Collect "sizes" to sort
            sizes = []
            if options['ALL']['sort_by_chain_size_analysis']:
                print("MAIN - 03b.i - SORT BY CHAIN SIZE")
                # The size of each chain of each complex is extracted from the main data structure
                for x in pdbi_chs:
                    pdbi, ch = x[:4], x[5]
                    chi = str_data[pdbi]['ENCOMPASS']['structure']['kchains'].index(ch)
                    s = len(str_data[pdbi]['ENCOMPASS']['structure']['chains'][chi]['sequence'])
            elif options['ALL']['sort_by_num_neighbors_analysis']:
                print("MAIN - 03b.i - SORT BY NUMBER OF NEIGHBORS")
                for x in pdbi_chs:
                    sizes.append((x, numtotneighs[x]))
            else:
                print("No sorting method selected")
                exit(1)
    
            # The chains are sorted by their size
            sizes = sorted(sizes, key=lambda x: x[1])
            print("SIZES", len(sizes))
            print(sizes)
            # The chains are divided in batches of at most 100 elements, where chains are distributed
            #  so that their sizes would be balanced (i.e. no batch contains only very short or very long chains)
            nbatches = 1 + len(sizes) // BATCH_SIZE
            do_onlys = []
            doc = []  # for testing
            for i in range(nbatches):
                do_onlys.append([])
                doc.append([])
            for ipp, pp in enumerate(sizes):
                x, s = pp
                do_onlys[ipp % nbatches].append(x)
                doc[ipp % nbatches].append((x, s))  # for testing
            for i in range(len(doc)):
                print(f"BATCH {i} {[(j[0], j[1]) for j in doc[i]]}")
                # print("BATCH", i, " ".join([str(j[1]) for j in doc[i]]))

        print("MAIN - 03 - ENTRIES", do_onlys)
        print(len(do_onlys))

        if options['ALL']['quick_test_analysis']:
            print("\n\n########### MAIN - 04 - TEST RUN")
            print("\n\n\n\n\n\n\nTHIS IS A TEST\n\n\n\n\n\n\n")
            do_onlys = do_onlys[:5]

        if not options['ALL']['no_structurewisetable_analysis']:
            print("\n\n########### MAIN - 05 - MAKE STRUCTUREWISE TABLE")
            make_structurewise_table(str_data, locations['SYSFILES']['structurewise_table'])

        print("\n\n########### MAIN - 06 - RUN LOCUSTS")
        fn_opt = run_locusts_on_analysis.default_fn_options(options)
        if not os.path.exists(locusts_tmp):
            os.mkdir(locusts_tmp)

        command = 'source {0}; python3 {1} <arg0> <arg1> <arg2> <arg3> <arg4> <arg5> '.format(
            options['ALL']['hpcactipy'],
            'EncoMPASS/' + os.path.basename(os.path.abspath(__file__))
            ) + sys.argv[-1][1:]

        # Figure parameters
        figsize = (8, 8)
        points_per_inch = 72.272  # Conversion factor
        pkl_fn = locations['FSYSPATH']['analyses'] + "sp.pkl"
        ntms_fn = locations['FSYSPATH']['analyses'] + "ntms.pkl"

        with open(ntms_fn, 'wb') as f:
            pickle.dump(extract_ntms(locations['SYSFILES']['structurewise_table']), f)

        print("SCATTERPLOTS BACKGROUND")
        scatterplots_background(locations, pdbi_chs, pkl_fn=pkl_fn)  # Creates pkl_fn!!

        print("straln_profiles_preprocess")
        straln_profiles_preprocess(options, locations)

        hpc_test_name = 'analysis_parallel_{0}'.format(timestamp)
        stepname = 'analysis_parallel'

        # copy pkl in analyses
        str_pkl = locations['FSYSPATH']['analyses'] + 'str_data.pkl'
        templ = locations['FSYSPATH']['analyses'] + 'str_data_template.json'
        shutil.copyfile(locations['FSYSPATH']['cache'] + 'str_data_straln.pkl', str_pkl)
        shutil.copyfile(locations['SYSFILES']['data_structure_template'], templ)

        workdir = locations['FSYSPATH']['ancache']  # must be rewritten in run_locusts starting from env_root
        # These arguments are relative to the workdir!!!
        fixed_options = [
            str_pkl,
            templ,
            locations['FSYSPATH']['analysis'],
            pkl_fn,
            ntms_fn
        ]

        run_locusts_on_analysis.run_locusts(
            fn_opt,
            locations,
            workdir,
            hpc_test_name,
            stepname,
            locusts_tmp,
            do_onlys,
            command,
            fixed_options
        )

    print(f"Warning: {locusts_tmp} is still there and occupies a lot of space! You chould think about removing it as soon as you check everything is in place.")
