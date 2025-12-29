from simulation.calculate_sigmas import *
from sys import argv
import warnings

# warnings.filterwarnings("error")
max_rep = 20


def process_dir_object(
    d: Path,
    calculate_pic,
    overwrite: bool = True,
    ils=False,
    eps=None,
    ngenes=1,
    max_replicates=None,
):
    tree_file = d / "parent_trees_4plus_pos_edges.phy"
    if ils:
        samp_file = d / f"samples_ngenes{ngenes}.npz"
        outdir = d
    else:
        samp_file = d / "no_ils/samples_ngenes1.npz"
        if ngenes != 1:
            raise ValueError("ngenes must be 1 unless ILS is True")
        outdir = d / "no_ils"
    outfile = Path(
        f"pic_js_sigmas{max_rep}_ngenes{ngenes}_full_v2"
        + (f"_{eps}" if eps else "")
        + ".npy"
    )

    if (
        (not overwrite and outfile.exists())
        or not samp_file.exists()
        or not tree_file.exists()
    ):
        return 1
    # print(samp_file)

    with open(tree_file, "r") as f:
        parent_trees = [
            dendropy.Tree.get_from_string(t, "newick", rooting="force-rooted")
            for t in f
        ]
    npzfile = np.load(samp_file)
    n = len(npzfile.files)
    samples = [v / ngenes for v in npzfile.values()]
    sigmas = np.empty((1, n, max_rep))
    for i, (t, s) in enumerate(zip(parent_trees, samples)):
        try:
            pic = PIC(
                t=t,
                chars=s,
                calculate_pic=calculate_pic,
                max_replicates=max_replicates,
                eps=eps,
            )
            #             we've verified that cherries/paired doesn't change under JS

            sigmas[:, i, :] = pic.estimate_sigmas()
        except ZeroDivisionError:
            print("error", i, args[0])

    outfile = (
        f"pic_js_sigmas{max_rep}_ngenes{ngenes}_full_v2"
        + (f"_{eps}" if eps else "")
        + ".npy"
    )
    np.save(outdir / outfile, sigmas[0])
    return 0


njobs = int(argv[1])
print(njobs, os.cpu_count(), parent_dir)


def process(dirname):
    return process_dir_object(
        dirname,
        overwrite=False,
        ngenes=ngenes,
        max_replicates=max_rep,
        ils=False,
        calculate_pic=JamesSteinPIC,
    )


# for h in (7,1,3,5):
#     for d in parent_dir.glob(f'*/*/h{h}') :
#         process(d)

# for d in parent_dir.glob(f"*/*/h7"):
#     process(d)

with Parallel(n_jobs=njobs) as parallel:
    skipped = parallel(
        delayed(process)(d) for h in (3, 7, 5, 1) for d in parent_dir.glob(f"*/*/h{h}")
    )
print("skipped:", sum(skipped))
