from simulation.calculate_sigmas import *
from sys import argv
import warnings

# warnings.filterwarnings("error")
max_rep = 20


def process_dir_object(
    d: Path,
    overwrite: bool = False,
    ils=False,
    eps: float =None,
    ngenes: int =1,
    calculate_pic=continuous.PhylogeneticIndependentContrasts,
    max_replicates=None,
):
    """Assumes ngenes=1 is equiv to no ILS"""
    tree_file = d / "parent_trees_4plus_pos_edges.phy"
    if ils:
        samp_file = d / f"samples_ngenes{ngenes}.npz"
        outdir = d
    else:
        samp_file = d / "no_ils/samples_ngenes1.npz"
        if ngenes != 1:
            raise ValueError("ngenes must be 1 unless ILS is True")
        outdir = d / "no_ils"
    outfile = outdir / (
        f"pic_dendropy_sigmas{max_rep}_ngenes{ngenes}_cherry{'_{eps:.2f}' if eps else ''}.npy"
    )

    if (
        (not overwrite and outfile.exists())
        or not samp_file.exists()
        or not tree_file.exists()
    ):
        return 1

    with open(tree_file, "r") as f:
        parent_trees = [
            dendropy.Tree.get_from_string(t, "newick", rooting="force-rooted")
            for t in f
        ]
    npzfile = np.load(samp_file)
    n = len(npzfile.files)
    samples = [v / ngenes for v in npzfile.values()]
    sigmas = np.empty((3, n, max_rep))
    for i, (t, s) in enumerate(zip(parent_trees, samples)):
        try:
            pic = PIC(
                t=t,
                chars=s,
                calculate_pic=calculate_pic,
                max_replicates=max_replicates,
                eps=eps,
            )
            res = np.vstack(
                (
                    pic.estimate_sigmas(),
                    pic.cherries_sigma(),
                    pic.paired_lineages_sigma(),
                )
            )
            sigmas[:, i, :] = res
        except ZeroDivisionError:
            print("error", i, d, args[0])
            print( pic.estimate_sigmas().shape, pic.cherries_sigmas().shape)
    for suffix, sig in zip(("full", "cherries", "paired"), sigmas):
        outfile = f"pic_dendropy_sigmas{max_rep}_ngenes{ngenes}_{suffix}" + (f"_{eps}" if eps else '') + ".npy"
        np.save(outdir / outfile, sig)
    print('wrote',outfile)
    return 0


def process(dirname, eps):
    return process_dir_object(
        dirname,
        overwrite=False,
        eps = eps,
        ngenes=ngenes,
        max_replicates=max_rep,
        ils=False,
        calculate_pic=continuous.PhylogeneticIndependentContrasts,
    )


# for h in (7,1,3,5):
#     for d in parent_dir.glob(f'*/*/h{h}') :
#         process(d)

# for d in parent_dir.glob(f"*/*/h7"):
#     process(d)


njobs = int(argv[1])
eps = float(argv[2]) if len(argv)>1 else None

print(njobs, os.cpu_count(), parent_dir)
import random

dirs = list(parent_dir.glob(f"*/*/h*"))
random.shuffle(dirs)

with Parallel(n_jobs=njobs, timeout = 9000) as parallel:
    skipped = parallel(
        delayed(process)(d, eps)  for d in dirs 
    )
print("skipped:", sum(skipped))
