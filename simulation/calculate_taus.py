from simulate_traits import phylogenetic_covariance,simulate_traits,combinations,np,Path
from dendropy import Tree
from scipy.stats import kendalltau
from scipy.spatial.distance import squareform,pdist
from joblib import Parallel, delayed
from sys import argv


def patristic_distances(tree):
    """We must sort the ns lexicographically, because that's how the sigmas were generated"""
    # Kendall's Tau will treat ((A,B),C) as different from ((B,A),C), so will always overestimate dist.  We use the patristic distances to account for this. BUT we will always have ties due to these symmetries, so need to use Kendall's Tau-b coefficient.
    ns = sorted(tree.taxon_namespace)
    # ntaxa = len(ns)
    # dmat = np.zeros((ntaxa, ntaxa))

    pdm=tree.phylogenetic_distance_matrix()
    return np.array([pdm.distance(*pair) for pair in combinations(ns,2)])


def get_phylo_cov(t:str):
    return phylogenetic_covariance(Tree.get_from_string(t,'newick')) 

def patristic_distances2(tree):

    
    pdm=tree.phylogenetic_distance_matrix()
    return pdm.distances()

def calculate_dist(tree, traits):
    if isinstance(tree,str):
        tree=Tree.get_from_string(tree,'newick')
    tree_dist = patristic_distances(tree)
    taus=np.empty(len(traits))
    for i,t in enumerate(traits):
        d = pdist(t[np.newaxis].T)
        taus[i] = kendalltau(tree_dist,d).statistic
    return taus

def main():
    njobs=int(argv[1]) if len(argv)>1 else 12

    for h in [1, 3, 5, 7]:
        for i, d in enumerate(
            Path('/N/project/phyloML/rate_timescaling/data/tree_sims').rglob(f'*/*/h{h}')):
            # [Path('/N/project/phyloML/rate_timescaling/data/tree_sims/l11_m5.5_h5723175.894/')]):
            new_dir = d/'no_ils'
            new_dir.mkdir(exist_ok=True)
            if not (new_dir / 'samples_ngenes1.npz').exists() or (new_dir/'taus.npy').exists() :continue
            print('processing no ILS (BD) tree from ',d)
            try:
                with open(d/'parent_trees_4plus.phy','r') as f:
                    parent_trees = f.readlines()
            except FileNotFoundError:
                continue
            with Parallel(njobs) as parallel:
                try:
                    sigmas=np.load(new_dir/f'sigmas_species_tree.npz').values()
                except:
                    sigmas = parallel(delayed(get_phylo_cov)(t) for t in parent_trees)
                    np.savez(new_dir/f'sigmas_species_tree.npz', *sigmas)
                # dists = *map(patristic_distances, parent_trees),
                # break

                # this will load the samples from the npz file, or generate them if they don't exist
                samples = simulate_traits([1], [sigmas], data_dir=new_dir, overwrite=False) #parent_trees x #nsamples x ntaxa

                if len(samples)!=len(parent_trees):
                    continue # or regenerate

                # for t,s in zip(parent_trees,samples):
                #     calculate_dist(t,s)
                taus = np.array(parallel(delayed(calculate_dist)(t,s) for t,s in zip(parent_trees,samples)))
            np.save(new_dir/'taus.npy',taus)

if __name__=='__main__':
    main()
