import os
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


def main(args):
    maindir = args[2]
    if not maindir.endswith("/"):
        maindir += "/"
    if not os.path.isdir(maindir+'/step3_interpretation/delta_scores/'+args[0]):
        os.mkdir(maindir+'/step3_interpretation/delta_scores/'+args[0])
    setup_pool(args[0], int(args[1]), args[2])


def setup_pool(cluster, workers,maindir):
    basedir = maindir+'/step3_interpretation/'
    modeldir = maindir+'/step2_alliteration/%s/models/'%(cluster)
    delta_pool = []
    for fold in range(10):
        effect_input = basedir + 'ism_inputs/'+cluster+'.effect.fasta'
        noneffect_input = basedir + 'ism_inputs/'+cluster+'.noneffect.fasta'
        kmer_scores = basedir + 'kmer_scores/'+cluster+'/fold'+str(fold)+'.scores'
        kmer_output = basedir + 'delta_scores/' + cluster + '/fold' + str(fold) + '.delta.scores'
        delta_pool.append((noneffect_input, effect_input, kmer_scores, kmer_output,maindir))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        merge=pool.map(get_delta, delta_pool)


def get_delta(inputs):
    os.system('perl ' + inputs[4] + '/step3_interpretation/deltasvm.pl ' + inputs[0] + ' ' + inputs[1] + ' ' + inputs[2] + ' ' + inputs[3])


if __name__ == "__main__":
    main(sys.argv[1:])
