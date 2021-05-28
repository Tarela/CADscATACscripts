import os
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


def main(args):
    maindir = args[2]
    if not maindir.endswith("/"):
        maindir += "/"
    if not os.path.isdir(maindir+'step3_interpretation/ism_scores/'+args[0]):
        os.mkdir(maindir+'step3_interpretation/ism_scores/'+args[0])
    setup_pool(args[0], int(args[1]), args[2], args[3])


def setup_pool(cluster, workers,maindir,gkmdir):
    basedir = maindir+'/step3_interpretation/'
    modeldir = maindir+'/step2_alliteration/%s/models/'%(cluster)
    test_pool = []
    for fold in range(10):
        effect_fasta = basedir + 'ism_inputs/' + cluster + '.effect.fasta'
        noneffect_fasta = basedir + 'ism_inputs/' + cluster + '.noneffect.fasta'
        #model = basedir + 'gkmsvm/Cluster' + cluster + '/fold' + str(fold) + '/train/train.output.model.txt'
        model = "%s%s.%s.model.txt"%(modeldir,cluster,fold)
        effect_output = basedir + 'ism_scores/' + cluster + '/fold' + str(fold) + '.effect.scores'
        noneffect_output = basedir + 'ism_scores/' + cluster + '/fold' + str(fold) + '.noneffect.scores'
        test_pool.append((effect_fasta, model, effect_output,gkmdir))
        test_pool.append((noneffect_fasta, model, noneffect_output,gkmdir))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        merge=pool.map(test_svm, test_pool)


def test_svm(inputs):
    os.system(inputs[3]+'/gkmpredict -T 16 ' + inputs[0] + ' ' + inputs[1] + ' ' + inputs[2])


if __name__ == "__main__":
    cluster = sys.argv[1]
    workers = sys.argv[2]
    maindir = sys.argv[3]
    gkmdir = sys.argv[4]
    main([cluster, workers,maindir,gkmdir ])
