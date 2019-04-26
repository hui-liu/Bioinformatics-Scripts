import sys
import os
import re
from collections import OrderedDict
import numpy as np
import pandas as pd
from ete3 import Tree
import fastcluster
from collections import OrderedDict


def phylogenetic_tree_to_cluster_format(tree, pairwise_estimates):
    id_map = {
        pairwise_estimates.index[i]: i for i in range(len(pairwise_estimates))}
    t = Tree(tree)
    # midpoint rooting
    midpoint = t.get_midpoint_outgroup()
    if not midpoint:  # midpoint = None when their are only two leaves
        midpoint = list(t.get_leaves())[0]
    t.set_outgroup(midpoint)

    # algorithm for getting cluster data structure
    n = len(id_map)
    out = []
    pairwise_distances = {}
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node.name = id_map[node.name]
            id_map[node.name] = node.name  # add identity map for renamed nodes
            # to id_map for line below
            pairwise_distances[node.name] = {
                id_map[x.name]: node.get_distance(x) for x in t.get_leaves()
            }
        else:
            node.name = n
            n += 1
            children = node.get_children()
            out.append(
                [children[0].name, children[1].name,
                 children[0].get_distance(children[1]),
                 len(node.get_leaves())])
    return np.array(out), pairwise_distances

def average_linkage_clustering(pairwise_estimates):
    clustering = fastcluster.average(pairwise_estimates)
    return clustering

def compute_weights(df, min_ks=0.005, max_ks=5):
    df = df[~df.index.duplicated()]  # for safety
    df["WeightOutliersIncluded"] = round(1 / df.groupby(['Family', 'Node'])[
        'Ks'].transform('count'), 5)
    df_ = df[df["Ks"] <= max_ks]
    df_ = df_[df_["Ks"] >= min_ks]
    df["WeightOutliersExcluded"] = np.zeros(len(df.index))
    df.loc[df_.index, "WeightOutliersExcluded"] = round(1 / df_.groupby(
            ['Family', 'Node'])['Ks'].transform('count'), 5)
    return df

def _weighting(pairwise_estimates):
    if pairwise_estimates is None:
        return None, None, None

    if pairwise_estimates['Ks'].shape[0] < 2:
        return None, None, None

    pairwise_distances = None

    # Average linkage clustering based on Ks
    clustering = average_linkage_clustering(pairwise_estimates['Ks'])
    return clustering, pairwise_distances

def _calculate_weighted_ks(clustering, pairwise_estimates, Family,
                           pairwise_distances=None):
    # None -> None
    if pairwise_estimates is None or clustering is None:
        return None

    # process the clustering structure to get weights
    leaves = pairwise_estimates['Ks'].shape[0]
    nodes = {i: [i] for i in range(leaves)}

    weights = {}
    out = set()
    for x in range(clustering.shape[0]):
        node_1, node_2, distance = clustering[x, 0], clustering[x, 1], \
                                   clustering[x, 2] * 2
        grouping_node = leaves + x
        nodes[grouping_node] = nodes[node_1] + nodes[node_2]
        for i in nodes[node_1]:
            for j in nodes[node_2]:
                if pairwise_distances:
                    distance = pairwise_distances[i][j]
                p1 = pairwise_estimates['Ks'].index[i]
                p2 = pairwise_estimates['Ks'].index[j]
                pair = "__".join(sorted([p1, p2]))
                weights[pair] = [
                    p1, p2, Family,
                    pairwise_estimates['Ka'].iloc[i, j],
                    pairwise_estimates['Ks'].iloc[i, j],
                    pairwise_estimates['Omega'].iloc[i, j],
                    pairwise_estimates['d4'].iloc[i, j],
                    grouping_node
                ]

                if pairwise_estimates['Ks'].iloc[i, j] > 5:
                    out.add(grouping_node)

    df = pd.DataFrame.from_dict(weights, orient='index')
    df.columns = ['Paralog1', 'Paralog2', 'Family',
                  'Ka', 'Ks', 'Omega', 'd4', 'Node']

    res = compute_weights(df)
    return res

def clusterlst(filename):
    clusters = OrderedDict()
    with open(filename, 'r') as f:
        for line in f:
            lsp = line.split()
            id, gr = lsp[0].split("|")[1], lsp[1]
            clusters.setdefault(gr, []).append(id)
    return clusters

def write_control(f, control_dict):
    for key in control_dict.keys():
        f.write('{0} = {1}\n'.format(key, control_dict[key]))

def parse_codeml_out(codeml_out):
    seqid = {}
    #paras = {}
    #count1 = 0
    #count2 = 0
    columns = set()
    with open(codeml_out, 'r') as f:
        fcont = f.read()
    for line in fcont.split("\n"):
        if "Reading seq #" in line:
            temp1 = line.split()
            seqid[temp1[-2].lstrip("#").rstrip(":")] = temp1[-1]
            columns.add(temp1[-1])
    results_dict = {
            'Ka': pd.DataFrame(
                np.zeros((len(list(columns)), len(list(columns)))),
                index=sorted(list(columns)),
                columns=sorted(list(columns))),
            'Ks': pd.DataFrame(
                np.zeros((len(list(columns)), len(list(columns)))),
                index=sorted(list(columns)),
                columns=sorted(list(columns))),
            'Omega': pd.DataFrame(
                np.zeros((len(list(columns)), len(list(columns)))),
                index=sorted(list(columns)),
                columns=sorted(list(columns))),
            'd4': pd.DataFrame(
                np.zeros((len(list(columns)), len(list(columns)))),
                index=sorted(list(columns)),
                columns=sorted(list(columns)))
        }
    for line in fcont.split("\n"):
        VS = re.findall('\s+(\d+)\s+vs.\s+(\d+)', line)
        if VS:
            #count1 += 1
            g1, g2 = [seqid[i] for i in VS[0]]
            #id = '__'.join(sorted([g1, g2]))
            #paras[count1] = [id]
        elif "four-fold sites" in line and 'd4' in line:
            #count2 += 1
            temp2 = [x for x in line.replace("=","").split()]
            w, dN, dS, d4 = temp2[1: -3: 2]
            results_dict['Ka'][g1][g2] = dN
            results_dict['Ka'][g2][g1] = dN
            results_dict['Ks'][g1][g2] = dS
            results_dict['Ks'][g2][g1] = dS
            results_dict['Omega'][g1][g2] = w
            results_dict['Omega'][g2][g1] = w
            results_dict['d4'][g1][g2] = d4
            results_dict['d4'][g2][g1] = d4
            #paras[count2] += [dN, dS, w, d4]
    return results_dict

def parseFasta(filename):
    fas = {}
    id = None
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                id = header.split()[0]
                fas[id] = []
            else:
                fas[id].append(line.rstrip().upper())
        for id, seq in fas.items():
            fas[id] = ''.join(seq)
    return fas

def pal2nal(pal, nuc_seqs):
    nal = {}
    nal_ = ''
    for pid, seq in pal.items():
        if pid not in nuc_seqs:
            print("{0} not in cds file".format(pid))
        else:
            nuc = nuc_seqs[pid]
            nal_ = ''
            j = 0
            for i in range(len(seq)):
                if seq[i] == '-':
                    nal_ += '---'
                else:
                    nal_ += nuc[j:j + 3]
                    j += 3
            nal[pid] = nal_
    return nal

def alndict2phy(alndict, cluster):
    num = str(len(alndict))
    l = str(len(list(alndict.values())[0]))
    outname = cluster +  ".phy"
    with open("TEMP/" + outname, 'w') as f:
        f.write(num + " " + l + "\n")
        for s in alndict:
            f.write(s + "    " + alndict[s] + "\n")
    return outname

def alignment(pep_dict, nuc_dict, cluster):
    fasfile = cluster + ".fasta"
    pepaln = cluster  + "." + "aln"
    os.system('mkdir TEMP')
    with open("TEMP/" + fasfile, "w") as f:
       for s in pep_dict:
           f.write(">" + s + "\n" + pep_dict[s] + "\n")
    os.system('/usr/local/bin/mafft --localpair --maxiterate 1000 --quiet --preservecase TEMP/%s > TEMP/%s' % (fasfile, pepaln))
    pepaln_dict = parseFasta("TEMP/" + pepaln)
    codon_aln = pal2nal(pepaln_dict, nuc_dict)
    phy = alndict2phy(codon_aln, cluster)
    return phy

def run_phyml(msafilename, phyml_path='phyml'):
    msa_phyml = "TEMP/" + msafilename
    tree_path = "TEMP/" + msafilename + '.nw'
    log = msa_phyml + '.log'
    command = " ".join([phyml_path, '-i', msa_phyml, '-q', '-d', 'aa', '>', log, '2>&1'])
    os.system(command)
    os.system(' '.join(['mv', '{}_phyml_tree*'.format(msa_phyml), tree_path]))
    return tree_path

def run_codeml(msafilename, cluster):
    control_file = cluster + '.ctrl'
    out_file = cluster + '.codeml'
    logfile = cluster + ".log"
    control = {
        'seqfile': msafilename, 'outfile': out_file,
        'noisy': 9, 'verbose': 0, 'runmode': -2, 'seqtype': 1,
        'CodonFreq': 2, 'clock': 0, 'aaDist': 0,
        'aaRatefile': 'dat/jones.dat', 'model': 0, 'NSsites': 0, 'icode': 0,
        'Mgene': 0, 'fix_kappa': 0,
        'kappa': 2, 'fix_omega': 0, 'omega': .4, 'fix_alpha': 1, 'alpha': 0,
        'Malpha': 0, 'ncatG': 8,
        'getSE': 0, 'RateAncestor': 1, 'Small_Diff': .5e-6, 'cleandata': 1,
        'method': 0
    }
    with open('TEMP/' + control_file, 'w') as f:
        write_control(f, control)
    os.chdir("TEMP")
    os.system('/home/liuhui/bin/paml4.9i/bin/codeml %s > %s' % (control_file, logfile))
    results_dict = parse_codeml_out(logfile)
    os.chdir("..")
    #os.system('rm 2ML.dN 2ML.dS 2ML.t 2NG.dN 2NG.dS 2NG.t rst rst1 rub')
    #os.system('rm %s.fasta %s.aln %s.phy %s.ctrl %s.codeml %s.log' % (id, id, id, id, id, id))
    #os.system('rm -r TEMP')
    return results_dict

clusters = clusterlst(sys.argv[1])
pep_dict = parseFasta(sys.argv[2])
nuc_dict = parseFasta(sys.argv[3])
out = open(sys.argv[4], 'w')
out.write("\t".join(['Pair', 'Gene1', 'Gene2', 'Family', 'Ka', 'Ks', 'Omega', 'd4', 'Node', 'WeightOutliersIncluded', 'WeightOutliersExcluded']) + "\n")
for c in clusters:
    p = clusters[c]
    pseqs_dict = {k:pep_dict[k] for k in p}
    nseqs_dict = {k:nuc_dict[k] for k in p}
    msa = alignment(pseqs_dict, nuc_dict, c)
    pairwise_estimates = run_codeml(msa, c)
    if len(p) > 2:
        tree_path = run_phyml(msa)
        clustering, pairwise_distances = phylogenetic_tree_to_cluster_format(
                tree_path, pairwise_estimates['Ks'])
        results =  _calculate_weighted_ks(clustering, pairwise_estimates, c)
        for i in results.index.values:
            lst = [i] + results.loc[i].tolist()
            out.write("\t".join(map(str, lst)) + "\n")
    elif len(p) == 2:
        clustering, pairwise_distances = _weighting(pairwise_estimates)
        results =  _calculate_weighted_ks(clustering, pairwise_estimates, c)
        for i in results.index.values:
            lst = [i] + results.loc[i].tolist()
            out.write("\t".join(map(str, lst)) + "\n")
    os.system('rm -r TEMP')
out.close()
