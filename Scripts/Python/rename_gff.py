#!/usr/bin/env python
# -*- coding: utf-8 -*-
# modified from https://www.jianshu.com/p/e2f13e2f42ba

import re
import sys
import argparse
import logging
import gffutils
from collections import defaultdict

def rename(args):
    seqid2name = dict()
    newseqid = []
    for line in open(args.change, 'r'):
        tem = line.strip().split()
        seqid2name[tem[0]] = tem
        newseqid.append(tem[0])

    db = gffutils.create_db(args.gff,':memory:',force=True, keep_order=True, merge_strategy="create_unique", sort_attribute_values = False)
    mRNA_children=("five_prime_UTR","three_prime_UTR","CDS","exon")
    idmap={
        "CDS" : "cds",
        "exon" : "exon",
        "five_prime_UTR":"5utrp",
        "three_prime_UTR":"3utrp"
    }

    lst = {}
    mRNA_lst = {}
    seqid = None
    for gene in db.features_of_type("gene",order_by=('seqid','start','end')):
        if gene.seqid != seqid:
            genenum = 0
        genenum += args.addnum
        seqid = gene.seqid

        genename = '{0}G{1:06}'.format(seqid2name[seqid][2],genenum)
        lst.setdefault(seqid, []).append('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={geneid}\n'.format(seqid=seqid2name[seqid][1], source=gene.source, featuretype=gene.featuretype, start=gene.start, end=gene.end, score=gene.score, strand=gene.strand, frame=gene.frame, geneid=genename))

        for t,mRNA in enumerate(db.children(gene,featuretype="mRNA",order_by=('seqid','start','end'))):
            mrna_num = t + 1
            mrnaid = '{genename}.{num}'.format(genename=genename,num=mrna_num)
            mRNA_lst.setdefault(seqid, []).append([mRNA.id, mrnaid])
            lst.setdefault(seqid, []).append('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={mrnaid};Parent={geneid}\n'.format(seqid=seqid2name[mRNA.seqid][1], source=mRNA.source, featuretype=mRNA.featuretype, start=mRNA.start, end=mRNA.end, score=mRNA.score, strand=mRNA.strand, frame=mRNA.frame, mrnaid=mrnaid, geneid=genename))

            numdict = defaultdict(int)
            for child in db.children(mRNA,featuretype=mRNA_children,order_by=("start",'end'), reverse=False):
                exonnum = len(list(db.children(mRNA, featuretype="exon")))
                cdsnum = len(list(db.children(mRNA, featuretype="CDS")))
                utr5num = len(list(db.children(mRNA, featuretype="five_prime_UTR")))
                utr3num = len(list(db.children(mRNA, featuretype="three_prime_UTR")))
                numdict[child.featuretype] += 1
                if gene.strand == "+":
                    child_id = '{genename}.{childid}{num}'.format(genename=mrnaid,childid=idmap[child.featuretype],num=numdict[child.featuretype])
                    lst.setdefault(seqid, []).append('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={child_id};Parent={mrnaid}\n'.format(seqid=seqid2name[child.seqid][1], source=child.source, featuretype=child.featuretype, start=child.start, end=child.end, score=child.score, strand=child.strand, frame=child.frame, child_id=child_id, mrnaid=mrnaid))
                else:
                    if child.featuretype == "exon":
                        child_id = '{genename}.{childid}{num}'.format(genename=mrnaid,childid=idmap[child.featuretype],num=exonnum-numdict[child.featuretype]+1)
                    elif child.featuretype == "CDS":
                        child_id = '{genename}.{childid}{num}'.format(genename=mrnaid,childid=idmap[child.featuretype],num=cdsnum-numdict[child.featuretype]+1)
                    elif child.featuretype == "five_prime_UTR":
                        child_id = '{genename}.{childid}{num}'.format(genename=mrnaid,childid=idmap[child.featuretype],num=utr5num-numdict[child.featuretype]+1)
                    elif child.featuretype == "three_prime_UTR":
                        child_id = '{genename}.{childid}{num}'.format(genename=mrnaid,childid=idmap[child.featuretype],num=utr3num-numdict[child.featuretype]+1)
                    lst.setdefault(seqid, []).append('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={child_id};Parent={mrnaid}\n'.format(seqid=seqid2name[child.seqid][1], source=child.source, featuretype=child.featuretype, start=child.start, end=child.end, score=child.score, strand=child.strand, frame=child.frame, child_id=child_id, mrnaid=mrnaid))

    f_out = open('%s.gene.gff3' % args.prefix, 'w')
    mRNA_out = open("gene_id_conversion.tsv", 'w')
    for i in newseqid:
        if i in lst:
            for j in lst[i]:
                f_out.write(j)
    f_out.close()
    for x in newseqid:
        if x in mRNA_lst:
            for y in mRNA_lst[x]:
                mRNA_out.write("\t".join(y) + "\n")
    mRNA_out.close()

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="rename gff3 file")

    parser.add_argument('-g', '--gff', required=True, help='gff3 file')
    parser.add_argument('-c', '--change', required=True, help='a file, correspondence between sequence name and gene name prefix')
    parser.add_argument('-a', '--addnum', type=int,default=1, help='diff in gene number, such as if addnum = 10, xx1G000010, xx1G000020')
    parser.add_argument('-p', '--prefix', default='result', help='prefix of output')
    args = parser.parse_args()

    rename(args)


if __name__ == "__main__":
    main()
