#!/usr/bin/env python


# Copyright (C) 2018  Shengwei Hou, housw2010@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import sys
import os
import argparse
from ete3 import NCBITaxa
from collections import Counter
from collections import OrderedDict


class Fasta(object):

    def __init__(self, header, seq):
        self.header = header
        self.seq = seq.upper()

    def get_gc_content(self):
        ret = 0
        if len(self.seq) == 0:
            return ret
        for c in self.seq:
            if c == "G" or c == "C":
                ret += 1
        return float(ret)/len(self.seq)

    def __str__(self):
        return ">"+self.header+"\n"+self.seq+"\n"


def parse_fasta(fasta_file):
    """
    :param fasta_file: input fasta file
    :return:           yield Fasta record as a generator
    """
    header = ""
    seq = []
    with open(fasta_file, "r") as ih:
        for line in ih:
            if line.startswith(">"):
                if header:
                    yield Fasta(header, "".join(seq))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        yield Fasta(header, "".join(seq))


def get_feature_tree_from_GLASSgo_sRNAs(GLAssgo_sRNA_file):
    """
    :param GLAssgo_sRNA_file: GLASSgo output sRNA file in fasta format
    :param required_taxid: []
    :return: an ete tree with features including ["sci_name", "taxid", "rank", "count"]
    """

    # get taxid2counts for all leaves
    taxids = []
    for sRNA in parse_fasta(GLAssgo_sRNA_file):
        header = sRNA.header
        if "taxID:" in header:
            taxid = header.strip().split("taxID:")[-1]
            # for taxid like this one: taxID:1639;1230340
            if ";" in taxid:
                taxid = taxid.split(";")[-1]
            taxids.append(taxid)
    taxid2counts = dict(Counter(taxids))

    # aggregate counts from leaves to parents
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(taxids, intermediate_nodes=True)
    print(tree.get_ascii(attributes=["sci_name"]))
    for leaf in tree.iter_leaves():
        leaf.add_features(count=taxid2counts.get(leaf.name, 0))
    for node in tree.traverse("postorder"):
        count = 0
        if node.is_leaf():
            continue
        else:
            children = node.get_children()
            for child in children:
                count += child.count
        node.add_features(count=count)

    return tree


def get_taxid2counts(ete3_feature_tree, rank=True):
    """
    :param ete3_feature_tree: the ete3 feature tree returned by get_feature_tree_from_GLASSgo_sRNAs
    :param rank: only count taxids whose rank is not 'no rank'
    :return: taxid2counts, taxid2annotation
    """

    taxid2counts = {}
    taxid2annotations = {}

    for node in ete3_feature_tree.traverse("preorder"):
        taxid = node.taxid
        count = node.count
        annotation = [node.sci_name, node.rank]
        if not rank:
            taxid2counts[taxid] = count
            taxid2annotations[taxid] = annotation
        else:
            if node.rank == 'no rank':
                continue
            taxid2counts[taxid] = count
            taxid2annotations[taxid] = annotation

    return taxid2counts, taxid2annotations


def main():

    # main parser
    parser = argparse.ArgumentParser(description="plot phylogenetic distribution of GLASSgo output")
    parser.add_argument("input_GLASSgo_Dir", help="input GLASSgo output directory")
    parser.add_argument('-r', '--rank', action='store_true', help='only calculate counts for taxids which have a taxonomy rank')
    parser.add_argument('-s', '--suffix', help="suffix of GLASSgo sRNA files [.fa]", default='.fa')
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory [./]", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    # show usage
    if len(sys.argv) < 2:
        sys.stderr.write("\nError: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    # parse args
    args = parser.parse_args()

    # input and output handling
    if not args.prefix:
        args.prefix = "GLASSgo_sRNA_distribution"
        output_count = os.path.join(args.out_folder, args.prefix+".tsv")

    if os.path.exists(output_count):
        if args.force:
            print("Warning: output file exists, will be overwritten!")
        else:
            sys.stderr.write("Error: output file detected, please remove it or use --force option to overwrite it")

    input_dir = os.path.abspath(args.input_GLASSgo_Dir)

    all_taxid2annotations = {} # {taxid1:anno1}
    sRNA_taxid2counts = OrderedDict({}) # {sRNA1:taxid2counts, sRNA2:taxid2counts}
    for sRNA in os.listdir(input_dir):
        if sRNA.endswith(args.suffix):
            tree = get_feature_tree_from_GLASSgo_sRNAs(os.path.join(input_dir, sRNA))
            taxid2counts, taxid2annotations = get_taxid2counts(tree, rank=args.rank)
            # get taxid2count for this sRNA
            sRNA_name = sRNA.lstrip('GLASSgo_output_').rstrip(args.suffix)
            sRNA_taxid2counts[sRNA_name] = taxid2counts
            # update all_taxid2annotations
            for taxid in taxid2counts.keys():
                if taxid not in all_taxid2annotations:
                    all_taxid2annotations[taxid] = taxid2annotations[taxid]

    with open(output_count, 'w') as count_oh:
        count_oh.write('TaxID' +"\t"+ 'sci_name' +"\t"+ 'rank' +"\t"+ "\t".join(sRNA_taxid2counts.keys()) + "\n")
        for taxid in sorted(all_taxid2annotations.keys()):
            count_line = [str(taxid) +"\t"+ "\t".join(all_taxid2annotations[taxid])]
            for sRNA, taxid2counts in sRNA_taxid2counts.items():
                sRNA_count = taxid2counts.get(taxid, 0)
                count_line.append(str(sRNA_count))
            count_oh.write("\t".join(count_line)+"\n")

if __name__ == "__main__":
    main()

