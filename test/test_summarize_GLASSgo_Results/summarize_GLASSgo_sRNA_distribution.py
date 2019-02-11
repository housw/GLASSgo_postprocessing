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


def get_taxid2counts_from_GLASSgo_sRNAs(GLAssgo_sRNA_file):
    """ a sample of GLASSgo output: 
        >HG001_02899
        CGTATACAAGGATAAAGCTTATAACAGTAGTAATTGTTGCTATCAAACGAACAACATATATTCTATTTTCAGATAGCAA
        >HE681097.1:c736762-736685 Staphylococcus aureus subsp. aureus HO 5096 0412 complete genome-p.c.VAL:96.25%-taxID:1074252
        ATATACAAGGATAAAGCTTATAACAGTAGTAATTGTTGCTATCAAACGAACAACATATATTCTATTTTCAGATAGCAA
    """

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
    return taxid2counts


def main():

    # main parser
    parser = argparse.ArgumentParser(description="plot phylogenetic distribution of GLASSgo output")
    parser.add_argument("input_GLASSgo_Dir", help="input GLASSgo output directory")
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
        out_file = os.path.join(args.out_folder, args.prefix+".tsv")

    if os.path.exists(out_file):
        if args.force:
            print("Warning: output file exists, will be overwritten!")
        else:
            sys.stderr.write("Error: output file detected, please remove it or use --force option to overwrite it")

    input_dir = os.path.abspath(args.input_GLASSgo_Dir)

    all_taxids = []
    sRNA_taxid2counts = OrderedDict({}) # {sRNA1:taxid2counts, sRNA2:taxid2counts}
    for sRNA in os.listdir(input_dir):
        if sRNA.endswith(args.suffix):
            taxid2counts = get_taxid2counts_from_GLASSgo_sRNAs(os.path.join(input_dir, sRNA))
            sRNA_name = sRNA.lstrip('GLASSgo_output_').rstrip(args.suffix)
            #print(sRNA_name, taxid2counts)
            for taxid in taxid2counts.keys():
                all_taxids.append(taxid)
            if sRNA_name not in sRNA_taxid2counts:
                sRNA_taxid2counts[sRNA_name] = taxid2counts
            else:
                raise Exception("Duplicate sRNA names are found: {name}".format(name=sRNA_name))

    with open(out_file, 'w') as oh:
        oh.write('TaxID' + "\t" + "\t".join(sRNA_taxid2counts.keys()) + "\n")
        for taxid in sorted(set(all_taxids)):
            line = [taxid]
            for sRNA, taxid2counts in sRNA_taxid2counts.items():
                sRNA_count = taxid2counts.get(taxid, 0)
                #print(taxid, sRNA, sRNA_count)
                line.append(str(sRNA_count))
            oh.write("\t".join(line)+"\n")


if __name__ == "__main__":
    main()

