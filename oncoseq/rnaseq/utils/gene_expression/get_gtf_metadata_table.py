'''
Created on Feb 24, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import collections

from bx.intervals.cluster import ClusterTree

from oncoseq.rnaseq.lib.gtf import GTFFeature

def get_gtf_metadata(gtf_file, 
                     omit_attrs=None,
                     group_by="gene_id", 
                     feature_type="exon"):
    if omit_attrs is None:
        omit_attrs = []
    # read gtf file and group by gene
    gene_feature_map = collections.defaultdict(lambda: [])
    gene_attrs_set = set()
    for feature in GTFFeature.parse(open(gtf_file)):
        if feature.feature_type != feature_type:
            continue
        feature_id = feature.attrs[group_by]
        gene_feature_map[feature_id].append(feature)
        gene_attrs_set.update(feature.attrs.keys())
    gene_attrs_set.difference_update(omit_attrs)
    gene_attrs_list = sorted(gene_attrs_set)
    metadata_fields = ["tracking_id", "locus", "strand", "num_exons", "transcript_length"] + gene_attrs_list
    metadata_inds = dict((x,i) for i,x in enumerate(metadata_fields))
    metadata_dict = {}
    # output metadata sorted by gene id
    for feature_id,features in gene_feature_map.iteritems():
        # collect attributes for this gene
        attrdict = collections.defaultdict(lambda: set())
        # cluster exons together for each gene
        cluster_tree = ClusterTree(0,1)
        for i,f in enumerate(features):
            cluster_tree.insert(f.start, f.end, i)
            for k,v in f.attrs.iteritems():
                if k in gene_attrs_set:
                    # some attributes have multiple values separated by a comma
                    attrdict[k].update(v.split(','))
        # determine larger exon clusters
        transcript_length = 0
        exon_clusters = []
        for start, end, indexes in cluster_tree.getregions():
            exon_clusters.append((start,end))
            transcript_length += (end - start)
        del cluster_tree
        chrom = features[0].seqid
        locus_start = min(e[0] for e in exon_clusters)
        locus_end = max(e[1] for e in exon_clusters)
        locus_string = "%s:%d-%d" % (chrom, locus_start, locus_end)
        strand = features[0].strand
        num_exons = len(exon_clusters)
        # make metadata row
        metadata = [feature_id, locus_string, strand, num_exons, transcript_length] + ['NA'] * len(gene_attrs_list)
        # get all attributes
        for k,vals in attrdict.iteritems():
            ind = metadata_inds[k]
            metadata[ind] = ','.join(map(str, sorted(vals)))
        metadata_dict[metadata[0]] = metadata
    return metadata_fields, metadata_dict

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--group-by', dest="group_by",
                        choices=["gene_id", "transcript_id"],
                        default="gene_id")
    parser.add_argument('--omit', dest="omit_list",
                        default="exon_number",
                        help="default=%(default)s")
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    if not os.path.exists(args.gtf_file):
        parser.error("gtf file %s not found" % (args.gtf_file))
    omit_list = args.omit_list.split(",")
    # get gene metadata
    logging.info("Reading metadata")
    metadata_fields, metadata_dict = \
        get_gtf_metadata(args.gtf_file, omit_attrs=omit_list)
    # write matrix file
    logging.info("Writing table")
    header_fields = list(metadata_fields)
    print '\t'.join(header_fields)
    for k in sorted(metadata_dict):
        fields = list(map(str, metadata_dict[k]))
        print '\t'.join(fields)
    return 0

if __name__ == '__main__':
    sys.exit(main())
