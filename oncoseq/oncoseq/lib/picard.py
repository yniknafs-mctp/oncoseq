'''
Created on Nov 2, 2011

@author: mkiyer
'''

def parse_metrics(line_iter):
    for line in line_iter:
        if not line: continue
        if line.startswith('## METRICS'): break
    keys = line_iter.next().strip().split('\t')
    results = []
    for line in line_iter:
        if not line: break
        if not line.strip(): break
        vals = line.strip().split('\t')
        vals[1:] = map(float, vals[1:])
        results.append(vals)
    return keys, results

def get_num_aligned_reads(alignment_summary_metrics_file):
    keys, results = parse_metrics(open(alignment_summary_metrics_file))
    category_index = keys.index("CATEGORY")
    pf_reads_aligned_index = keys.index("PF_READS_ALIGNED")
    num_aligned_reads = 0
    for fields in results:
        category = fields[category_index]
        if category == "PAIR":
            continue
        num_aligned_reads += fields[pf_reads_aligned_index]
    return num_aligned_reads

#def _picard_parse_hist(line_iter):
#    for line in line_iter:
#        if not line: break
#        if not line.strip(): break
#        fields = line.strip().split('\t')
#        bin = int(fields[0])
#        freq = float(fields[1])
#        yield bin, freq
#
#def _picard_parse_cycle_qualities(input_file):
#    fileh = open(input_file)
#    for line in fileh:
#        if not line: continue
#        if not line.strip(): continue
#        if line.startswith('#'): continue
#        if line.startswith('CYCLE'): break
#    results = list(_picard_parse_hist(fileh))
#    fileh.close()
#    return results
#
#def _picard_parse_index_stats(input_file):
#    results = []
#    for line in open(input_file):
#        if line.startswith("NoCoordinateCount"):
#            ref = "Unaligned"
#            length = 0
#            count = int(line.split()[1])
#        else:
#            fields = line.strip().split('\t')
#            ref = fields[0].split()[0]
#            length = int(fields[1])
#            count = int(fields[2].split()[1])
#        results.append((ref, length, count))
#    return results
#
#def _picard_parse_quality_distribution(input_file):
#    fileh = open(input_file)
#    for line in fileh:
#        if not line: continue
#        if not line.strip(): continue
#        if line.startswith('#'): continue
#        if line.startswith('QUALITY'): break
#    results = list(_picard_parse_hist(fileh))
#    fileh.close()
#    return results
#
#def _picard_parse_metrics(line_iter):
#    for line in line_iter:
#        if not line: continue
#        if line.startswith('## METRICS'): break
#    keys = line_iter.next().strip().split('\t')
#    results = []
#    for line in line_iter:
#        if not line: break
#        if not line.strip(): break
#        vals = line.strip().split('\t')
#        vals[1:] = map(float, vals[1:])
#        results.append(vals)
#    return keys, results