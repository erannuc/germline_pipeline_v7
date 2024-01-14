from graph import Graph_struct
import argparse
from pysam import VariantFile
import numpy as np

'''
my_instance = Graph_struct(5)
my_instance.add_edge(1, 0)
my_instance.add_edge(2, 3)
my_instance.add_edge(3, 0)
print("1-->0")
print("2-->3")
print("3-->0")
conn_comp = my_instance.connected_components()
print("The connected components are :")
print(conn_comp)
'''

def representative(c, index2details):
    # choose a representative variant from the connected component c:

    if len(c) == 1:
        return 0
    else:
        return np.ndarray.argmin(np.array([index2details[i]['pb'] for i in c])) # representative is the variant with minimal p_value

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="cluster variants on same haplotype", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', help='input file in standard vcf v4.3 format with input variants. could be gzipped or not')
    parser.add_argument('-o', help='output vcf file in standard vcf v4.3 format with output variants, including tags for cluster representative variants')

    args = parser.parse_args()
    vcf_in = VariantFile(args.i)  # auto-detect input format
    vcf_out = VariantFile(args.i)  # auto-detect input format
    V = {str(i): [] for i in range(1,23)} # dict of lists of all variants. key is chromosome
    index2details = {}
    vcount = 0
    for rec in vcf_in:
        chr = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0]
        ga = rec.info.get('GA')[0]
        pb = rec.info.get('PB')[0]
        op = rec.info.get('OP')[0]

        print(rec)
        print(op)
        V[chr].append((vcount, chr, pos, ref, alt, ga, pb, op))
        index2details[vcount] = {'chr': chr, 'pos': pos, 'ref': ref, 'alt': alt, 'ga': ga, 'pb': pb, 'op': op}
        vcount += 1

    variants_graph = Graph_struct(vcount)

    # now populate the graph with edges
    for c in range(1, 23):
        c = str(c)
        for v in range(len(V[c])):
            for u in range(v+1, len(V[c])):
                if abs(V[c][v][2] - V[c][u][2]) < 10000 and abs(V[c][v][5] - V[c][u][5]) < 0.01:
                    variants_graph.add_edge(V[c][v][0], V[c][u][0])

    conn_comp = variants_graph.connected_components()
    print(conn_comp)
    print(len(conn_comp))

    # choose representative from each cluster

    representative_op_coors = {}
    for c in conn_comp:
        print(c)
        print(len(c))
        print("^^^^")
        rep_cluster_index = representative(c, index2details)
        rep = c[rep_cluster_index]
        print("^^^^", rep_cluster_index, rep)
        representative_op_coors[(index2details[rep]['chr'], index2details[rep]['op'])] = len(c)

    # another pass to print representative variants
    vcf_in = VariantFile(args.i)  # auto-detect input format
    vcf_out = VariantFile(args.o, 'w', header=vcf_in.header)  # auto-detect input format
    for rec in vcf_in:
        chr = rec.chrom
        op = rec.info.get('OP')[0]
        rec.info['RP'] = [1]
        print(chr, op)
        if (chr, op) in representative_op_coors:
            rec.info['RP'] = [1]
            rec.info['CS'] = [representative_op_coors[(chr, op)]]
            print("rep", (chr, op), representative_op_coors[(chr, op)])
        else:
            rec.info['RP'] = [0]
        vcf_out.write(rec)

