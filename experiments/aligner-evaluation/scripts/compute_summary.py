from argparse import ArgumentParser, RawTextHelpFormatter
from gzip import GzipFile
from google.protobuf.internal.decoder import _DecodeVarint32
from vg_pb2 import Alignment
from numpy import cumsum
from joblib import Parallel, delayed
from edlib import align
from bisect import bisect
from re import split


def load_graph(gfa_graph):

    vertex_labels, edges = dict(), dict()

    for line in open(gfa_graph).read().split('\n')[:-1]:
        if line[0] == 'S':
            str_id, label = line[1:].strip().split()
            vertex_labels[str_id] = label
        if line[0] == 'L':
            tail_str_id, _, head_str_id, _, _ = line[1:].strip().split()
            tail_id, head_id = tail_str_id, head_str_id
            if tail_id not in edges:
                edges[tail_id] = list()
            edges[tail_id].append(head_id)

    return vertex_labels, edges


def get_read_info(read_label, read_header):

    s, t = map(int, read_header.split()[1].split(',')[-1].split('-'))
    is_reverse_comp = '-strand' in read_header.split()[1].split(',')

    return read_label, is_reverse_comp, s, t


def load_reads_and_ref(fastq, fasta, path):

    fastq_lines = open(fastq).read().split('\n')[:-1]

    ref = open(fasta).readlines()[-1].strip() if fasta else ''
    ref_rev_comp = ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}[b] for b in ref[::-1])
    path = [s.strip().split()[1] for s in open(path).readlines()] if path else list()

    if path:
        reads = {
            read_header.strip().split()[0][1:]: get_read_info(read_label.strip(), read_header.strip())
            for read_header, read_label in zip(fastq_lines[::4], fastq_lines[1::4])
        }
    else:
        reads = {
            read_header.strip().split()[0][1:]: (read_label.strip(), False, -1, -1)
            for read_header, read_label in zip(fastq_lines[::4], fastq_lines[1::4])
        }

    return reads, ref, ref_rev_comp, path


def read_gam(gam_filename):

    with open(gam_filename, 'rb') as f:
        buf = GzipFile(fileobj=f).read()
        n = 0
        while n < len(buf):
            an, n = _DecodeVarint32(buf, n)
            for i in range(an):
                msg_len, n = _DecodeVarint32(buf, n)
                msg_buf = buf[n:n + msg_len]
                n += msg_len
                aln = Alignment()
                aln.ParseFromString(msg_buf)
                yield aln


def parse_gam(raw_gam, vertex_labels):

    name = raw_gam.name.split()[0]
    rev_cnt = 0
    path = list()
    idx, n = 0, len(raw_gam.path.mapping)
    first_node_off, last_node_off = 0, 0

    seqs = list()

    for x in raw_gam.path.mapping:

        node_name = x.position.name
        if node_name =='':
            node_name = str(x.position.node_id)
        ll = vertex_labels[node_name]
        original_len = len(ll)

        if x.position.is_reverse:
            rev_cnt += 1
        if idx == 0:
            if rev_cnt > 0:
                first_node_off = original_len - x.position.offset
            else:
                first_node_off = x.position.offset

        if idx == n - 1:
            suma = sum(i.from_length for i in x.edit)
            if rev_cnt > 0:
                last_node_off = original_len - suma - (x.position.offset if idx == 0 else 0)
            else:
                last_node_off = suma + (x.position.offset if idx == 0 else 0)

        if idx == 0 and idx == n - 1:
            if rev_cnt > 0:
                ll = ll[last_node_off:first_node_off]
            else:
                ll = ll[first_node_off:last_node_off]
        elif idx == 0:
            if rev_cnt > 0:
                ll = ll[:first_node_off]
            else:
                ll = ll[first_node_off:]
        elif idx == n - 1:
            if rev_cnt > 0:
                ll = ll[last_node_off:]
            else:
                ll = ll[:last_node_off]

        if rev_cnt > 0:
            ll = ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}[b] for b in ll[::-1])

        path.append(node_name)
        seqs.append(ll)
        idx += 1

    seq = ''.join(seqs)

    return name, seq, len(path), rev_cnt, len(seq), path, first_node_off, last_node_off


def load_gam(gam_filename, vertex_labels):

    ret = dict()

    for raw_gam in read_gam(gam_filename):
        a = parse_gam(raw_gam, vertex_labels)
        if a[0] not in ret:
            ret[a[0]] = a
        else:
            # take the longer one when multiple alignments
            if a[-4] > ret[a[0]][-4]:
                ret[a[0]] = a

    return ret


def read_gaf(gaf_filename):

    for line in open(gaf_filename, 'r').read().split('\n')[:-1]:
        items = line.split('\t')
        raw_path = items[5]
        path = split('<|>', raw_path)[1:]
        rev = '<' in raw_path

        yield items[0], path, rev, int(items[7]), int(items[8])


def parse_gaf(raw_gaf, vertex_labels):

    id, path, rev, f_o, l_o = raw_gaf

    if rev:
        rev_cnt = len(path)
        first_node_off, last_node_off = len(vertex_labels[path[0]]) - f_o, len(vertex_labels[path[-1]]) - l_o
    else:
        rev_cnt = 0
        first_node_off, last_node_off = f_o, l_o


    seqs = list()
    n = len(path)

    for idx, node_id in enumerate(path):

        ll = vertex_labels[node_id]
        original_length = len(ll)

        if idx == 0 and idx == n - 1:
            if rev_cnt > 0:
                ll = ll[last_node_off:first_node_off]
            else:
                ll = ll[first_node_off:last_node_off]
        elif idx == 0:
            if rev_cnt > 0:
                ll = ll[:first_node_off]
            else:
                ll = ll[first_node_off:]
        elif idx == n - 1:
            if rev_cnt > 0:
                ll = ll[last_node_off:]
            else:
                ll = ll[:last_node_off]

        if rev_cnt > 0:
            ll = ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}[b] for b in ll[::-1])
            if idx < n - 1:
                last_node_off += original_length
        else:
            if idx < n - 1:
                last_node_off -= original_length

        seqs.append(ll)

    seq = ''.join(seqs)

    return id, seq, n, rev_cnt, len(seq), path, first_node_off, last_node_off


def load_gaf(gaf_filename, vertex_labels):

    ret = dict()

    for raw_gaf in read_gaf(gaf_filename):
        a = parse_gaf(raw_gaf, vertex_labels)
        if a[0] not in ret:
            ret[a[0]] = a
        else:
            # take the longer one when multiple alignments
            if a[-4] > ret[a[0]][-4]:
                ret[a[0]] = a

    return ret


def compute_overlap(read, node):
    return max(0, min(read[1], node[1]) - max(read[0], node[0]))


def compute_summary(args):

    vertex_labels, edges = load_graph(args.graph)
    reads, ref_seq, ref_seq_rev_comp, ref_path = load_reads_and_ref(args.fastq, args.fasta, args.path)
    node_limits = list(cumsum([len(vertex_labels[v]) for v in ref_path]))

    for alignment_filename, metric in zip(args.alignments, args.metrics):
        if alignment_filename.endswith('.gam'):
            alignments = load_gam(alignment_filename, vertex_labels)
        elif alignment_filename.endswith('.gaf'):
            alignments = load_gaf(alignment_filename, vertex_labels)


        csv = open(metric, 'w')
        csv.write('id,len_aln,ed_read,ed_true,overlap,len_truth,len_read\n')

        def compute_edit(read_id):

            read_label, reverse_read, s, t = reads[read_id]
            t = min(t, len(ref_seq))
            true_seq = ''

            if s > 0 and t > 0:  # SIMULATED READ

                if reverse_read:
                    true_seq = ref_seq_rev_comp[s: t]
                    s, t = len(ref_seq) - t, len(ref_seq) - s  ## Transform to original coordinates
                else:
                    true_seq = ref_seq[s: t]

                first_node_index, last_node_index = bisect(node_limits, s), bisect(node_limits, t - 1)
                read_path = ref_path[first_node_index:last_node_index + 1]
                read_path_intervals = [[0, len(vertex_labels[n])] for n in read_path]
                read_path_intervals[0][0] = s - (0 if first_node_index == 0 else node_limits[first_node_index - 1])
                read_path_intervals[-1][-1] = t - (0 if last_node_index == 0 else node_limits[last_node_index - 1])
                read_path_nodes = {read_path[i]: read_path_intervals[i] for i in range(len(read_path))}

            aln_seq = ''
            overlap = 0

            if read_id in alignments:

                a = alignments[read_id]
                aln_seq = a[1]

                if s > 0 and t > 0:

                    first_node_off = a[-2]
                    last_node_off = a[-1]
                    reverse_alignment = a[3] > 0

                    if reverse_alignment == reverse_read:  # Only consider overlap if same direction paths

                        alignment_path = a[5]
                        alignment_path_intervals = [[0, len(vertex_labels[n])] for n in alignment_path]
                        if reverse_alignment:
                            alignment_path_intervals[0][-1] = first_node_off
                            alignment_path_intervals[-1][0] = last_node_off
                        else:
                            alignment_path_intervals[0][0] = first_node_off
                            alignment_path_intervals[-1][-1] = last_node_off


                        alignment_path_nodes = {alignment_path[i]: alignment_path_intervals[i] for i in
                                                range(len(alignment_path))}
                        overlap = sum(compute_overlap(read_path_nodes[x], alignment_path_nodes[x]) for x in a[5] if
                                      x in read_path_nodes)

            row = list()
            row.append(read_id)
            row.append(len(aln_seq))
            row.append(align(read_label, aln_seq, mode='NW')['editDistance'])
            row.append(align(true_seq, aln_seq, mode='NW')['editDistance'])
            row.append(overlap)
            row.append(len(true_seq))
            row.append(len(read_label))

            return row

        reads_ids = [id for id in reads]
        reads_n = len(reads_ids)
        block_size = 200
        n_blocks = reads_n // block_size + 1

        def compute_edit_kernel(tid):
            l, r = tid * reads_n // n_blocks, (tid + 1) * reads_n // n_blocks
            tmp = []
            for i in range(l, r):
                tmp.append(compute_edit(reads_ids[i]))
            return tmp

        tmp_list = Parallel(n_jobs=args.threads, prefer="threads")(
            delayed(compute_edit_kernel)(t_idx) for t_idx in range(n_blocks))
        processed_list = list()
        for c in tmp_list:
            processed_list += c

        for row in processed_list:
            csv.write(','.join(list(map(str, row))) + '\n')
        csv.close()


if __name__ == '__main__':

    parser = ArgumentParser(
        description='''
               Computes distances between alignments and input reads.
               If read are simulated it also computes overlaps and distances to truth.            
            ''',
        formatter_class=RawTextHelpFormatter
    )

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-g', '--graph', type=str, help='Input gfa file', required=True)
    requiredNamed.add_argument('-fq', '--fastq', type=str, help='Input fastq file', required=True)
    requiredNamed.add_argument(
        '-als', '--alignments', type=str, help='Output gam/gaf files (with extension, each)', required=True, nargs='+'
    )
    requiredNamed.add_argument(
        '-mts', '--metrics', type=str, help='Output csv files with metrics', required=True, nargs='+'
    )

    parser.add_argument('-p', '--path', type=str, help='Output path file (node ids of selected path)')
    parser.add_argument('-fa', '--fasta', type=str, help='Output fasta of original path')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads', default=30)

    compute_summary(parser.parse_args())
