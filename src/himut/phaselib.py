import time
import pysam
import queue
import bisect
import natsort
import numpy as np
import himut.cslib
import himut.bamlib
import himut.haplib
import multiprocessing as mp
from scipy.stats import binom_test
from collections import defaultdict
from typing import Dict, List, Tuple


def get_edges(
    chrom: str,
    bam_file: str,
    min_bq: int,
    min_mapq: int,
    hpos_lst: List[int],
    hetsnp_lst: List[Tuple[str, int, str, str]],
    hetsnp2hidx: Dict[Tuple[str, int, str, str], int],
) -> Tuple[List[Tuple[int, int]], Dict[Tuple[int, int], np.ndarray]]:

    # counter = 0
    edge2counts = defaultdict(lambda: np.zeros(4))
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for line in alignments.fetch(chrom):
        ccs = himut.bamlib.BAM(line)
        if not ccs.is_primary:
            continue

        if ccs.mapq < min_mapq:
            continue

        _idx = bisect.bisect_right(hpos_lst, ccs.tstart)
        _jdx = bisect.bisect_right(hpos_lst, ccs.tend)
        if _idx == _jdx:
            continue
        elif _jdx - _idx == 1:
            continue

        ccs.cs2tpos2qbase()
        hetsnp_subset_lst = [hetsnp_lst[_kdx] for _kdx in range(_idx, _jdx)]
        for idx, snpidx in enumerate(hetsnp_subset_lst):
            i = hetsnp2hidx[snpidx]
            ibase, ibq = ccs.tpos2qbase[snpidx[0]]
            if ibq < min_bq:
                continue
            idx_state = 0 if ibase == snpidx[1] else 1
            for snpjdx in hetsnp_subset_lst[idx + 1 :]:
                j = hetsnp2hidx[snpjdx]
                jbase, jbq = ccs.tpos2qbase[snpjdx[0]]
                if jbq < min_bq:
                    continue
                jdx_state = 0 if jbase == snpjdx[1] else 1
                if not idx_state and not jdx_state:  # cis1
                    edge2counts[(i, j)][0] += 1
                elif idx_state and jdx_state:  # cis2
                    edge2counts[(i, j)][1] += 1
                elif not idx_state and jdx_state:  # trans1
                    edge2counts[(i, j)][2] += 1
                elif idx_state and not jdx_state:  # trans2
                    edge2counts[(i, j)][3] += 1
    edge_lst = natsort.natsorted(list(edge2counts.keys()))
    return edge_lst, edge2counts


def table2binom_test(table: np.ndarray) -> float:
    total_count = np.sum(table)
    total_cis_count = table[0] + table[1]
    p_value = binom_test(total_cis_count, total_count, p=0.5, alternative="two-sided")
    return p_value


def build_graph(
    edge_lst: List[Tuple[int, int]], edge2counts: Dict[Tuple[int, int], int]
) -> Dict[int, List[int]]:

    graph = defaultdict(list)
    for (i, j) in edge_lst:
        total_count = np.sum(edge2counts[i, j])
        if total_count == 0:
            continue
        graph[i].append(j)
        graph[j].append(i)
    return graph


def get_phased_graph(
    edge_lst: List[Tuple[int, int]],
    edge2counts: Dict[Tuple[int, int], np.ndarray],
    min_p_value: float,
    min_phase_proportion: float,
) -> Dict[int, List[int]]:

    phased_graph = {}
    edge2p_value = {}
    graph = build_graph(edge_lst, edge2counts)
    for i in natsort.natsorted(graph.keys()):
        phased_edges = []
        phase_consistent_edge_count = 0
        for j in graph[i]:
            if i < j:
                p_value = table2binom_test(edge2counts[(i, j)])
                edge2p_value[(i, j)] = p_value
            else:
                p_value = edge2p_value[(j, i)]
            if p_value < min_p_value:
                phased_edges.append(j)
                phase_consistent_edge_count += 1
        total_edge_count = len(graph[i])
        phase_consistent_edge_proportion = (
            phase_consistent_edge_count / total_edge_count
        )
        if phase_consistent_edge_proportion > min_phase_proportion:
            phased_graph[i] = phased_edges
    return phased_graph


def build_haplotype_block(
    edge_lst: List[Tuple[int, int]],
    edge2counts: Dict[Tuple[int, int], np.ndarray],
    min_p_value: float,
    min_phase_proportion: float,
) -> List[List[Tuple[int, int]]]:

    bfs = [{}]
    bfs_idx = 0
    seen = set()
    graph = get_phased_graph(edge_lst, edge2counts, min_p_value, min_phase_proportion)
    for nodeA in graph:  # bfs
        if nodeA in seen:
            continue
        seen.add(nodeA)
        q = queue.Queue()

        edge_lst = []
        bfs[bfs_idx][nodeA] = "0"
        for nodeB in graph[nodeA]:
            k1 = min(nodeA, nodeB)
            k2 = max(nodeA, nodeB)
            table = edge2counts[(k1, k2)]
            total_cis_count = table[0] + table[1]
            total_trans_count = table[2] + table[3]
            if total_cis_count > total_trans_count:
                bfs[bfs_idx][nodeB] = "0"
            else:
                bfs[bfs_idx][nodeB] = "1"
            edge_lst.append((nodeA, nodeB))
        q.queue = queue.deque(edge_lst)

        while not q.empty():
            nodeA, nodeB = q.get()
            if nodeB in seen:
                continue
            if not nodeB in graph:
                continue

            seen.add(nodeB)
            k3 = min(nodeA, nodeB)
            k4 = max(nodeA, nodeB)
            table = edge2counts[(k3, k4)]
            total_cis_count = table[0] + table[1]
            total_trans_count = table[2] + table[3]
            hap_idx = bfs[bfs_idx][nodeA]
            if total_cis_count > total_trans_count:
                if hap_idx == "0":
                    bfs[bfs_idx][nodeB] = "0"
                else:
                    bfs[bfs_idx][nodeB] = "1"
            else:
                if hap_idx == "0":
                    bfs[bfs_idx][nodeB] = "1"
                else:
                    bfs[bfs_idx][nodeB] = "0"

            for nodeC in graph[nodeB]:
                if nodeC in seen:
                    continue
                q.put((nodeB, nodeC))
        bfs_idx += 1
        bfs.append({})

    hblock_lst = []
    for idx in range(len(bfs)):
        hblock = []
        for node in bfs[idx]:
            hblock.append((node, bfs[idx][node]))
        if len(hblock) == 0 or len(hblock) == 1:
            continue
        hblock = natsort.natsorted(hblock)
        hblock_lst.append(hblock)
    return hblock_lst


def get_hblock_statistics(
    hblock_lst: List[List[Tuple[int, int]]],
    hetsnp_lst: List[Tuple[str, int, str, str]],
) -> Tuple[int, int, int, int, int]:

    block_len_lst = []
    block_distance_lst = []
    for hblock in hblock_lst:
        block_len = len(hblock)
        posidx = hetsnp_lst[hblock[0][0]][0]
        posjdx = hetsnp_lst[hblock[-1][0]][0]
        block_distance = posjdx - posidx
        block_len_lst.append(block_len)
        block_distance_lst.append(block_distance)

    num_hetsnps = len(hetsnp_lst)
    num_blocks = len(hblock_lst)
    if len(block_len_lst) == 0:
        largest_block = 0
        smallest_block = 0
        longest_block = 0
        shortest_block = 0
    else:
        largest_block = max(block_len_lst)
        smallest_block = min(block_len_lst)
        longest_block = max(block_distance_lst)
        shortest_block = min(block_distance_lst)
    return (
        num_hetsnps,
        num_blocks,
        smallest_block,
        largest_block,
        shortest_block,
        longest_block,
    )


def get_hblock(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    vcf_file: str,
    min_bq: int,
    min_mapq: int,
    min_p_value: float,
    min_phase_proportion: float,
    chrom2hblock_lst: Dict[str, List[List[Tuple[int, int]]]],
) -> List[List[Tuple[int, Tuple[int, int]]]]:

    hetsnp_lst, _, hetsnp2hidx = himut.vcflib.load_hetsnps(vcf_file, chrom, chrom_len)
    hpos_lst = [hetsnp[0] for hetsnp in hetsnp_lst]
    edge_lst, edge2counts = get_edges(
        chrom, bam_file, min_bq, min_mapq, hpos_lst, hetsnp_lst, hetsnp2hidx
    )
    chrom2hblock_lst[chrom] = build_haplotype_block(
        edge_lst, edge2counts, min_p_value, min_phase_proportion
    )


def get_chrom_hblock(
    bam_file: str,
    vcf_file: str,
    region: str,
    region_list: str,
    min_bq: int,
    min_mapq: int,
    min_p_value: float,
    min_phase_proportion: float,
    threads: int,
    version: str,
    out_file: str,
) -> Tuple[
    Dict[str, List[List[Tuple[int, int]]]], Dict[str, Tuple[int, int, int, int, int]]
]:
    cpu_start = time.time() / 60
    print("himut is phasing hetsnps with {} threads".format(threads))
    himut.util.check_num_threads(threads)
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, _ = himut.util.load_loci(region, region_list, tname2tsize)
    himut.util.check_phaser_input_exists(
        bam_file, vcf_file, out_file, chrom_lst, tname2tsize
    )
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2hblock_lst = manager.dict()
    get_hblock_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            vcf_file,
            min_bq,
            min_mapq,
            min_p_value,
            min_phase_proportion,
            chrom2hblock_lst,
        )
        for chrom in chrom_lst
    ]
    p.starmap(get_hblock, get_hblock_arg_lst)
    p.close()
    p.join()
    print("himut finished phasing hetsnps")

    print("himut is returning phased hetsnps")
    himut.vcflib.dump_phased_hetsnps(
        bam_file,
        vcf_file,
        region,
        region_list,
        tname2tsize,
        min_bq, 
        min_mapq,
        min_p_value,
        min_phase_proportion,
        threads,
        chrom_lst,
        chrom2hblock_lst,
        version,
        out_file,
    )
    print("himut finished returning phased hetsnps")
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("himut haplotype phasing took {} minutes".format(duration))
    himut.util.exit()
