#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Create by ygidtu@gmail.com at 2019.01.14

Collect SJ.out.tab results and calculate the PSI
"""
import sys
import argparse as ap
import os
import logging
import re
import glob

try:
    import pysam
except ImportError as err:
    pysam = None

from tqdm import tqdm

from multiprocessing import Pool, cpu_count, Manager


def set_logging(log_name):
    u"""
    set logging handler
    Created by Zhang yiming at 2018.11.13
    :return:
    """
    sh = logging.StreamHandler()

    # [%(module)s:%(lineno)d]
    formatter = logging.Formatter(
        fmt="[%(asctime)s] - [%(levelname)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    sh.setFormatter(formatter)
    sh.setLevel(logging.INFO)

    log = logging.getLogger(log_name)
    log.setLevel(logging.DEBUG)
    log.addHandler(sh)
    return log


logger = set_logging("")


"""
Code to collect the results from SJ.out.tab file
"""


def __read_sj_out_tab__(infile, genome_fa, uniq=False):
    u"""
    @2018.12.25
    read data from sj out tab
    :param infile:
    :param genome_fa:
    :param uniq: bool
    :return:
    """
    strand = {
        "1": "+",
        "2": "-",
        "0": "."
    }
    motifs = {
        "0": None,
        "1": "GT/AG",
        "2": "CT/AC",
        "3": "GC/AG",
        "4": "CT/GC",
        "5": "AT/AC",
        "6": "CT/AT"
    }
    uniq_count = {}
    multi_count = {}
    motif = {}

    with open(infile) as r:
        for line in r:
            lines = line.split()
            try:
                strand_ = strand[lines[3]]
            except KeyError:
                strand_ = __determine_strand__(
                    chromosome=lines[0],
                    start=int(lines[1]),
                    end=int(lines[2]),
                    genome_fa=genome_fa
                )

            junctions = "%s:%s-%s:%s" % (lines[0], lines[1], lines[2], strand_)
            uniq_count[junctions] = int(lines[6])

            if not uniq:
                multi_count[junctions] = int(lines[7]) + int(lines[6])

                motif[junctions] = {
                    "motif_code": int(lines[4]),
                    "motif_seq": motifs[lines[4]],
                    "annotation": lines[5],
                    "overhang": int(lines[8])
                }

    if uniq:
        return uniq_count
    return uniq_count, multi_count, motif


def __determine_strand__(chromosome, start, end, genome_fa):
    u"""
    check if the base around junctions are in GT/AG(+); CT/AC(-); GC/AG(+); CT/GC(-); AT/AC(+); GT/AT(-)
    :param chromosome:
    :param start:
    :param end:
    :param genome_fa:
    :return:
    """
    if pysam is None:
        return "."

    if genome_fa and os.path.exists(genome_fa):
        index = genome_fa + ".fai"
        if not os.path.exists(index):
            index = None
        elif os.path.getctime(index) < os.path.getctime(genome_fa):
            os.remove(index)
            index = None

        if not index:
            logger.info("Create index for %s" % genome_fa)
            pysam.faidx(genome_fa)

        with pysam.FastaFile(genome_fa) as r:

            for i in (start, end):
                sequence = ""
                for j in r.fetch(reference=chromosome, start=i - 1, end=i + 1):
                    sequence += j

                if re.search(r"(GT|AT|GC)", sequence) or re.search(r"AG", sequence):
                    return "+"
                elif re.search(r"(CT)", sequence) or re.search(r"GC|AT", sequence):
                    return "-"

    return "."


def collect_sj_out_tab(input_directory, output_file, genome_fa):
    u"""
    @2018.12.25
    Collect SJ.out.tab file data, and merge it into json
    :return:
    """
    files = glob.glob(os.path.join(input_directory, "*SJ.out.tab"))

    data = {
        "uniq_count": {},
        "multi_count": {},
        "motif": {},
    }

    header = ["junctions"]
    for idx, f in enumerate(sorted(files)):
        logger.info("%d/%d Reading %s" % (idx + 1, len(files), f))
        label = re.sub(r"[_\.]?SJ.out.tab", "", os.path.basename(f))

        header.append(label)

        uniq_count, multi_count, motif = __read_sj_out_tab__(f, genome_fa=genome_fa)

        data["uniq_count"][label] = uniq_count
        data["multi_count"][label] = multi_count

        for k, v in motif.items():
            if k in data["motif"].keys():
                v["overhang"] = max(data["motif"][k]["overhang"], v["overhang"])
            data["motif"][k] = v

    rows = set()
    for i in header[1:]:
        for j in data["uniq_count"][i].keys():
            rows.add(j)

    rows = sorted(list(rows))

    def write_counts_to_file(rows, header, file_name, data):
        with open(file_name, "w+") as w:
            w.write("\t".join(header) + "\n")
            for i in rows:
                tmp_row_data = [i]
                for j in header[1:]:
                    try:
                        tmp_row_data.append(str(data[j][i]))
                    except KeyError:
                        tmp_row_data.append("NA")

                w.write("\t".join(tmp_row_data) + "\n")

    logger.info("Write uniq counts")
    write_counts_to_file(
        rows,
        header,
        output_file + "_uniq_counts.tsv",
        data["uniq_count"]
    )

    logger.info("Write multi counts")
    write_counts_to_file(
        rows,
        header,
        output_file + "_multi_counts.tsv",
        data["multi_count"]
    )

    logger.info("Write motif")
    with open(output_file + "_motif_overhang.tsv", "w+") as w:
        w.write("junctions\tannotation\tmotif_code\tmotif_seq\toverhang\n")
        for i in rows:
            try:
                w.write("%s\t%s\t%s\t%s\t%s\n" % (
                    i,
                    str(data["motif"][i]["annotation"]),
                    str(data["motif"][i]["motif_code"]),
                    str(data["motif"][i]["motif_seq"]),
                    str(data["motif"][i]["overhang"])
                ))
            except KeyError:
                continue

"""
Calculate PSI
"""


class GenomicLoci(object):
    u"""
    Basic class of exon, gene and transcripts
    Created by Zhang yiming at 2018.11.13
    """

    __slots__ = [
        "chromosome",
        "start",
        "end",
        "strand"
    ]

    def __init__(self, chromosome, start, end, strand):
        u"""
        init this class by chromosome, start, end and strand
        :param chromosome: str
        :param start: int
        :param end: int
        :param strand: str
        """
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        if self.end < self.start:
            raise ValueError("%d should bigger than %d" % (self.end, self.start))

        if self.strand not in ("+", "-", "."):
            raise ValueError("%s error, strand should be + or -" % self.strand)

    @classmethod
    def create_from_string(cls, string):
        u"""
        Create genomic loci from string,
        :param string: chr1:100-1000:+
        :return: GenomicLoci
        """
        chromosome, sites, strand = string.split(":")
        sites = [int(x) for x in sites.split("-")]

        return cls(chromosome, sites[0], sites[1], strand)

    def __lt__(self, other):
        u"""
        less than, this is different with self.is_upstream()

        2:100-200 < 1:100-200
        1:101-200 < 1:100-200
        1:100-200 < 1:100-199

        :param other: another GenomicLoci or it's children class
        :return: Boolean
        """
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome

        if self.start != other.start:
            return self.start < other.start

        return self.end < other.end

    def __gt__(self, other):
        u"""
        greater than, this is different with self.is_downstream()

        2:100-200 > 1:100-200
        1:101-200 > 1:100-200
        1:100-200 > 1:100-199

        :param other: another GenomicLoci or it's children class
        :return: Boolean
        """
        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome

        if self.start != other.start:
            return self.start > other.start

        return self.end > other.end

    def __eq__(self, other):
        u"""
        check two GenomicLoci is same
        Based on __hash__, so different classes normally will not be the equal, even the sites are same
        :param other: another GenomicLoci
        :return: Boolean
        """
        return self.__hash__() == other.__hash__()

    def __and__(self, other):
        u"""
        check whether two GenomicLoci have any overlap
        :param other:
        :return:
        """
        if not (isinstance(other, GenomicLoci) or issubclass(other, GenomicLoci)):
            raise ValueError("other should be GenomicLoci or it's subclass, not %s" % type(other))

        if self.chromosome == other.chromosome and self.strand == other.strand or \
                (self.chromosome == "all" or other.chromosome == "all"):
            return self.start <= other.end and self.end >= other.start
        return False

    def __add__(self, other):
        u"""
        merge two overlapped GenomicLoci together


        eg: 1:100-200 + 1:150-300 -> 1:100-300

            but, 1:100-200 + 1:201-300 will raise the Value Error

        :param other: another GenomicLoci or it's children class
        :return: new GenomicLoci
                - this function will NOT change any data inside existed object
                - this function will NOT keep any further data inside of child class
                    eg: Transcript will lost all the message of it's id, host genes;
                        Exon will lost all the messages of it's id and host information too
        """
        if self.chromosome != other.chromosome or \
                self.strand != other.strand:
            raise ValueError("Both should be same chromosome and same strand")

        return GenomicLoci(
            chromosome=self.chromosome,
            strand=self.strand,
            start=min(self.start, other.start),
            end=max(self.end, other.end)
        )

    def __hash__(self):
        u"""
        generate hash
        :return: hash
        """
        return hash((self.chromosome, self.start, self.end, self.strand))

    def __str__(self):
        u"""
        convert this class into string
        :return:
        """
        return "%s:%d-%d:%s" % (self.chromosome, self.start, self.end, self.strand)

    @property
    def length(self):
        u"""
        provide a length property
        :return:
        """
        return self.end - self.start

    def is_upstream(self, other):
        u"""
        check if self is upstream of other

        1:100-200 is upstream of 1:201-300

        :param other: another GenomicLoci or it's subclass
        :return: Boolean
        """
        if not (isinstance(other, GenomicLoci) or issubclass(other, GenomicLoci)):
            raise ValueError("other should be GenomicLoci or it's subclass, not %s" % type(other))

        if self.chromosome != other.chromosome and not (self.chromosome == "all" or other.chromosome == "all"):
            return self.chromosome < other.chromosome

        return self.end < other.start

    def is_downstream(self, other):
        u"""
        check if self is downstream of other

        1:201-300 is down stream of 1:100-200

        :param other: another GenomicLoci or it's subclass
        :return: Boolean
        """
        if not (isinstance(other, GenomicLoci) or issubclass(other, GenomicLoci)):
            raise ValueError("other should be GenomicLoci or it's subclass, not %s" % type(other))

        if self.chromosome != other.chromosome and not (self.chromosome == "all" or other.chromosome == "all"):
            return self.chromosome > other.chromosome

        return self.start > other.end

    def overlap_level(self, other, by_narrow=False, by_all=True):
        u"""
        calculate overlap level between two GenomicLoci or it's subclasses, based on self loci
        :param other: another GenomicLoci or it's subclass
        :param by_narrow: calculate the overlap_level by the smaller region
                        eg:1:100-200 + 1:150-300
                            False -> the overlap level is (200 - 150)/(200 - 100) -> 50/100
                            True -> the narrow one is 1:100-200 (200 - 100 < 300 - 150),
                                    then overlap level is (200-150)/(200-100) -> 50/100

        :param by_all: calculate the overlap_level based on the total region of the loci,
                        eg: 1:100-200 + 1:150-300
                            False -> the overlap level is (200 - 150)/(200 - 100) -> 50/100
                            True -> the overlap level is (200 - 150)/(300-100) -> 50/200

        :return: Float, is not any overlap return None
        """
        if not (isinstance(other, GenomicLoci) or issubclass(other, GenomicLoci)):
            raise ValueError("other should be GenomicLoci or it's subclass, not %s" % type(other))

        overlap = self.overlap_distance(other)

        if overlap is not None:
            if by_narrow is True:
                return overlap / min([self.end - self.start, other.end - other.start])

            if by_all is True:
                return overlap / (max([self.end, other.end]) - min(self.start, other.start))
            return overlap / (self.end - self.start)
        return overlap

    def overlap_distance(self, other):
        u"""
        calculate the overlap distance (bp) between two GenomicLoci or it's subclasses
        :param other: another GenomicLoci or it's subclass
        :return: Int, if not overlap return None
        """
        if not (isinstance(other, GenomicLoci) or issubclass(other, GenomicLoci)):
            raise ValueError("other should be GenomicLoci or it's subclass, not %s" % type(other))

        if self.chromosome == "all" or other.chromosome == "all":
            return min(self.end, other.end) - max(self.start, other.start)

        if self.chromosome == other.chromosome and \
                self.start <= other.end and \
                self.end >= other.start:
            return min(self.end, other.end) - max(self.start, other.start)
        return None


class DirectedGraph(object):
    u"""
    Directed and weighted graph  class
    """

    __slots__ = [
        "__same_starts__",
        "__same_ends__",
        "__edges__"
    ]

    def __init__(self):
        u"""
        init this class
        """
        self.__edges__ = {}
        self.__same_starts__ = {}
        self.__same_ends__ = {}

    @property
    def edges(self):
        u"""
        return weighted edges
        :return:
        """
        return self.__edges__

    @property
    def edges_in_loci(self):
        u"""
        return list of GenomicLoci of edges
        :return:
        """
        edges = []
        for i in self.edges.keys():
            tmp = GenomicLoci(
                chromosome="all",
                start=i[0],
                end=i[1],
                strand="+"
            )
            edges.append(tmp)
        return edges

    @property
    def neighbors(self):
        u"""
        get all the neighbors
        :return: list of list, [[1, 2...], [2, 3...]]
        """
        return [self.get_neighbors(x) for x in self.nodes]

    @property
    def nodes(self):
        u"""
        get all the nodes inside this graph
        :return:
        """
        data = set(self.__same_starts__.keys()) | set(self.__same_ends__.keys())
        return sorted(list(data))

    def add_edge(self, source, target, weight=0):
        u"""
        add new edge to this graph
        :param source:
        :param target:
        :param weight:
        :return:
        """
        new_edge = (source, target)
        self.__edges__[new_edge] = weight

        source_node = set() if source not in self.__same_starts__.keys() else self.__same_starts__[source]
        target_node = set() if target not in self.__same_ends__.keys() else self.__same_ends__[target]

        source_node.add(target)
        target_node.add(source)

        self.__same_starts__[source] = source_node
        self.__same_ends__[target] = target_node

    def get_weight(self, target):
        u"""
        get weight of edge
        :param target: (int, int)
        :return: weight, int
        """
        if not isinstance(target, (int, int, tuple)):
            raise ValueError("Target should be tuple of two int, not %s" % type(target))
        return self.__edges__.get(target)

    def get_neighbors(self, target):
        u"""
        get neighbors by node
        :param target: int node
        :return: list of int (neighbors)
        """
        if not isinstance(target, int):
            raise ValueError("Target should be int, not %s" % type(target))

        res = []
        if target in self.__same_starts__.keys():
            res += self.__same_starts__[target]

        if target in self.__same_ends__.keys():
            res += self.__same_ends__[target]

        return res

    def get_same_starts(self, target):
        u"""
        return different end sites with same start
        :param target: The `same` start site
        :return: list of int (different end sites)
        """
        if not isinstance(target, int):
            raise ValueError("Target should be int, not %s" % type(target))
        return list(self.__same_starts__[target]) if target in self.__same_starts__.keys() else []

    def get_same_ends(self, target):
        u"""
        return different start sites with same start
        :param target: The `same` end site
        :return: list of int (different start sites)
        """
        if not isinstance(target, int):
            raise ValueError("Target should be int, not %s" % type(target))
        return list(self.__same_ends__[target]) if target in self.__same_ends__.keys() else []


class JunctionGraph(object):
    u"""
    create junction graph based on networkX
    Created by Zhang yiming at 2018.11.14
    """

    __class__ = [
        "graphs",
        "distance_tolerance",
        "__process__"
        "__overlap_percentage__"
    ]

    def __init__(
            self,
            process=1,
            alternative_threshold=5,
    ):
        u"""
        init this class
        :param alternative_threshold: default 5, filter junctions pairs with sum of counts
        :param process: cpu usage
        """
        # graphs separated by chromosome and strand
        self.graphs = {}
        self.__process__ = process
        self.alternative_threshold = alternative_threshold

    def add_junctions(self, data, threshold=None):
        u"""
        add junction to this graph
        :param data: dict -> {GenomicLoci, int}
        :param threshold: filtered junctions
        :param quite: whether to show progressbar
        :return:
        """
        if not isinstance(data, dict):
            raise ValueError("Should be dict")

        if threshold and not isinstance(threshold, (GenomicLoci, set)):
            raise ValueError("threshold should be set of GenomicLoci")

        for site, count in data.items():

            if threshold and site in threshold:
                continue

            key = "%s#%s" % (site.chromosome, site.strand)

            tmp_graph = self.graphs[key] if key in self.graphs.keys() else DirectedGraph()

            tmp_graph.add_edge(site.start, site.end, count)

            self.graphs[key] = tmp_graph

    def clear(self):
        u"""
        clear all the data in this class, make this class recycle
        :return:
        """
        self.graphs.clear()

    @staticmethod
    def __get_all_psi_per_graph__(graph, key, alternative_threshold):
        u"""
        get PSI table of single graph
        :param graph: DiGraph
        :param key: chromosome#strand
        :param alternative_threshold: int, used to filter the junction pairs with low sum of counts
        :return:
        """
        if not isinstance(graph, DirectedGraph):
            raise ValueError("graph should be DirectedGraph, not %s" % type(graph))

        data_same_start = {}
        data_same_end = {}
        chromosome, strand = key.split("#")
        for i in graph.nodes:
            neighbors = graph.get_neighbors(i)

            same_starts, same_ends = {}, {}

            for j in neighbors:

                if j > i:
                    tmp_loci = GenomicLoci(
                        chromosome=chromosome,
                        strand=strand,
                        start=i,
                        end=j
                    )
                    weight = graph.get_weight((i, j))
                    same_starts[tmp_loci] = weight if weight is not None else 0
                elif j < i:
                    tmp_loci = GenomicLoci(
                        chromosome=chromosome,
                        strand=strand,
                        start=j,
                        end=i
                    )
                    weight = graph.get_weight((j, i))
                    same_ends[tmp_loci] = weight if weight is not None else 0

            total = sum(same_starts.values())
            for key, value in same_starts.items():
                if total < alternative_threshold:
                    continue
                data_same_start[key] = "%.3f" % (value / total)

            total = sum(same_ends.values())
            for key, value in same_ends.items():
                if total < alternative_threshold:
                    continue
                data_same_end[key] = "%.3f" % (value / total)

        return data_same_start, data_same_end

    def get_all_psi(self):
        u"""
        calculate each PSI of same starts and same ends
        :return: same starts PSI DataFrame, same ends PSI DataFrame
        """
        same_start, same_end = {}, {}

        for key, graph in self.graphs.items():
            res = self.__get_all_psi_per_graph__(
                graph, key, self.alternative_threshold
            )

            same_start.update(res[0])
            same_end.update(res[1])

        #     res_by_process = []
        #     with Pool(processes=self.__process__) as pool:
        #         for key, graph in self.graphs.items():
        #             res_by_process.append(
        #                 pool.apply_async(
        #                     self.__get_all_psi_per_graph__,
        #                     args=(graph, key, self.alternative_threshold, )
        #                 )
        #             )
        #
        #         for i in res_by_process:
        #             i = i.get()
        #             same_start.update(i[0])
        #             same_end.update(i[1])

        return same_start, same_end


class CalPSI(object):
    """
    Class to handle the PSI calculation issues
    """

    def __init__(self, args):
        u"""
        init this class with argparse arguments
        """
        if os.path.exists(args.input_dir):
            self.input_directory = os.path.abspath(args.input_dir)

        output_dir = os.path.dirname(args.output_prefix)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        self.output_prefix = os.path.abspath(args.output_prefix)

        self.single_threshold = self.__check_thresholds__(args.single_threshold, 2)
        self.overall_threshold = self.__check_thresholds__(args.overall_threshold, 100)
        self.alternative_threshold = self.__check_thresholds__(args.alternative_threshold, 5)
        self.processes = self.__check_thresholds__(args.processes, cpu_count())
        self.genome = args.genome_fasta

    @staticmethod
    def __check_thresholds__(for_check, default):
        u"""
        if the input parameter is out of range, eg: < 0, return default
        :param for_check: int
        :param default: int
        :return: int
        """

        if for_check < 0:
            logger.warning("negative value %d was passed, using default value now" % for_check)
            return default

        return for_check

    @staticmethod
    def __read_single__(args):
        u"""
        Read single SJ.out.tab file
        :return:
        """
        strand = {
            "1": "+",
            "2": "-",
            "0": "."
        }

        args["lock"].acquire()
        logger.info("%d/%d Reading %s" % (args["idx"], args["total"], args["file_path"]))
        args["lock"].release()
        label = re.sub(r"[_\.]?SJ.out.tab", "", os.path.basename(args["file_path"]))

        # {sample: {GenomicLoci: int}}
        # {GenomicLoci: int}  -> filter low abundance junctions across samples
        tmp, tmp_threshold = {}, {}
        with open(args["file_path"]) as r:
            for line in r:
                lines = line.split()
                try:
                    strand_ = strand[lines[3]]
                except KeyError:
                    strand_ = __determine_strand__(
                        chromosome=lines[0],
                        start=int(lines[1]),
                        end=int(lines[2]),
                        genome_fa=args["genome"]
                    )

                junction = GenomicLoci(
                    chromosome=lines[0],
                    start=int(lines[1]),
                    end=int(lines[2]),
                    strand=strand_
                )

                count = int(lines[6])

                if count < args["single_threshold"]:
                    continue

                tmp[junction] = int(lines[6])

                tmp_threshold[junction] = int(lines[6]) + tmp_threshold.get(junction, 0)

        return {label: tmp}, tmp_threshold

    def read(self, processes):
        u"""
        Read junctions
        :return: {sample label: {junction: count}}
        """

        files = glob.glob(os.path.join(self.input_directory, "*SJ.out.tab"))
        m = Manager()
        lock = m.Lock()
        args = []
        for idx, f in enumerate(files):
            args.append({
                "file_path": f,
                "single_threshold": self.single_threshold,
                "idx": idx,
                "total": len(files),
                "genome": self.genome,
                "lock": lock
            })

        with Pool(processes) as pool:
            res = list(pool.imap_unordered(self.__read_single__, args))

        data = {}
        threshold = set()
        for i, j in res:
            data.update(i)

            for key, value in j.items():
                if value < self.overall_threshold:
                    threshold.add(key)

        return data, threshold

    @staticmethod
    def __write_to_file__(data, output, header=None):
        u"""
        write data to file
        :param data: {infile: {junction: value}}
        :param output: path to output file
        :param header:
        :return:
        """
        samples = sorted(data.keys())
        junctions = set()

        for i in samples:
            for j in data[i]:
                junctions.add(j)

        logger.info("Write to %s" % output)
        with open(output, "w+") as w:
            if header:
                w.write(header + "\t%s\n" % "\t".join(samples))

            junctions = sorted(list(junctions))

            for row in junctions:
                tmp_cols = []
                for col in samples:
                    cell = data[col]

                    if row in cell.keys():
                        tmp_cols.append(str(cell[row]))
                    else:
                        tmp_cols.append("NA")

                w.write("%s\t%s\n" % (row, "\t".join(tmp_cols)))

    @staticmethod
    def __calculate_single__(args):
        u"""

        :param args:
        :return:
        """
        psi_table_start = {}
        psi_table_end = {}

        # args["lock"].acquire()
        # logger.info("%d/%d Calculating %s" % (args["idx"], args["total"], args["key"]))
        # args["lock"].release()

        graph = JunctionGraph(
            process=1,
            alternative_threshold=args["alternative_threshold"]
        )

        graph.add_junctions(args["value"], args["threshold"])

        tmp = graph.get_all_psi()
        tmp_dict = {}
        for k, v in tmp[0].items():
            tmp_dict[k] = v
        psi_table_start[args["key"]] = tmp_dict

        tmp_dict = {}
        for k, v in tmp[1].items():
            tmp_dict[k] = v
        psi_table_end[args["key"]] = tmp_dict

        return psi_table_start, psi_table_end

    def calculate(self):
        u"""
        Calculate PSI
        :return:
        """
        splice_junctions, threshold = self.read(self.processes)

        psi_table_start = {}
        psi_table_end = {}

        m = Manager()
        lock = m.Lock()

        args = []
        for idx, key in enumerate(splice_junctions.keys()):
            args.append({
                "key": key,
                "value": splice_junctions[key],
                "idx": idx + 1,
                "lock": lock,
                "total": len(splice_junctions),
                "threshold": threshold,
                "alternative_threshold": self.alternative_threshold
            })

        with Pool(processes=self.processes) as pool:
            res = list(tqdm(pool.imap(self.__calculate_single__, args), total=len(args)))

        for i, j in res:
            psi_table_start.update(i)
            psi_table_end.update(j)

        self.__write_to_file__(
            data=splice_junctions,
            output=self.output_prefix + ".filtered_counts",
            header="junction"
        )

        # print(psi_table_start)
        self.__write_to_file__(
            data=psi_table_start,
            output=self.output_prefix + ".same_start_psi_table",
            header="junctions"
        )

        self.__write_to_file__(
            data=psi_table_end,
            output=self.output_prefix + ".same_end_psi_table",
            header="junctions"
        )


"""
Main function
"""


def main(args):
    u"""
    Main function
    :param args:
    :return:
    """
    if not args.collect:
        CalPSI(args).calculate()
    else:
        collect_sj_out_tab(
            input_directory=args.input_dir,
            output_file=args.output_prefix,
            genome_fa=args.genome_fasta
        )

    logger.info("Done")


if __name__ == '__main__':

    parser = ap.ArgumentParser(
        """
        Suite of functions for collect SJ counts and calculate PSI \n
        
        WARNING: pysam required for collect -g/--genome-fasta \n

        """,
        formatter_class=ap.ArgumentDefaultsHelpFormatter,

    )

    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='{version}'.format(version="0.1")
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        type=str,
        required=True,
        help="Path to the directory contains SJ.out.tab files"
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        required=True,
        help="The prefix of the output file"
    )

    parser.add_argument(
        "-g",
        "--genome-fasta",
        default=None,
        type=str,
        help="Path to the genome fasta file [optional]"
    )

    parser.add_argument(
        "-c",
        "--collect",
        action="store_true",
        help="If this enabled, only merge all SJ.out.tab into a huge counts table"
    )

    parser.add_argument(
        "-st",
        "--single-threshold",
        type=int,
        default=2,
        help="Filter low abundance junctions in single sample"
    )

    parser.add_argument(
        "-ot",
        "--overall-threshold",
        type=int,
        default=100,
        help="Filter low abundance junctions across all samples"
    )

    parser.add_argument(
        "-at",
        "--alternative-threshold",
        type=int,
        default=5,
        help="Filter low abundance junctions combinations"
    )

    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Use how many CPU, negative number means use all CPU"
    )

    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        args = parser.parse_args(sys.argv[1:])
        print(args)

        main(args)


