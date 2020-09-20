"""
##############################################################################
#
#   Merge short sequence motifs into a bigger consensus sequences.
#   Permutes the motif list 100x in order to escape dependancy on the order.
#   Optimal clustering is the one with the highest standard deviation
#   in the pairwise distances between the cardinalities of the clusters.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 20-09-2020
#   LICENSE: GPL_v3.0
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np
import random
import copy


def parse_arguments():
    """Parser of the command-line arguments."""
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-v",
        "--verbosity",
        dest="verbosity",
        choices=("DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"),
        default="ERROR",
        help="Verbosity/Log level. Defaults to ERROR",
    )
    parser.add_argument(
        "-l", "--logfile", dest="logfile", help="Store log to this file."
    )
    parser.add_argument(
        "--in",
        dest="infile",
        required=True,
        help="Path to the textfile with sequence motifs.",
    )
    parser.add_argument(
        "--overlap",
        dest="overlap",
        required=True,
        help="Minimal number of bases to overlap for the merge.",
    )
    return parser


##############################################################################


class MergedMotif:
    def __init__(self, consensus, motifs):
        self.consensus = consensus
        self.motifs = motifs


def slide(arg1, arg2):
    """Slide one string over the other and calculate the match scores."""

    s1 = arg1.consensus
    s2 = arg2.consensus
    # prepare the alignment strings:
    len1 = len(s1)
    len2 = len(s2)
    for i in range(len2 - 1):
        s1 = "-" + s1 + "-"
    for i in range(len1 + len2 - 2):
        s2 = s2 + "-"
    aln_len = len(s1)
    match = False

    for i in range(1, len1 + len2):
        # compute the alignment score:
        score = 0
        for pos in range(aln_len):
            if s1[pos] != "-" and s2[pos] != "-":
                if s1[pos] == s2[pos]:
                    score += 1  # match
                else:
                    score -= 100  # mismatch

        # if we find a sufficient match finish now
        # else: scan further
        if score >= min(int(options.overlap), len1, len2):
            match = True
            break
        else:
            s2 = "-" + s2[: len(s2) - 1]

    if match:
        # prune the "-" from both ends
        for pos in range(-1, -aln_len - 1, -1):
            if not (s1[pos] == "-" and s2[pos] == "-"):
                break
        s1 = s1[: pos + 1]
        s2 = s2[: pos + 1]
        aln_len = len(s1)
        for pos in range(aln_len):
            if not (s1[pos] == "-" and s2[pos] == "-"):
                break
        s1 = s1[pos:]
        s2 = s2[pos:]
        aln_len = len(s1)

        # merge into the consesnsus sequence
        consensus = ""
        for pos in range(aln_len):
            if s1[pos] == "-":
                consensus = consensus + s2[pos]
            else:
                consensus = consensus + s1[pos]
        return (score, MergedMotif(consensus, arg1.motifs + arg2.motifs))
    else:
        return (-1, None, None, None)


def cluster_motifs(S):
    """Cluster a set of motifs based on shared subsequences"""

    loop = True

    while loop and len(S) > 1:

        loop = False

        best_match = (None, None, 0, None)  # s1,s2,score,merged
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                score_match = slide(S[i], S[j])
                if score_match[0] > best_match[2]:
                    best_match = (S[i], S[j], score_match[0], score_match[1])

        # if the best match is still to be merged - merge it and loop again
        if best_match[2] > 0:
            S.remove(best_match[0])
            S.remove(best_match[1])
            S.append(best_match[3])
            loop = True

    # if we converged to one cluster then this is an optimal solution
    if len(S) == 1:
        return (S, len(S[0].motifs), 0)

    else:
        # get the cardinality of the clusters
        clusters_cardinality = tuple([len(i.motifs) for i in S])
        clusters_distances = []
        for i in range(len(clusters_cardinality)):
            for j in range(i + 1, len(clusters_cardinality)):
                clusters_distances.append(
                    abs(clusters_cardinality[i] - clusters_cardinality[j])
                )
        clusters_distances = tuple(clusters_distances)

        # print(clusters_cardinality)
        # print(clusters_distances)
        # print(np.std(clusters_distances))

        # clustering, number of clusters, stdev(distances)
        return (S, len(clusters_cardinality), np.std(clusters_distances))


def main():
    """Main body of the script."""

    # read the list of motifs and generate the initial set
    with open(options.infile, "r") as f:
        motifs = f.read().splitlines()
    S_init = [MergedMotif(m, [m]) for m in motifs]

    # permute randomly 1000x times and cluster the motifs
    clusterings = []
    for _ in range(1000):
        S = copy.copy(S_init)
        random.shuffle(S)
        clusterings.append(cluster_motifs(S))

    # first sieve: less clusters is better
    min_number_of_clusters = min([i[1] for i in clusterings])
    clusterings = list(filter(lambda x: x[1] == min_number_of_clusters, clusterings))

    # second sieve: higher std is better
    max_std = max([i[2] for i in clusterings])
    clusterings = list(filter(lambda x: x[2] == max_std, clusterings))
    # if there is >1 then these are the same clusterings that originate from different
    # order of the initial list - does not matter, pick the first one
    S = clusterings[0][0]

    # print the cluster, consensus and alignments
    for i in S:
        print(i.motifs)
        print(i.consensus)
        if len(i.motifs) == 1:
            print(i.consensus)
        else:
            # print in the order of alignment: position of the 1st base
            m_to_print = []
            for m in i.motifs:
                pos = i.consensus.find(m)
                aln_string = ""
                for xx in range(pos):
                    aln_string += "-"
                aln_string += m
                suffix = len(i.consensus) - len(aln_string)
                for xx in range(suffix):
                    aln_string += "-"
                m_to_print.append(aln_string)
            for m in sorted(m_to_print)[::-1]:
                print(m)
        print("")


##############################################################################

if __name__ == "__main__":

    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()

        # set up logging during the execution
        formatter = logging.Formatter(
            fmt="[%(asctime)s] %(levelname)s - %(message)s",
            datefmt="%d-%b-%Y %H:%M:%S",
        )
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger = logging.getLogger("logger")
        logger.setLevel(logging.getLevelName(options.verbosity))
        logger.addHandler(console_handler)
        if options.logfile is not None:
            logfile_handler = logging.handlers.RotatingFileHandler(
                options.logfile, maxBytes=50000, backupCount=2
            )
            logfile_handler.setFormatter(formatter)
            logger.addHandler(logfile_handler)

        # execute the body of the script
        start_time = time.time()
        logger.info("Starting script")
        main()
        seconds = time.time() - start_time

        # log the execution time
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        logger.info(
            "Successfully finished in {hours}h:{minutes}m:{seconds}s",
            hours=int(hours),
            minutes=int(minutes),
            seconds=int(seconds) if seconds > 1.0 else 1,
        )
    # log the exception in case it happens
    except Exception as e:
        logger.exception(str(e))
        raise e
