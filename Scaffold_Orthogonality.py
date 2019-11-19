#!/usr/bin/env python
# -*- coding: utf-8 -*-3
# 2019 Dietzlab (TUM), [Elija Feigl, Floris Engehardt]
""" This script analyses orthogonality of two DNA-Origami scaffold strands.
    Multiple criteria for orthogonality of the two sequences can be specified
    to determine the level of orthogonality. The script is inteded to be easily
    applicable for a broad untrained audience, weherefore it only relies on
    standard packages and elaborate input processing.
"""

import argparse
import attr
from pathlib import Path
from math import nan

_author__ = "Elija Feigl, Floris Engelhardt"
__copyright__ = "Copyright 2019, Dietzlab (TUM)"
__credits__ = ["TU München"]
__license__ = "None"
__version__ = "0.1"

REVERSE = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
}

# interactive python debugger
# import ipdb; ipdb.set_trace()


@attr.s(slots=True)
class Project(object):
    folder: Path = attr.ib()
    sc1_path: str = attr.ib()
    sc2_path: str = attr.ib()
    n: int = attr.ib()
    is_linear: bool = attr.ib(default=False)
    get_rev_compl: bool = attr.ib(default=False)


def parse_scaffolds(project: dict) -> (str, str):
    def proc_infile(file) -> str:
        sc_str = file.readline().upper().strip("\n")
        is_not_ATGC = not set(sc_str).issubset(set("ATGC"))
        if is_not_ATGC:
            print("faulty sequence")
            exit(0)
        return sc_str

    path_sc1 = project.folder / project.sc1_path
    path_sc2 = project.folder / project.sc2_path

    with open(path_sc1) as sc1_file:
        sc1 = proc_infile(file=sc1_file)
    with open(path_sc2) as sc2_file:
        sc2 = proc_infile(file=sc2_file)
    return sc1, sc2


def check_ortho(sc1: str,
                sc2: str,
                project: Project,
                ) -> (int, int):
    """ Check orthogonality of 2 dna-sequences and output relevant data.
        TODO: (EF) we could try using difflib package to improve performance.
    """
    def complement(segment: str) -> str:
        """ get reverse complementary sequence"""
        for base in segment:
            base = REVERSE[base]
        return segment

    def circularise_sc(project: Project, sc1: str, sc2: str) -> (str, str):
        if project.is_linear:
            sc1_full = sc1
            sc2_full = sc2
        else:
            sc1_full = sc1 + sc1[:project.n]
            sc2_full = sc2 + sc2[:project.n]
        return sc1_full, sc2_full

    def is_match(seg1: str, seg2: str, project: Project) -> bool:
        return seg1 == seg2 and len(seg1) == project.n

    count = 0
    count_RC = 0 if project.get_rev_compl else nan
    n_count, n_count_RC = list(), list()
    sc1_full, sc2_full = circularise_sc(project=project, sc1=sc1, sc2=sc2)

    for i in range(len(sc1)):
        repeat_count, repeat_count_RC = 0, 0
        sc1_seg = sc1_full[i: (i + project.n)]
        for i in range(len(sc2)):
            sc2_seg = sc2_full[i: (i + project.n)]
            if is_match(seg1=sc1_seg, seg2=sc2_seg, project=project):
                count += 1
                repeat_count += 1

            if project.get_rev_compl:
                sc2_seg_RC = complement(sc2_full[(i + project.n):i:-1])
                if is_match(seg1=sc1_seg, seg2=sc2_seg_RC, project=project):
                    count_RC += 1
                    repeat_count_RC += 1

        if repeat_count > 1:
            n_count.append(repeat_count)
        if repeat_count_RC > 1:
            n_count_RC.append(repeat_count_RC)

    count_corrected = count - (sum(n_count) - len(n_count))
    count_revcompl_corrected = (count_RC
                                - (sum(n_count_RC) - len(n_count_RC))
                                )
    output = {"count": count,
              "count_revcompl": count_RC,
              "count_corrected": count_corrected,
              "count_revcompl_corrected": count_revcompl_corrected,
              "n_count": n_count,
              "n_count_revcompl": n_count_RC,
              }
    return output


def write_output(project: dict, data: dict) -> None:
    """ either write uoput to a file or to std-out.
        TODO: (EF) not implemented.
    """
    filename = "ups"
    print("write data to file: ", filename)
    raise NotImplementedError


def get_description() -> str:
    return """evaluate scaffold orthogonality
              """


def proc_input() -> dict:
    parser = argparse.ArgumentParser(
        description=get_description(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-f", "--folder",
                        help="input folder",
                        type=str,
                        default="./",
                        )
    parser.add_argument("-s", "--scaffold1",
                        help="name of scaffold1 file, expects .txt",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        )
    parser.add_argument("-c", "--scaffold2",
                        help="name of scaffold2 file, expects .txt",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        )
    parser.add_argument("-l", "--is_linear",
                        help="scaffolds are not circular",
                        action="store_true"
                        )
    parser.add_argument("-r", "--rev_compl",
                        help="also count reverse complementary sequences",
                        action="store_true"
                        )
    parser.add_argument("-n", "--segment_length",
                        help="segment length",
                        type=int,
                        default=7,
                        )
    args = parser.parse_args()
    project = Project(
        folder=Path(args.folder),
        sc1_path=args.scaffold1,
        sc2_path=args.scaffold2,
        is_linear=args.is_linear,
        get_rev_compl=args.rev_compl,
        n=args.segment_length,
    )
    return project


def main():
    project = proc_input()
    sc1, sc2 = parse_scaffolds(project=project)
    output = check_ortho(
        sc1=sc1,
        sc2=sc2,
        project=project,
    )
    for name, value in output.items():
        print(name, ": ", value)
    # write_output(project=project, data=output)


if __name__ == "__main__":
    main()
