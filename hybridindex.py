#!/usr/bin/env python3
import sys
import random
import argparse
from collections import defaultdict


class Marker():
    def __init__(self, chrom, pos, ref):
        self.chrm = chrom
        self.pos = int(pos)
        self.ref = ref
        self.id = f"{self.chrm}_{pos}"

    def __repr__(self):
        return f"{self.chrm} {self.pos} {self.ref}"

    def make_new(self):
        '''
        Since the main marker is shared across all samples, this allows a keep way to, without requiring the copy module
        make a copy of the object that can then have the alleles in it
        '''
        return Marker(self.chrm, self.pos, self.ref)

    def set_alleles(self, counts):
        self.alleles = counts
        A, T, C, G, N, d = map(int, self.alleles)

        if sum([A, T, C, G, N, d]) == 0:
            self.chosen = "No coverage"
        else:
            choice_str = "A"*A+"T"*T+"C"*C+"G"*G+"N"*N+"D"*d
            chosen_allele = random.choice(choice_str)
            self.chosen = chosen_allele


class Sample():
    def __init__(self, id):
        self.id = id
        self.markers = []
        self.equal_to_ref = 0
        self.diffr_to_ref = 0

    def add_marker(self, string, marker):
        nm = marker.make_new()
        nm.set_alleles(string.split(":"))
        self.markers.append(nm)
        if nm.chosen == nm.ref:
            self.equal_to_ref += 1
        else:
            if nm.chosen == "No coverage":
                pass
            else:
                self.diffr_to_ref += 1

    def __repr__(self):
        return f"Sample {self.id}: With {len(self.markers)} markers"

    def have_choices_in(self):
        for i in self.markers:
            if sum(map(int, i.alleles)) != max(map(int, i.alleles)):
                print(i, i.alleles)


def choose_markers(markers, distance=200000, mode="dynamic"):
    spread_markers = {}
    if mode == "dynamic":

        for chr_marker in markers.values():
            window = [chr_marker[0]]
            current_pos = chr_marker[0].pos
            current_index = 0
            # print(current_pos)
            w_dist = distance//2
            # print("starting to look for windows")
            # spread_markers[chr_marker[0].id] = chr_marker[0]
            t_index = current_index + 1
            while t_index != current_index:
                t_index = current_index
                # print(current_pos)
                chosen = False
                for j, m in enumerate(chr_marker[current_index+1:]):

                    if m.pos - current_pos > w_dist:

                        if not chosen:

                            chosen_marker = random.choice(window)
                            current_pos = chosen_marker.pos
                            spread_markers[chosen_marker.id] = chosen_marker
                            chosen = True

                        if m.pos-current_pos > distance:

                            current_index = current_index+j+1
                            window = [chr_marker[current_index]]
                            current_pos = chr_marker[current_index].pos
                            chosen = False

                            break

                        else:

                            continue
                    else:
                        # print(f"{m} was added to the window")
                        window.append(m)

                # print(window)
                # input()
            chosen_marker = random.choice(window)
            spread_markers[chosen_marker.id] = chosen_marker
        return spread_markers

    elif mode == "static":
        for chr_marker in markers.values():
            spread_markers[chr_marker[0].id] = chr_marker[0]
            current_pos = chr_marker[0].pos

            for m in chr_marker[1:]:
                if m.pos - current_pos > distance:
                    spread_markers[m.id] = m
                    current_pos = m.pos

        return spread_markers

    elif mode == "all":
        for chr_marker in markers.values():
            for m in chr_marker:
                spread_markers[m.id] = m
        return spread_markers
    else:
        raise "Mode is not dynamic,static or all"


if __name__ == "__main__":
    '''
    usage: hybridindex.py [-h] -r replicates -d distance -m mode -ss set_initial_seed -sf sync_file -mf marker_file -of output_file

    Hybrid Index

    Required arguments:
    -sf sync file with the counts
    -mf allele frequency files of the control group
    -of output file

    Optional arguments:
    -r number of replicates (default = 100)
    -d mininum distance, in basepairs, between two markers, if mode = "all" then this argument is ignored (default = 200000)
    -m mode to select the markers, can be dynamic, static, or all. In all, all markers are considered, in static the first marker after d is considered, in dynamic, a window with half distance is used to select the next marker
    '''
    parser = argparse.ArgumentParser(description="Sync parser for hybrid index")
    parser._action_groups = []  # remove all argument groups from the list
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-sf', metavar='sync', type=str, required=True,
                          help='Sync file')
    required.add_argument('-mf', metavar='marker',  type=str, required=True,
                          help="Path to the file with allele frequencies")
    required.add_argument('-of', metavar='output', type=str, required=True,
                          help="Path to the output file")
    optional.add_argument('-r', metavar='n_reps', type=int, default=100,
                          help="Number of repetitions")
    optional.add_argument('-d', metavar='distance', type=int, default=200000,
                          help='mininum distance, in basepairs, between two markers, if mode = "all" then this argument is ignored (default = 200000)')
    optional.add_argument('-m', metavar='mode', type=str, default="dynamic",
                          help="mode to select the markers, can be dynamic, static, or all. In all, all markers are considered. In static the first marker after d is considered. In dynamic, a window with half distance is used to select the next marker")
    optional.add_argument('-ss', metavar='setseed', type=str, default=10000,
                          help="Set seed")

    args = parser.parse_args()
    # random.seed(10001)
    sync_file = args.sf
    marker_file = args.mf
    out_file = args.of
    n_reps = args.r
    distance = args.d
    mode = args.m
    initial_seed = args.ss
    random.seed(initial_seed)
    seeds = [random.randint(1, 100000) for _ in range(n_reps)]
    for rep in range(n_reps):
        print(f"processing rep {rep} of {n_reps}")
        random.seed(seeds[rep])

        markers = defaultdict(list)
        with open(marker_file) as mark:
            for line in mark:
                chrm, pos, ref, alt, prob_dom, prob_wild = line.strip().split("\t")
                prob_dom = float(prob_dom)
                prob_wild = float(prob_wild)
                if prob_dom >= prob_wild:
                    markers[chrm].append(Marker(chrm, pos, ref))
                else:
                    markers[chrm].append(Marker(chrm, pos, alt))

        spread_markers = choose_markers(markers, distance, mode)
        print(f"Number of markers eligible {len(spread_markers)}")

        with open(sync_file) as sync:
            for marker_id, line in enumerate(sync):
                data = line.strip().split("\t")
                chrm, pos, ref = data[0: 3]
                # current_marker = Marker(chrm, pos, ref)
                snp_id = f"{chrm}_{pos}"
                individuals = data[3:]
                if marker_id == 0:  # if reading the first line from the sync file, maker a list with all sample objects
                    samples = [Sample(i) for i in range(len(individuals))]

                elif spread_markers.get(snp_id):

                    for i, allele_counts in enumerate(individuals):

                        samples[i].add_marker(allele_counts, spread_markers[snp_id])

        if rep == 0:
            final_samples = []
            for s in samples:
                final_samples.append([s.id, s.equal_to_ref, s.diffr_to_ref])
        else:
            for i, s in enumerate(samples):
                final_samples[i] += [s.equal_to_ref, s.diffr_to_ref]

    with open(out_file, "w") as out:
        out.write("Sample_ID\t{}\tmean_equal\tmean_different\n".format(
            "\t".join(["rep_{}_equal\trep_{}_diff".format(x, x) for x in range(n_reps)])))
        for sample in final_samples:
            data = sample[1:]
            equal = [x for i, x in enumerate(data) if i % 2 == 0]
            print(equal)
            print(sum(equal))
            print(len(equal))
            equal = sum(equal)/len(equal)
            diffr = [x for i, x in enumerate(data) if i % 2 == 1]
            diffr = sum(diffr)/(len(diffr))
            out.write("{}\t{}\t{}\n".format("\t".join([str(x) for x in sample]), equal, diffr))
