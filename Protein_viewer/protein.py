from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PPBuilder
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqIO.PdbIO import _res2aacode

class Protein:
    L_MIN, L_MAX = 1, 30
    L_RANGE = range(L_MIN, L_MAX + 1)

    def __init__(self, path):
        self.parser = PDBParser()
        self.path = Path(path)
        self.structure = self.parser.get_structure(self.path.stem, self.path)
        self._chains = None
        self._kmers = None

    @property
    def chains(self):
        if self._chains is None:
            self._chains = self.get_chains()
        return self._chains

    def get_chains(self):
        return {chain.id: chain for model in self.structure for chain in model}

    @property
    def kmers(self):
        if self._kmers is None:
            self._kmers = self.get_kmers()
        return self._kmers

    @property
    def methods(self):
        return ["Euclidean distance"]

    def get_kmers(self):
        kmers_dict = {}
        for chain_id, chain in self.chains.items():
            residues = [r for r in chain.get_residues() if r.get_id()[0] == " "]
            chain_kmers_dict = defaultdict(list)

            for i, _ in enumerate(residues):
                for kmer_len in self.L_RANGE:
                    kmers = []
                    if i + kmer_len > len(residues):
                        break
                    for res in residues[i : i + kmer_len]:
                        kmers.append(res)
                    chain_kmers_dict[kmer_len].append(kmers)
            kmers_dict[chain_id] = chain_kmers_dict
        return kmers_dict

    def get_kmer_pairs(self, chain_id_1, chain_id_2):
        kmer_pairs = defaultdict(list)
        if chain_id_1 not in self.kmers.keys() or chain_id_2 not in self.kmers.keys():
            return f"Wrong chain id, available chains: {self.kmers.keys()}"

        for chain_kmer_len in self.kmers[chain_id_1]:
            for chain_kmer_1 in self.kmers[chain_id_1][chain_kmer_len]:
                for chain_kmer_2 in self.kmers[chain_id_2][chain_kmer_len]:
                    kmer_pairs[chain_kmer_len].append((chain_kmer_1, chain_kmer_2))

        return kmer_pairs

    def get_seq_from_residues(self, residues):
        seq = ""
        for res in residues:
            amino_acid_letter = _res2aacode(res)
            seq += amino_acid_letter
        return seq

    @staticmethod
    def get_distance_between_residues(res_1, res_2):
        if "CA" not in res_1 or "CA" not in res_2:
            return float("NaN")
        coord_1 = res_1["CA"].get_coord()
        coord_2 = res_2["CA"].get_coord()

        return np.linalg.norm(coord_1 - coord_2) #return euclidian distance between coords

    def get_kmer_strings_with_distances(self, chain_id_1, chain_id_2, kmer_len):
        kmer_pairs = self.get_kmer_pairs(chain_id_1, chain_id_2)
        kmer_string_with_distances = defaultdict(list)

        for kmer_pair in kmer_pairs[kmer_len]:
            first_kmer_str = self.get_seq_from_residues(kmer_pair[0])
            second_kmer_str = self.get_seq_from_residues(kmer_pair[1])
            residue_pairs = [(res, kmer_pair[1][i]) for i, res in enumerate(kmer_pair[0])]
            distances = [self.get_distance_between_residues(res_1, res_2) for res_1, res_2 in residue_pairs]
            kmer_string_with_distances[kmer_len].append(
                (
                    first_kmer_str,
                    second_kmer_str,
                    distances,
                    np.mean(distances).round(2),
                )
            )

        return kmer_string_with_distances[kmer_len]

    def get_alpha(self):
        ppb = PPBuilder()
        v = dict()
        for chain_id in self.chains:
            vectors = []
            pp = ppb.build_peptides(self.chains[chain_id])
            for el in pp:
                alpha_carbons = el.get_ca_list()
                for atom in alpha_carbons:
                    residue = atom.get_parent()
                    residue_id = residue.get_full_id()
                    resseq = residue_id[3][1]
                    vectors.append([*atom.get_vector(), resseq])
            vectors = np.array(vectors)
            v[chain_id] = vectors

        return v