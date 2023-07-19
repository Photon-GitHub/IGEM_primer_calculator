import math
import random

from Bio import Seq, SeqUtils


def is_hairpin(seq: Seq, min_stem_length: int = 2, min_loop_length: int = 3):
    """
    Checks if the given sequence is a hairpin.

    :param seq: The genetic sequence
    :param min_stem_length: The minimum number of base pairs in the stem
    :param min_loop_length: The minimum number of base pairs in the loop

    :return: True if the sequence is a hairpin, False otherwise
    """
    # Check if the sequence is long enough to form a hairpin
    if len(seq) < 2 * min_stem_length + min_loop_length:
        return False

    # Iterate through possible stem lengths
    for stem_length in range(min_stem_length, (len(seq) - min_loop_length) // 2 + 1):
        # Extract the stem and its reverse complement
        stem = seq[:stem_length]
        rev_complement = str(stem.reverse_complement())

        # Check for reverse complement in the remaining sequence
        for reverse_stem_start in range(stem_length + min_loop_length, len(seq) - stem_length + 1):
            if seq[reverse_stem_start:reverse_stem_start + stem_length] == rev_complement:
                return True

    return False


def check_overlap(plasmid: Seq, target_sequence: Seq, target_overlap_len: int, min_stem_length: int, min_loop_length: int) -> Seq:
    seq = target_sequence[0:target_overlap_len]
    return seq if seq not in plasmid and not is_hairpin(seq, min_stem_length, min_loop_length) else None


def calculate_overlap(plasmid: Seq, target_sequence: Seq, min_stem_length: int, min_loop_length: int) -> Seq:
    """Calculate the overlap between the plasmid and the target sequence."""
    for i in range(2, len(target_sequence)):
        seq = check_overlap(plasmid, target_sequence, i, min_stem_length, min_loop_length)
        if seq is not None:
            return seq

    return None


def random_seq(ats: int, cgs: int) -> Seq:
    if ats <= 0 or cgs <= 0:
        raise ValueError("ats and cgs must be greater than 0")

    nucleotides = []

    while ats > 0 or cgs > 0:
        # Choose nucleotide set based on available ATs and CGs
        available_nucleotides = 'ATCG' if ats > 0 and cgs > 0 else 'AT' if ats > 0 else 'CG'

        # Select a random nucleotide and append it to the list
        chosen_nucleotide = random.choice(available_nucleotides)
        nucleotides.append(chosen_nucleotide)

        # Decrement the corresponding count
        if chosen_nucleotide in "AT":
            ats -= 1
        else:
            cgs -= 1

    # Combine the nucleotides into a string and return as a Seq object
    return Seq.Seq(''.join(nucleotides))


def calculate_primer(plasmid: Seq, target_sequence: Seq, primer_length: int, cg_ratio: float, min_stem_length: int, min_loop_length: int) -> (Seq, Seq, Seq):
    """Calculate a primer for a given sequence."""
    overlap = calculate_overlap(plasmid, target_sequence, min_stem_length, min_loop_length)

    if overlap is None:
        print("No overlap possibility found. Are the sequences correct?")
        return None, None, None

    # Calculate the number of CG nucleotides already in the overlap
    cg_in_overlap = sum(1 for nt in overlap if nt in "CGcg")

    # Calculate how many additional CG nucleotides are needed
    total_cg_needed = math.ceil(cg_ratio * primer_length)
    additional_cg_needed = total_cg_needed - cg_in_overlap

    # If additional_cg_needed is negative, set it to 0.
    if additional_cg_needed < 0:
        print("Warning: The overlap already contains too many CG nucleotides to reach the target CG ratio. Consider increasing the primer length.")
        additional_cg_needed = 0

    # Remaining length to fulfill primer_length
    remaining_length = primer_length - len(overlap)

    # Calculate the number of additional random nucleotides needed
    additional_at_needed = remaining_length - additional_cg_needed

    rnd_seq = random_seq(additional_at_needed, additional_cg_needed)
    primer = rnd_seq + overlap

    return rnd_seq, overlap, primer


def calculate_primers(plasmid: Seq, target_sequence: Seq, primer_length: int, cg_ratio: float, min_stem_length: int, min_loop_length: int) -> (Seq, Seq):
    if primer_length < 6:
        print("Primer length must be at least 6")
        return None

    if cg_ratio < 0 or cg_ratio > 1:
        print("CG ratio must be between 0 and 1")
        return None

    target_cropped_len = min(len(target_sequence), 30)

    fw_rnd, fw_overlap, fw_primer = calculate_primer(plasmid, target_sequence, primer_length, cg_ratio, min_stem_length, min_loop_length)
    if fw_rnd is None or fw_overlap is None or fw_primer is None:
        print("No forward primer found")
        return None

    print(f"Found forward primer: {fw_primer} (len: {len(fw_primer)}) with overlap {fw_overlap} (len: {len(fw_overlap)})"
          f" and random part {fw_rnd} (len: {len(fw_rnd)}). CG ratio: {SeqUtils.gc_fraction(fw_primer)}")
    print(" ")
    print(fw_primer)
    print(" " * len(fw_rnd) + str(target_sequence[:target_cropped_len]))
    print(" ")

    rv_plasmid = plasmid[::-1]
    rv_target_sequence = target_sequence[::-1]

    rv_rnd, rv_overlap, rv_primer = calculate_primer(rv_plasmid, rv_target_sequence, primer_length, cg_ratio, min_stem_length, min_loop_length)
    if rv_rnd is None or rv_overlap is None or rv_primer is None:
        print("No reverse primer found")
        return None

    rv_rnd = rv_rnd[::-1]
    rv_overlap = rv_overlap[::-1]
    rv_primer = rv_primer[::-1]

    print(f"Found reverse primer: {rv_primer} (len: {len(rv_primer)}) with overlap {rv_overlap} (len: {len(rv_overlap)})"
          f" and random part {rv_rnd} (len: {len(rv_rnd)}). CG ratio: {SeqUtils.gc_fraction(rv_primer)}")
    print(" ")
    print(str(target_sequence[len(target_sequence) - target_cropped_len:]))
    print(" " * (target_cropped_len - len(rv_overlap)) + str(rv_primer))

    return fw_primer, rv_primer


def verify_primer(plasmid, rest_of_sequence, overlap, min_stem_length, min_loop_length, reverse=False):
    if overlap in plasmid:
        print("The overlap is in the plasmid.")
        return False

    if rest_of_sequence in plasmid:
        print("The rest of the sequence is in the plasmid.")
        return False

    primer = rest_of_sequence + overlap if reverse else overlap + rest_of_sequence
    primer_name = "reverse" if reverse else "forward"
    if is_hairpin(primer, min_stem_length, min_loop_length):
        print(f"The {primer_name} primer is a hairpin.")
        return False

    cg_ratio = SeqUtils.gc_fraction(primer)
    if not 0.4 <= cg_ratio <= 0.6:
        print(f"Warning: The {primer_name} primer has a CG ratio of {cg_ratio}, outside the recommended range of 0.4 to 0.6.")
        print(f"Otherwise, the {primer_name} primer is valid.")
        return False
    else:
        print(f"The {primer_name} primer has a CG ratio of {cg_ratio}.")
        print(f"The {primer_name} primer is valid.")
        return True


def verify_forward_primer(plasmid, rest_of_sequence, overlap, min_hairpin_len=4):
    return verify_primer(plasmid, rest_of_sequence, overlap, min_hairpin_len, False)


def verify_reverse_primer(plasmid, rest_of_sequence, overlap, min_hairpin_len=4):
    return verify_primer(plasmid, rest_of_sequence, overlap, min_hairpin_len, True)


if __name__ == "__main__":
    # Calculate the primer for the given sequence
    #
    # plasmid:          This should be the complete plasmid DNA sequence (Note: This assumes that target_sequence is NOT in the plasmid)
    # target_sequence:  The complete DNA sequence of the target the primer should be designed for
    # primer_length:    The target length of the primer
    # cg_ratio:         The target ratio of C and G in the primer
    # min_stem_length:  The minimum number of base pairs in the stem (the hairpin sequence that is the same in the forward and reverse direction)
    # min_loop_length:  The minimum number of base pairs in the loop (the hairpin sequence that loops around without matching anything)

    plasmid = Seq.Seq(
        "GGAAATTGTAAACGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTGTAATACGACTCACTATAGGGCGAATTGGAGCTCCACCGCGGTGGCGGCCGCTCTAGAACTAGTGGATCCCCCGGGCTGCAGGAATTCGATATCAAGCTTATCGATACCGTCGACCTCGAGGGGGGGCCCGGTACCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATTCCGAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTG")

    target_sequence = Seq.Seq(
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGAGCCACGGCGTCCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTTCAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGCGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAGCCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA")

    fw_primer, rv_primer = calculate_primers(plasmid=plasmid,
                                             target_sequence=target_sequence,
                                             primer_length=20,
                                             cg_ratio=0.55,
                                             min_stem_length=2,
                                             min_loop_length=3)

    # Methods to verify the primers:
    # If you have a primer, identify the overlap, all other nucleotides are the rest of the sequence.
    # min_stem_length:  The minimum number of base pairs in the stem (the hairpin sequence that is the same in the forward and reverse direction)
    # min_loop_length:  The minimum number of base pairs in the loop (the hairpin sequence that loops around without matching anything)
    #
    # verify_forward_primer(plasmid=plasmid, rest_of_sequence=Seq.Seq("PASTE THE REST OF THE SEQUENCE HERE"), overlap=Seq.Seq("PASTE THE OVERLAP WITH THE TARGET SEQUENCE HERE"), min_stem_length=2, min_loop_length=3)
    # verify_reverse_primer(plasmid=plasmid, rest_of_sequence=Seq.Seq("PASTE THE REST OF THE SEQUENCE HERE"), overlap=Seq.Seq("PASTE THE OVERLAP WITH THE TARGET SEQUENCE HERE"), min_stem_length=2, min_loop_length=3)
