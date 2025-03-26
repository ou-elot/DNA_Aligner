import numpy
import sys

def read_sequences(filepath):
    """
    Reads a file containing two DNA sequences separated by a newline.

    :param filepath: Path to the file
    :return: Tuple (seq1, seq2)
    """
    with open(filepath, "r") as file:
        lines = file.read().strip().split("\n")
        if len(lines) != 2:
            raise ValueError("File must contain exactly two lines for two DNA sequences.")
        return lines[0], lines[1]  # Returns as (sequence1, sequence2)

def global_alignment(seq1: str, seq2: str, match_reward: int, mismatch_penalty: int, gap_opening: int, gap_extension: int):
  """
  Compute the affine alignment of two strings based on match reward, mismatch penalty, 
  gap opening penalty, and gap extension penalty.
  """
  
