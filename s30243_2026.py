# s30243
# 2026-05-06
# program to generate dna sequence files in fasta format
# imlpemnts features
# 1. Batch mode , i.e. the user specifies the number of sequences to generate (e.g. 5); the program generates them in a loop, each with a unique ID (e.g. Seq_001, Seq _002, etc.) and saves them all to one FASTA file ( multi -FASTA).
# 2. Configurable nucleotide distribution, meaning the user enters the percentage of each nucleotide (e.g., A=30, C=20, G=20, T=30). The program validates that the sum equals 100%.
# 3. Searching for motifs .. after generating the sequence, the program searches for the motif specified by the user (e.g. "ATG") and lists all the positions of occurrences (indexing from 1, according to the biological convention).
# 4. In silico transcription - the program generates an mRNA sequence (T to U replacement ) and saves it as a separate record in a FASTA file.

import random
from collections import defaultdict
from typing import Callable

NUCLEOTIDES = ("A", "T", "G", "C")


def generate_sequence(length: int, weights: list[float] = None) -> str:
    """
    Generates a random DNA sequence of a given length.
    If weights are provided, ensures exact nucleotide counts based on percentages.
    """
    if weights is None:
        return "".join(random.choices(NUCLEOTIDES, k=length))

    counts = []
    for w in weights[:-1]:
        counts.append(int(length * w / 100))
    counts.append(length - sum(counts))

    pool = []
    for n, c in zip(NUCLEOTIDES, counts):
        pool.extend([n] * c)

    random.shuffle(pool)
    return "".join(pool)

def calculate_stats(sequence: str) -> dict:
    """
    Calculates nucleotide frequencies and GC-content for a given sequence.
    """
    result = defaultdict(int)
    length = len(sequence)
    result["n"] = length

    for nucleotide in sequence:
        result[nucleotide] += 1

    result["GC"] = result["G"] + result["C"]

    for k, v in result.items():
        result[k] = float(v) / length

    return dict(result)

def print_stats(stats: dict):
    """
    Prints the calculated sequence statistics to the console.
    """
    print("--- Sequence statistic ---")
    print(f"N = {stats['n']}")
    for n in NUCLEOTIDES:
        print(f"{n}: {stats[n]:5.2f}")

    print(f"GC-content: {stats['GC']:5.2f}")

def sequence_lines(sequence: str, line_len: int = 80):
    """
    Splits a sequence string into a list of lines of a specific length.
    """
    lines = []
    for i in range(0, len(sequence), line_len):
        lines.append(sequence[i:i + line_len])
    return lines


def validate_length(value: str, range: tuple[int, int]):
    """
    Validates that a string value is an integer within the specified range.
    """
    try:
        parsed = int(value)
    except Exception as e:
        raise Exception("Integer required.") from e
    if range[0] <= parsed <= range[1]:
        return parsed
    raise Exception("Value must be in range [1, 100 000]")

def validate_id(value: str):
    """
    Validates that a string is a non-empty sequence identifier without internal whitespace.
    """
    stripped = value.strip()
    if not stripped:
        raise Exception("Empty string given")
    if " " in stripped:
        raise Exception("Whitespaces inside id are illegal")

    return stripped

def prompt_valid_input(text: str, validation_function: Callable):
    """
    Repeatedly prompts the user for input until it passes a validation function.
    """
    while True:
        try:
            val = input(f"{text}:\n")
            validated = validation_function(val)
            return validated
        except Exception as e:
            print(f"ERROR: {str(e)}")

def insert_name(sequence: str, name: str) -> str:
    """
    Inserts a lowercase name at a random position within a DNA sequence.
    """
    seq_len = len(sequence)

    rand_pos = random.randint(0, seq_len)
    result = sequence[:rand_pos] + name.lower() + sequence[rand_pos:]
    return result


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """
    Formats sequence data into a standard FASTA string.
    """
    result = f">{seq_id}{' - ' + description if description else ''}\n"
    result += "\n".join(sequence_lines(sequence, line_width))
    return result

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000):
    """
    Prompts for and validates a positive integer within a specified range.
    """
    return prompt_valid_input(prompt, lambda x: validate_length(x, (min_val, max_val)))

def save_file(name: str, content: str):
    """
    Saves a string content to a file with the given name.
    """
    with open(name, "w") as file:
        file.write(content)

def configure_weights() -> list[float]:
    """
    Interactive prompt to configure custom nucleotide distribution percentages.
    """
    while True:
        try:
            print("Enter nucleotide distribution percentages:")
            a = float(input("  A: "))
            t = float(input("  T: "))
            g = float(input("  G: "))
            c = float(input("  C: "))
            weights = [a, t, g, c]
            if abs(sum(weights) - 100.0) < 0.001:
                return weights
            print(f"ERROR: Sum must be 100% (currently {sum(weights)}%)")
        except ValueError:
            print("ERROR: Invalid input. Please enter numbers.")

def batch_mode():
    """
    Executes the program in batch mode, generating multiple sequences into a multi-FASTA file.
    """
    count = validate_positive_int("Number of sequences to generate")
    length = validate_positive_int("Length of each sequence")

    weights = None
    if prompt_valid_input("Use custom distribution? (y/n)", str.lower) == "y":
        weights = configure_weights()

    filename = prompt_valid_input("Output filename (e.g., batch.fasta)", str.strip)
    if not filename:
        filename = "batch.fasta"

    batch_content = []
    for i in range(1, count + 1):
        seq_id = f"Seq_{i:03d}"
        sequence = generate_sequence(length, weights=weights)
        content = format_fasta(seq_id, "", sequence)
        batch_content.append(content)

    save_file(filename, "\n".join(batch_content))
    print(f"Generated {count} sequences to {filename}")


def transcribe(sequence: str) -> str:
    """
    Performs in silico transcription by replacing Thymine (T) with Uracil (U).
    """
    return sequence.replace("T", "U")

def search_motif(sequence: str, motif: str) -> list[int]:
    """
    Searches for a specific motif in a sequence and returns all 1-based start positions.
    """
    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)  # 1-based indexing
        start = pos + 1
    return positions

def single_mode():
    """
    Executes the program in single mode, generating one sequence with custom options.
    """
    length = validate_positive_int("Sequence length")

    weights = None
    if prompt_valid_input("Use custom distribution? (y/n)", str.lower) == "y":
        weights = configure_weights()

    seq_id = prompt_valid_input("Sequence identifier", validate_id)
    description = prompt_valid_input("Description", str.strip)
    name = prompt_valid_input("Input user's name", str.strip)

    sequence = generate_sequence(length, weights=weights)
    sequence = insert_name(sequence, name)
    content = format_fasta(seq_id, description, sequence)

    if prompt_valid_input("Transcribe to mRNA? (y/n)", str.lower) == "y":
        rna_seq = transcribe(sequence)
        content += "\n" + format_fasta(f"{seq_id}_mRNA", f"mRNA transcript of {seq_id}", rna_seq)

    motif = prompt_valid_input("Enter motif to search (or press Enter to skip)", str.strip).upper()
    if motif:
        positions = search_motif(sequence, motif)
        if positions:
            print(f"Motif '{motif}' found at positions: {', '.join(map(str, positions))}")
        else:
            print(f"Motif '{motif}' not found.")

    stats = calculate_stats(sequence)
    save_file(f"{seq_id}.fasta", content)
    print_stats(stats)


def main():
    """
    Main entry point of the script. Handles mode selection.
    """
    mode = prompt_valid_input("Choose mode: [s]ingle or [b]atch", lambda x: x.lower() if x.lower() in ["s", "b"] else exec('raise Exception("Invalid mode")'))

    if mode == "s":
        single_mode()
    else:
        batch_mode()


if __name__ == '__main__':
    main()
