import random
import sys
from collections import defaultdict

NUCLEOTIDES = ("A", "T", "G", "C")


def generate_sequence(length: int) -> str:
    return "".join(random.choices(NUCLEOTIDES, k=length))

def calculate_stats(sequence: str) -> dict:
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
    print(f"--- Sequence statistic ---")
    print(f"N = {stats["n"]}")
    for n in NUCLEOTIDES:
        print(f"{n}: {stats[n]:5.2f}")

    print(f"GC-content: {stats["GC"]:5.2f}")

def sequence_lines(sequence: str, line_len: int = 80):
    lines = []
    for i in range(0, len(sequence), line_len):
        lines.append(sequence[i:i + line_len])
    return lines


def validate_length(value: str, range: tuple[int, int]):
    try:
        parsed = int(value)
    except:
        raise Exception("Integer required.")
    if range[0] <= parsed <= range[1]:
        return parsed
    raise Exception("Value must be in range [1, 100 000]")

def validate_id(value: str):
    stripped = value.strip()
    if not stripped:
        raise Exception("Empty string given")
    if " " in stripped:
        raise Exception("Whitespaces inside id are illegal")

    return stripped

def prompt_valid_input(text: str, validation_function: callable):
    while True:
        try:
            val = input(f"{text}:\n")
            validated = validation_function(val)
            return validated
        except Exception as e:
            print(f"ERROR: {str(e)}")

def insert_name(sequence: str, name: str) -> str:
    seq_len = len(sequence)

    rand_pos = random.randint(0, seq_len)
    result = sequence[:rand_pos] + name.lower() + sequence[rand_pos:]
    return result


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    result = f">{seq_id}{" - " + description if description else ""}\n"
    result += "\n".join(sequence_lines(sequence, line_width))
    return result

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000):
    return prompt_valid_input(prompt, lambda x: validate_length(x, (min_val, max_val)))

def save_file(name: str, content: str):
    with open(name, "w") as file:
        file.write(content)

def main():
    length = validate_positive_int("Sequence length")
    id = prompt_valid_input("Sequence identifier", validate_id)
    description = prompt_valid_input("Description", str.strip)
    name = prompt_valid_input("Input name", str.strip)

    sequence = generate_sequence(length)
    sequence = insert_name(sequence, name)
    content = format_fasta(id, description, sequence)
    stats = calculate_stats(sequence)
    save_file(f"{id}.fasta", content)
    print_stats(stats)





if __name__ == '__main__':
    main()

