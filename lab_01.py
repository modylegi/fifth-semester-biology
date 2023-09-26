import random


def generate_dna_string(length, gc_content):
    num_gc = int(length * gc_content / 100)
    num_at = length - num_gc
    dna_string = 'G' * (num_gc // 2) + 'C' * (num_gc // 2) + 'A' * (num_at // 2) + 'T' * (num_at // 2)
    dna_list = list(dna_string)
    random.shuffle(dna_list)
    return ''.join(dna_list)


# ... (rest of the code remains the same)


def reverse_complement(dna_string):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(dna_string))


def find_orfs(dna_string):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    for frame in range(3):
        for i in range(frame, len(dna_string) - 2, 3):
            codon = dna_string[i:i + 3]
            if codon == start_codon:
                end = i + 3
                while end < len(dna_string) - 2:
                    if dna_string[end:end + 3] in stop_codons:
                        if (end - i) >= 30:  # Minimum length of 10 amino acids
                            orf_dna = dna_string[i:end + 3]
                            orf_protein = translate(orf_dna)
                            orfs.append((orf_dna, orf_protein, i + 1, end + 3, frame + 1))
                        break
                    end += 3
    return orfs


def translate(dna_string):
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    for i in range(0, len(dna_string), 3):
        codon = dna_string[i:i + 3]
        protein += codon_table.get(codon, '_')
    if protein[-1] == "_":
        protein = protein[:-1] + "*"
    return protein


# Input validation
while True:
    try:
        length = int(input("Enter the length of the DNA string (between 100 and 1000): "))
        if length < 100 or length > 1000:
            raise ValueError("Length out of range")
        gc_content = int(input("Enter the GC-content (between 20 and 80): "))
        if gc_content < 20 or gc_content > 80:
            raise ValueError("GC-content out of range")
        break
    except ValueError as e:
        print(e)

dna_string = generate_dna_string(length, gc_content)
reverse_dna_string = reverse_complement(dna_string)

print("\nGenerated DNA String:")
print(dna_string)

orfs = find_orfs(dna_string)
reverse_orfs = find_orfs(reverse_dna_string)

if not orfs and not reverse_orfs:
    print("\nNo ORFs found.")
else:
    print("\nORFs found:")
    for orf in orfs:
        print(f"\nORF (Forward Strand): {orf[0]}")
        print(f"Translated Protein: {orf[1]}")
        print(f"DNA Strand Range: {orf[2]}-{orf[3]}")
        print(f"Strand: Forward")
        print(f"Reading Frame: {orf[4]}")

    for orf in reverse_orfs:
        print(f"\nORF (Reverse Strand): {orf[0]}")
        print(f"Translated Protein: {orf[1]}")
        print(f"DNA Strand Range: {orf[2]}-{orf[3]}")
        print(f"Strand: Reverse")
        print(f"Reading Frame: {orf[4]}")
