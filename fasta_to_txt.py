
from Bio import SeqIO
import argparse


def seq2kmer(seq, k):
    """
    convert original sequence to kmers

    :param
        seq (str): original sequence.
        k (int): kmer of length k specified.

    :returns
        kmers(str): kmers separated by space

    """
    kmer = [seq[x:x + k] for x in range(len(seq) + 1 - k)]
    kmers = " ".join(kmer)

    return kmers


def fasta_to_sequence(input_fasta):
    """
    read fasta file & make fasta list

    Args
        input_fasta (str): input_fasta path

    Returns:
        input_fasta_list (list): fasta information list

    """
    input_fasta_list = []
    for record in SeqIO.parse(input_fasta, 'fasta'):
        input_fasta_list.append(str(record.seq))

    return input_fasta_list


def make_kmer_list(input_fasta_list, k):
    """
    make kmer information list

    Args
        input_fasta_list (list): fasta information list
        k (int): slect kmer unit

    Returns:
        kmer_seq_list (list): kmer information list

    """
    kmer_seq_list = []
    for line in input_fasta_list:
        kmer_seq = seq2kmer(line, k)
        kmer_seq_list.append(kmer_seq)

    return kmer_seq_list


def kmer_list_to_txt(kmer_seq_list, output_file):
    """
    make txt input for BERT training

    Args:
        output_file (str): output file name
        kmer_seq_list (list): kmer information list

    Returns:
        kmer input format txt

    """
    with open(output_file, 'w') as f:
        f.write('\n'.join(kmer_seq_list))

        print('kmer_list_to_txt done')


def main(input_fasta, output_file, k):
    input_fasta_list = fasta_to_sequence(input_fasta)

    kmer_seq_list = make_kmer_list(input_fasta_list, k)

    kmer_list_to_txt(kmer_seq_list, output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--fasta', help='input fasta file')
    parser.add_argument('-o', '--outputfile', help='Output file')
    parser.add_argument('-k', '--kmer', help='kmer unit')
    args = parser.parse_args()

    if args.fasta:
        main(args.fasta, args.outputfile, args.kmer)
