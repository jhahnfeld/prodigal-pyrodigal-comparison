#!/usr/bin/env python3

import os
import csv
import subprocess as sp
from pathlib import Path
from collections import namedtuple
from argparse import ArgumentParser

import pyrodigal
from xopen import xopen
from Bio import SeqIO


def prodigal_train(genome: Path, train_file: Path, prodigal_bin: Path, transl_table: int = 11):
    """
    Create a Prodigal training file.
    :param genome: Path to the genome
    :param train_file: Path to the training file
    :param prodigal_bin: Path to the prodigal binary
    :param transl_table: Translation table (11)
    """
    cmd = [
        prodigal_bin,
        '-i', str(genome),
        '-g', str(transl_table),
        '-t', str(train_file),
        '-c'
    ]
    proc = sp.run(
        cmd,
        env=os.environ.copy(),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        print('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        raise Exception(f'prodigal error! error code: {proc.returncode}')


def pyrodigal_train(genome: Path, transl_table: int = 11) -> pyrodigal.OrfFinder:
    """
    Train Pyrodigal
    :param genome: Path to the genome
    :param transl_table: Translation table (11)
    :return: pyrodigal.OrfFinder instance
    """
    sequences = []
    with xopen(genome, mode='r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            sequences.append(str(record.seq))

    orf_finder = pyrodigal.OrfFinder(closed=True)
    orf_finder.train(*sequences, translation_table=transl_table)

    return orf_finder


def prodigal_predict(genome: Path, train_file: Path, gff_file: Path, prodigal_bin: Path, transl_table: int = 11):
    """
    Predict CDS with Prodigal
    :param genome: Path to the genome (fasta)
    :param train_file: Path to the training file
    :param gff_file: Path to GFF output file
    :param prodigal_bin: Path to the prodigal binary
    :param transl_table: Translation table
    """
    cmd = [
        prodigal_bin,
        '-i', str(genome),
        '-g', str(transl_table),
        '-t', str(train_file),
        '-c',
        '-f', 'gff',
        '-o', str(gff_file)
    ]
    proc = sp.run(
        cmd,
        env=os.environ.copy(),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        print('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        raise Exception(f'prodigal error! error code: {proc.returncode}')


def parse_gff_output(gff_path: Path, pyrodigal_processed: bool = False) -> set[namedtuple]:
    """
    Parse the Prodigal output.
    :param gff_path: Path to the GFF output file
    :param pyrodigal_processed: Pyrodigal output
    :return: Set of namedtuples of CDS coordinates
    """
    cdss: set[namedtuple] = set()
    Gene: namedtuple = namedtuple('Gene', ['id',
                                           'start',
                                           'stop',
                                           'strand',
                                           'edge'])
    with gff_path.open() as fh:
        for line in fh:
            if line[0] != '#':
                (contig_id, inference, _, start, stop, score, strand, _, annotations_raw) = line.strip().split('\t')
                gff_annotations = split_gff_annotation(annotations_raw)

                if pyrodigal_processed:
                    c = contig_id.count('_')
                    contig_id = '_'.join(contig_id.split('_')[0:c])

                strand = 1 if strand == '+' else -1
                start, stop = int(start), int(stop)

                gene: namedtuple = Gene(contig_id,
                                        int(start),
                                        int(stop),
                                        strand,
                                        gff_annotations['partial'].startswith('1') or gff_annotations['partial'].endswith('1'))
                cdss.add(gene)
    return cdss


def split_gff_annotation(annotation_string: str) -> dict[str, str]:
    """
    Split the annotation string and save its content in a dictionary.
    :param annotation_string: Prodigal GFF annotation string
    :return: Splitted annotations
    """
    annotations = {}
    for expr in annotation_string.split(';'):
        if expr != '':
            key, value = expr.split('=')
            annotations[key] = value
    return annotations


def pyrodigal_predict(orf_finder: pyrodigal.OrfFinder, genome: Path) -> set[namedtuple]:
    """
    Predict all genes using Pyrodigal.
    :param orf_finder: Trained Pyrodigal OrfFinder instance
    :param genome: Path to the contig (fasta)
    :return:
    """
    cdss: set[namedtuple] = set()
    Gene: namedtuple = namedtuple('Gene', ['id',
                                           'start',
                                           'stop',
                                           'strand',
                                           'edge'])
    with xopen(genome, mode='r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            for prediction in orf_finder.find_genes(bytes(record.seq)):
                gene: namedtuple = Gene(record.id,
                                        prediction.begin,
                                        prediction.end,
                                        prediction.strand,
                                        prediction.partial_begin or prediction.partial_end)
                cdss.add(gene)
    return cdss


def compare_mismatches(prodigal_genes: set[namedtuple], pyrodigal_genes: set[namedtuple]) -> list:
    """
    Compare prodigal/pyrodigal mismatches to find related hits.
    :param prodigal_genes: Mismatched Prodigal genes
    :param pyrodigal_genes: Mismatched Pyrodigal genes
    :return: Related hits
    """
    ordered_pairs: list = []

    Gene: namedtuple = namedtuple('Gene', ['id',
                                           'start',
                                           'stop',
                                           'strand',
                                           'edge'])
    empty_gene: namedtuple = Gene(None,
                                  None,
                                  None,
                                  None,
                                  None)

    longer, shorter = (prodigal_genes, pyrodigal_genes) if len(prodigal_genes) >= len(pyrodigal_genes) else (pyrodigal_genes, prodigal_genes)
    prodigal_more_hits = True if len(prodigal_genes) >= len(pyrodigal_genes) else False

    for main_cds in longer:
        for candidate_cds in shorter:
            if main_cds.id == candidate_cds.id and main_cds.strand == candidate_cds.strand:  # matching contig AND strand
                if main_cds.stop == candidate_cds.stop:  # matching stop
                    mismatch = {'prodigal': main_cds,
                                'pyrodigal': candidate_cds,
                                'strand': main_cds.strand} if prodigal_more_hits else {'prodigal': candidate_cds,
                                                                                       'pyrodigal': main_cds,
                                                                                       'strand': main_cds.strand}
                    ordered_pairs.append(mismatch)
                    break
        else:
            mismatch = {'prodigal': main_cds,
                        'pyrodigal': empty_gene,
                        'strand': main_cds.strand} if prodigal_more_hits else {'prodigal': empty_gene,
                                                                               'pyrodigal': main_cds,
                                                                               'strand': main_cds.strand}
            ordered_pairs.append(mismatch)

    return ordered_pairs


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Compare CDS predictions of Prodigal to the predictions of Pyrodigal and save the differences in a TSV file.')
    parser.add_argument('--genome', '-g', type=Path, nargs='+', action='append', help='Input genomes (/some/path/*.fasta)')
    parser.add_argument('--prodigal', '-p', type=Path, default=Path('/vol/sorf/bakta/Prodigal/bin/prodigal'),
                        help='Path to a newly compiled Prodigal binary.')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output path (default="./comparison")')
    args = parser.parse_args()
    return args.genome, args.prodigal, args.output


def main():
    genomes, prodigal_bin, out_path = parse_arguments()
    out_path.resolve()

    tsv_file = open(out_path.joinpath('mismatches.tsv'), mode='w')
    tsv_writer = csv.writer(tsv_file, delimiter='\t', lineterminator='\n')
    header: list[str] = ['genome', 'strand', 'prodigal_start', 'prodigal_stop', 'prodigal_edge', 'pyrodigal_start', 'pyrodigal_stop', 'pyrodigal_edge']
    tsv_writer.writerow(header)

    for genome in genomes[0]:
        genome.resolve()
        prefix: str = str(genome).split('/')[-1].split('.')[0]

        tmp_path = out_path.joinpath(f'tmp/{prefix}')
        try:
            os.makedirs(tmp_path)
        except FileExistsError:
            pass

        # Train prodigal
        prodigal_train_file = tmp_path.joinpath(f'{prefix}.prodigal.trn')
        prodigal_train(genome=genome,
                       train_file=prodigal_train_file,
                       prodigal_bin=prodigal_bin)
        # Train pyrodigal
        orf_finder = pyrodigal_train(genome)

        # Predict prodigal
        prodigal_gff_file = tmp_path.joinpath(f'{prefix}.prodigal.gff')
        prodigal_predict(genome=genome,
                         train_file=prodigal_train_file,
                         gff_file=prodigal_gff_file,
                         prodigal_bin=prodigal_bin)
        prodigal_genes: set[dict] = parse_gff_output(gff_path=prodigal_gff_file)
        # Predict pyrodigal
        pyrodigal_genes: set[dict] = pyrodigal_predict(orf_finder=orf_finder,
                                                       genome=genome)

        print(f'Hits genome={prefix}: prodigal={len(prodigal_genes)}, pyrodigal={len(pyrodigal_genes)}, equal={prodigal_genes == pyrodigal_genes}')

        if prodigal_genes != pyrodigal_genes:
            difference_prodigal = prodigal_genes.difference(pyrodigal_genes)
            difference_pyrodigal = pyrodigal_genes.difference(prodigal_genes)
            ordered_pairs = compare_mismatches(difference_prodigal, difference_pyrodigal)

            for mismatch in ordered_pairs:
                # genome, strand, prodigal_start, prodigal_stop, prodigal_edge, pyrodigal_start, pyrodigal_stop, pyrodigal_edge
                line = [prefix,
                        mismatch['strand'],
                        mismatch['prodigal'].start,
                        mismatch['prodigal'].stop,
                        mismatch['prodigal'].edge,
                        mismatch['pyrodigal'].start,
                        mismatch['pyrodigal'].stop,
                        mismatch['pyrodigal'].edge]

                tsv_writer.writerow(line)

    tsv_file.close()


if __name__ == '__main__':
    main()


# EOF
