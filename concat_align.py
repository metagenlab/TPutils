#! /usr/bin/env python

# concatenation of protein/nucleotide alignments
# sequence concatenated based on fasta header id
# missing taxons replaced by gaps
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 02.2015
# -----------------------------------------------------


def multiple_alignments2concatenated_alignments(fasta_files, out_name):

    from Bio import AlignIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Align import MultipleSeqAlignment
    # identification of all distinct fasta headers id (all unique taxons ids) in all fasta
    # storing records in all_seq_data (dico)
    taxons = []
    all_seq_data = {}
    for one_fasta in fasta_files:
        all_seq_data[one_fasta] = {}
        with open(one_fasta) as f:
            alignment = AlignIO.read(f, "fasta")
        for record in alignment:
            if record.id not in taxons:
                taxons.append(record.id)
            all_seq_data[one_fasta][record.id] = record


    # building dictionnary of the form: dico[one_fasta][one_taxon] = sequence
    concat_data = {}

    for one_fasta in fasta_files:
        for taxon in taxons:
            # check if the considered taxon is present in the record
            if taxon not in all_seq_data[one_fasta]:
                # if taxon absent, create SeqRecord object "-"*len(alignments): gap of the size of the alignment
                seq = Seq("-"*len(all_seq_data[one_fasta][all_seq_data[one_fasta].keys()[0]]))
                all_seq_data[one_fasta][taxon] = SeqRecord(seq, id=taxon)
            if taxon not in concat_data:
                concat_data[taxon] = all_seq_data[one_fasta][taxon]
            else:
                concat_data[taxon] += all_seq_data[one_fasta][taxon]

    # concatenating the alignments, writing to fasta file
    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
    with open(out_name, "w") as handle:
        AlignIO.write(MSA, handle, "fasta")

if __name__ == '__main__':
    import argparse



    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--fasta',type=str,help="input fasta files ", nargs="+")
    parser.add_argument("-o",'--out_name',type=str,help="output msa name (default: msa.fa)", default="msa.fa")
    args = parser.parse_args()

    multiple_alignments2concatenated_alignments(args.fasta, args.out_name)




