#!/usr/bin/python


def blast_all_vs_all(input_fasta):

    import shell_command

    cmd1 = 'formatdb -i %s -n seq.db' % input_fasta

    a,b,c = shell_command.shell_command(cmd1)
    print (a,b,c)

    cmd2 = 'blastp -db seq.db -query %s -outfmt 6 -num_threads 6 -out blastall.out' % input_fasta

    a,b,c = shell_command.shell_command(cmd2)
    print (a,b,c)

    return "blastall.out"

def run_silix(input_fasta, blastall_file, identity=0.9, overlap=0.8):
    import shell_command

    if identity >1 or overlap>1:
        raise('Identity and overlab should be between 0 and 1')

    cmd = 'silix  %s %s -f FAM -i %s -r %s --net > seq.fnodes' % (input_fasta,
                                                blastall_file,
                                                identity,
                                                overlap)
    a,b,c = shell_command.shell_command(cmd)

def run_hifix(fasta_file,  silix_nodes, silix_net):
    import shell_command

    cmd1 = 'hifix %s %s %s > seq_HFX.fnodes' % (fasta_file,
                                                silix_net,
                                                silix_nodes)

    a,b,c = shell_command.shell_command(cmd1)
    print (a,b,c)


def main(fasta_file):
    blastall_file = blast_all_vs_all(fasta_file)
    #blastall_file = 'blastall.out'
    run_silix(fasta_file, blastall_file, 0.9, 0.8)
    run_hifix(fasta_file, 'seq.fnodes', 'blastall.net')

main('VF_all.fasta')
