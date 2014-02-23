package Genome::Model::Tools::Fasta::Diff;

use strict;
use warnings;

use Genome;
use Command;
use Bio::SeqIO;
use File::Temp;

class Genome::Model::Tools::Fasta::Diff(
    is => 'Command',
    has => [  ## kdiff3 supports 3 files, so we only support 3 files.
        file1 => {
            shell_args_position => 1,
            doc => 'fasta file to diff'
        },
        file2 => {
            shell_args_position => 2,
            doc => 'fasta file to diff'
        },
        file3 => {
            is_optional => 1,
            shell_args_position => 3,
            doc => 'fasta file to diff'
        }
    ]
);

sub help_brief {
    "use KDiff3 to show differences between fasta data files";
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS 
For now, it only works with fasta files with one section
EOS
}


sub execute {
    my $self = shift;

$DB::single = $DB::stopper;
    my @flat_files;
    foreach my $fasta_name ( $self->file1, $self->file2, $self->file3 ) {
        next unless $fasta_name;
        my $fasta = Bio::SeqIO->new(-file => $fasta_name, -format => 'fasta');
        unless ($fasta) {
            $self->error_message("Can't open fasta file $fasta_name");
            return;
        }

        my($fh,$filename) = File::Temp::tempfile;
        my $fseq = $fasta->next_seq();

        $self->debug_message("Flattening $fasta_name...\n");
        for (my $i = 1; $i <= $fseq->length; $i++) {
            $fh->print($fseq->subseq($i,$i),"\n");
        }
        $fh->close();
        push(@flat_files, $filename);
    }

    my $cmdline = 'kdiff3 '. join(' ',@flat_files);
    if(system $cmdline) {
        $self->error_message("Problem running kdiff3");
        return;
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
