package Genome::Model::Tools::Annotate::Adaptor::IndelContig;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::Annotate::Adaptor::IndelContig {
    is => 'Command',
    has => [
    contig_fasta_file => {
        type => 'String',
        is_optional => 0,
        doc => 'File of contigs from gmt validation build-remapping-contigs.',
        default => '',
    },
    contig_count_file => {
        type => 'String',
        is_optional => 0,
        doc => 'File of sites to convert from gmt validation count-contigs or gmt validation combine-counts.',
    },
    output_file => {
        is => 'Text',
        is_input => 1,
        is_output => 1,
        is_optional => 1,
        doc => "Store output in the specified file instead of sending it to STDOUT."
    },
    ]
};

sub execute {
    my $self=shift;
    my $contigs = $self->contig_fasta_file;
    my $sites = $self->contig_count_file;

    my $ofh;
    if($self->output_file) {
        $ofh = IO::File->new($self->output_file,"w");
    }
    else {
        $ofh = IO::File->new_from_fd(fileno(STDOUT),"w");
    }
    unless($ofh) {
        $self->error_message("Unable to open handle for output");
    }

    my $fh = IO::File->new($sites);
    unless($fh) {
        $self->error_message("Couldn't open site list " . $sites);
        return;
    }

    my %sites_to_report;
    while(my $line = $fh->getline) {
        chomp $line;
        $line =~ s/"//g;    #remove quoting. Shouldn't cause big problems. I hope.
        my @fields = split /\t/, $line;
        $sites_to_report{$fields[0]} = 1;   #add contig id to set of ids to parse out
    }
    $fh->close;

    my $cfh = IO::File->new($contigs);
    unless($cfh) {
        $self->error_message("Couldn't open fasta of contigs " . $contigs);
        return;
    }

    while(my $line = $cfh->getline) {
        next unless $line =~ /^>/;
        chomp $line;
        my @fields = split /\s+/,$line; #this should break down the fields
        next if $fields[0] =~ /contig/i;
        $fields[0] =~ s/^>//;   #remove the leading bracket
        if(exists($sites_to_report{$fields[0]})) {
            my ($ref_count_chr, $ref_count_start, $ref_count_stop,$ref,$var,$type) = $fields[4] =~ /Anno:([^.]+)[.]([0-9]+)[.]([0-9]+)[.](\S+)[.](\S+)[.](\S+)/;
            if($ref_count_start > $ref_count_stop) {
                warn "Reference coordinates to count make no sense. Swapping start and stop.";
                print STDERR $line,"\n";
                ($ref_count_start,$ref_count_stop) = ($ref_count_stop, $ref_count_start);
            }
            print $ofh "$ref_count_chr\t$ref_count_start\t$ref_count_stop\t$ref\t$var\t$type\t$fields[0]\n";
        }
    }

    $cfh->close;

    return 1;
}


1;

sub help_brief {
    "Converts contigs from gmt validation build-remapping-contigs output into a gmt annotate transcript-variants friendly input format",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt annotate adaptor indel-contig ...    
EOS
}

sub help_detail {
    <<'HELP';
Converts contigs from gmt validation build-remapping-contigs output into a gmt annotate transcript-variants friendly input format
HELP
}
