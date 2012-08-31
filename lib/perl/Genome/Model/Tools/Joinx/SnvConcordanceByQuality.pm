package Genome::Model::Tools::Joinx::SnvConcordanceByQuality;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Joinx::SnvConcordanceByQuality {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file_a => {
            is => 'Text',
            is_input => 1,
            doc => 'A sorted bed file containing snvs',
            shell_args_position => 1,
        },
        input_file_b => {
            is => 'Text',
            is_input => 1,
            doc => 'A sorted bed file used to compute concordance (% of input_file_a in input_file_b)',
            shell_args_position => 2,
        },
        output_file => {
            is => 'Text',
            is_input => 1,
            doc => 'The output file',
        },
    ],
};

sub help_brief {
    "Returns concordance information by quality."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt joinx snv-concordance-by-quality a.bed b.bed --output-file report.txt
EOS
}

sub help_detail {
    "For each quality level found in the input bed files (phred scale, 0-255),\n".
    "the number of intersections between A&B and the total number of snvs present\n".
    "at that quality level in A will be computed and displayed."
}

sub execute {
    my $self = shift;
    return unless -e $self->input_file_a and -e $self->input_file_b;

    unless (-s $self->input_file_a and -s $self->input_file_b) {
        IO::File->new($self->output_file, "w") or die "Failed to create file " . $self->output_file;
        return 1;
    }
    my $cmd = $self->joinx_path . ' snv-concordance-by-quality ' . $self->input_file_a . ' ' . $self->input_file_b . ' > ' . $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file_a, $self->input_file_b],
        output_files => [$self->output_file],
    );
    return 1;
}

sub parse_results_file {
    my $file = shift;

    my $fh = IO::File->new("<$file") || die "Failed to open $file";
    my $results = {
        hit => {},
        all => {},
        total_snvs => 0,
        total_hits => 0,
        total_concordance=> 0,
    };
    while (my $l = <$fh>) {
        if ($l =~ /^([0-9]+): ([0-9]+)\/([0-9]+).*/) {
            my ($qual, $hits, $all) = (int($1), int($2), int($3));
            $results->{all}->{$qual} = $all;
            $results->{hit}->{$qual} = $hits;
            my $c = 0;
            $c = $hits*100.0/$all if $all != 0;
            $results->{concordance}->{$qual} = $c;
        } elsif ($l =~ /Overall concordance: ([0-9.]+)/) {
            $results->{total_concordance} = $1;
        } elsif ($l =~ /Total Snvs: ([0-9]+)/) {
            $results->{total_snvs} = $1;
        } elsif ($l =~ /Hits: ([0-9]+)/) {
            $results->{total_hits} = $1;
        }
    }
    $fh->close;
    return $results; 
}

1;
