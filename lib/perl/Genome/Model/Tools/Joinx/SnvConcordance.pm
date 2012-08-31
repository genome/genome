package Genome::Model::Tools::Joinx::SnvConcordance;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Joinx::SnvConcordance {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file_a => {
            is => 'Text',
            doc => 'A sorted bed file containing snvs',
            shell_args_position => 1,
        },
        input_file_b => {
            is => 'Text',
            doc => 'A sorted bed file used to compute concordance (% of input_file_a in input_file_b)',
            shell_args_position => 2,
        },
    ],
    has_optional_input => [
        depth => {
            is => 'Boolean',
            doc => 'Use read depth instead of quality in report',
        },
        output_file => {
            is => 'Text',
            doc => 'The output file (defaults to stdout)',
        },
    ],
};

sub help_brief {
    "Returns concordance information for 2 snv lists."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt joinx snv-concordance --input-file-a a.bed --input-file-b b.bed --output-file report.txt
EOS
}

sub flags {
    my $self = shift;
    my @flags;
    push(@flags, "--depth") if $self->depth;
    return @flags;
}



sub execute {
    my $self = shift;
    my $output = "-";
    $output = $self->output_file if (defined $self->output_file);
    my $flags = join(" ", $self->flags);
    my $cmd = $self->joinx_path . " snv-concordance $flags " .
        $self->input_file_a . ' ' .
        $self->input_file_b
        . " -o $output";

    my %params = (
        cmd => $cmd,
        input_files => [$self->input_file_a, $self->input_file_b],
    );
    Genome::Sys->shellcmd(%params);

    return 1;
}

sub parse_results_file {
    my $file = shift;

    my $fh = IO::File->new("<$file") || die "Failed to open $file";
    my $results = {};
    my $category;
    my $match_type;
    while (my $l = <$fh>) {
		chomp $l;
        if ($l =~ /^([^\t]+)\t([0-9]+)$/) {
            $category = $1;
            my $total = $2;
            $results->{$category}{total} = $total;
        } elsif ($l =~ /^\t([^\t]+)$/) {
            $match_type = $1;
        } elsif ($l =~ /^\t\t([^\t]+)\t([0-9]+)\t([0-9.]+)$/) {
            my $hit_type = $1;
            my $count = $2;
            my $qual = $3;
            $results->{$category}{hits}{$match_type}{$hit_type}{count} = $count;
            $results->{$category}{hits}{$match_type}{$hit_type}{qual} = $qual;
        } else {
            die "Parse error in results file: '$l'";
        }
    }
    $fh->close;
    return $results;
}

1;
