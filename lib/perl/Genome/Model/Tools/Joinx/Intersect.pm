package Genome::Model::Tools::Joinx::Intersect;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Joinx::Intersect {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file_a => {
            is => 'Text',
            doc => 'Sorted bed file "A"',
            shell_args_position => 1,
        },
        input_file_b => {
            is => 'Text',
            doc => 'Sorted bed file "B"',
            shell_args_position => 2,
        },
    ],
    has_optional_input => [
        output_file => {
            is => 'Text',
            doc => 'The output file (defaults to stdout)',
        },
        miss_a_file => {
            is => 'Text',
            doc => 'Write misses from input "a" to this file',
        },
        miss_b_file => {
            is => 'Text',
            doc => 'Write misses from input "b" to this file',
        },
        first_only => {
            is => 'Boolean',
            default => 0,
            doc => 'Notice only the first thing to hit in b, not the full intersection',
        },
        output_both => {
            is => 'Boolean',
            default => 0,
            doc => 'concatenate intersecting lines in output',
        },
        exact_pos => {
            is => 'Boolean',
            default => 0,
            doc => 'require exact position matches (do not count overlaps)',
        },
        exact_allele => {
            is => 'Boolean',
            default => 0,
            doc => 'require exact allele match. implies --exact-pos',
        },
        iub_match => {
            is => 'Boolean',
            default => 0,
            doc => 'when using --exact-allele, this enables expansion and partial matching of IUB codes',
        },
        dbsnp_match => {
            is => 'Boolean',
            default => 0,
            doc => 'Special mode to match alleles to dbsnp and reverse compliment if they do not match the first time',
        },
    ],
};

sub help_brief {
    "Compute intersection (and optionally difference) of 2 bed files."
}

sub help_synopsis {
    my $self = shift;
    "gmt joinx intersect a.bed b.bed [--output-file=n.bed]"
}

sub flags {
    my $self = shift;

    my @flags;
    my @bool_flags = (
        'first_only',
        'output_both',
        'exact_pos',
        'exact_allele',
        'iub_match',
        'dbsnp_match',
    );
    for my $bf (@bool_flags) {
        if ($self->$bf) {
            my $tmp = "--$bf";
            $tmp =~ tr/_/-/;
            push(@flags, $tmp);
        }
    }
        
    push(@flags, "--miss-a " . $self->miss_a_file) if defined $self->miss_a_file;
    push(@flags, "--miss-b " . $self->miss_b_file) if defined $self->miss_b_file;

    return @flags;
}

sub execute {
    my $self = shift;
    my $output = "-";
    $output = $self->output_file if (defined $self->output_file);
    my $flags = join(" ", $self->flags);
    my $cmd = $self->joinx_path . " intersect $flags " .
        $self->input_file_a . ' ' . 
        $self->input_file_b .
        " -o $output";

    my %params = (
        cmd => $cmd,
        #adukes-sometimes files are empty in pipelines, and shellcommand chokes on existing but empty input files
        #input_files => [$self->input_file_a, $self->input_file_b],
        allow_zero_size_output_files=>1,
    );
    $params{output_files} = [$output] if $output ne "-";
    Genome::Sys->shellcmd(%params);

    return 1;
}

1;
