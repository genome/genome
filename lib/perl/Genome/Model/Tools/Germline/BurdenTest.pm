package Genome::Model::Tools::Germline::BurdenTest;

use strict;
use warnings;

use Carp qw/confess/;
use Data::Dumper;
use File::Basename qw/dirname/;
use Genome;
use POSIX qw( WIFEXITED );

my $base_R_commands = join("/", dirname(__FILE__), "gstat/burdentest/burdentest.R");

class Genome::Model::Tools::Germline::BurdenTest {
    is => 'Command::V2',
    has_input => [
        option_file => {
            is => 'Text',
            doc => 'File specifying informatin for the burden test to run',
        },
        analysis_data_type => {
            is => 'Text',
            doc => 'Whether or not the trait is binary or quantitative',
            valid_values => [ qw( B Q ) ],
        },
        phenotype_name => {
            is => 'Text',
            doc => 'Name of the trait to test in the phenotype file from the options ile',
        },
        gene_name => {
            is => 'Text',
            doc => 'Name of the gene to analyze',
        },
        permutations => {
            is => 'Integer',
            doc => 'The number of permutations to perform for calculating the p-value',
            default => 10000,
        },
        seed => {
            is => 'Integer',
            doc => 'The seed for the random number generator used to do the permutations',
            default => 123,
        },
        p_value_permutations => {
            is => 'Integer',
            doc => 'The number of permutations to perform on the p-values',
            default => 1,
        },
        covariates => {
            is => 'Text',
            doc => "String of '+' separated names of items from the phenotype file to be used as covariates",
            is_optional => 1,
        },
    ],
};

sub help_brief {
    "Runs one instance of the burden test."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt germline burden-test
EOS
}

sub execute {
    my $self = shift;
    $DB::single = 1;

    #TODO check input files and directory here

    my $rscript = $base_R_commands;
    unless(-e $rscript) {
        die $self->error_message("$rscript for running burden test doesn't exist or is empty");
    }

    my $option_file = $self->option_file;
    unless(-s $option_file) {
        die $self->error_message("$option_file does not exist or is empty");
    }

    my $phenotype = $self->phenotype_name;
    my $datatype = $self->analysis_data_type;
    my $gene = $self->gene_name;
    my $permutations = $self->permutations;
    my $seed = $self->seed;
    my $pval_permutations = $self->p_value_permutations;



    my $cmd = "Rscript $rscript $option_file $datatype $phenotype $gene $permutations:$seed:$pval_permutations";
    if($self->covariates) {
        $cmd .= " " . $self->covariates;
    }

    print "R:\n$cmd\n";
    WIFEXITED(system $cmd) or confess "Couldn't run: $cmd ($?)";

    return 1;
}

1;
