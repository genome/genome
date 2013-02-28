package Genome::Model::Tools::Tabix::ListChromosomes;

use strict;
use warnings;

use Genome;
use List::Util qw/min/;

class Genome::Model::Tools::Tabix::ListChromosomes {
    is => 'Genome::Model::Tools::Tabix',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Tabix indexed file to query',
            shell_args_position => 1,
        },
        suppress_output => {
            is => 'Boolean',
            doc => "When set, no output is displayed",
            default_value => 0,
        },
    ],
    has_output => [

        chromosomes => {
            is => "ARRAY",
            doc => "Return value, a list of chromosomes",
            is_many => 1,
            is_optional => 1,
        }
    ]
};

sub help_brief {
    return "List chromosomes in a tabix indexed file.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt tabix list-chromosomes in.vcf.gz
EOS
}

sub cmdline {
}

sub execute {
    my $self = shift;
    my $tabix = $self->tabix_path;
    my $input = $self->input_file;
    my $cmd = "$tabix -l $input";
    my @result = qx{$cmd};
    my $ret = $?;
    if ($ret != 0) {
        $self->error_message("Failed to fetch chromosome list (via tabix) for $input");
        return 0;
    }
    if (!$self->suppress_output) {
        print join("", @result);
    }

    chomp @result;
    $self->chromosomes(\@result);

    return 1;
}

1;
