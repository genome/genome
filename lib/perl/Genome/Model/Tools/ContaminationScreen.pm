package Genome::Model::Tools::ContaminationScreen;

use strict;
use warnings;

use Genome;    
use Workflow;
use IO::File;
use Bio::SeqIO;;

class Genome::Model::Tools::ContaminationScreen {
    is => 'Command',
    is_abstract => 1,
    has => [
            input_file =>   {
                                doc => 'file of reads to be checked for contamination',
                                is => 'String',
                                is_input => 1,
                                is_optional => 1,
                            },
            output_file =>  {
                                doc => 'file to write contaminations to',
                                is => 'String',
                                is_output => 1,
                                is_optional => 1,
                            },
            database =>     {
                                doc => 'alignment database', 
                                is => 'String',
                                is_input => 1,
                                is_optional => 1,
                            },
            summary_file => {
                                doc => 'statistical report',
                                is => 'String',
                                is_output => 1,
                                is_optional => 1,
                            },
         ],
};

sub help_brief {
    'Tools for working with solexa metagnomic assemblies.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools metagenomic assembly ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return $self;
}

sub execute
{
    die("implement execute in inheriting class");
}

sub _resolve_directory
{
    my $self = shift;
    my $class_name = ref $self;
    my @parts = split("::", $class_name);
    return '/gsc/var/tmp/fasta/' . $parts[@parts - 1];
}

1;
