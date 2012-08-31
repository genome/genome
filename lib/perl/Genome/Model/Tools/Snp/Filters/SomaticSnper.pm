package Genome::Model::Tools::Snp::Filters::SomaticSnper;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Statistics::R;

class Genome::Model::Tools::Snp::Filters::SomaticSnper
{
    is => 'Command',
    has => [
        subject_metric_file => 
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'File of experimental metrics for the model',
        },
        control_snp_file =>
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'File of experimental metrics for the normal',
        },
        basedir => 
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Basename for the output file',
        },
       ref_seq_id => 
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Chromosome name or something',
        },
        output_file =>
        {
            type => 'String',
            is_output => 1,
            doc => 'the name of the output file',
            calculate => q|
                        return $self->basedir . "/subject_only" . $self->ref_seq_id . ".csv";
                       |
      
        },
    ],
};

sub execute {
    my $self=shift;
    unless( -f $self->subject_metric_file) {
        $self->error_message("Cannot find " . $self->snp_metric_file);
        return;
    }
    unless( -f $self->control_snp_file ) {
        $self->error_message("Cannot find " . $self->snp_metric_file);
        return;
    }  

    unless( -d $self->basedir) {
        mkdir($self->basedir);
        unless( -d $self->basedir) {
            $self->error_message("Unable to create or access basedir: " . $self->basedir);
            return;
        }
    } 

    my $cmd = "gmt snp intersect ";
     $cmd .= "--headers1=1 --f1-only-output=/tmp/tabbed_output";
     $cmd .= " " . $self->subject_metric_file . " " . $self->control_snp_file;
    my $rv = system($cmd . "> /dev/null");
    if(($rv /= 256) != 0 ) {
        $self->error_message("System execution failure what! Command: $cmd\nReturn Value: $rv");
        return;
    }
    $rv = system("cat /tmp/tabbed_output | sed 's/	/,/g' > " . $self->output_file);
   if(($rv /= 256) != 0 ) {
        $self->error_message("System execution failure what! Command: $cmd\nReturn Value: $rv");
        return;
    }
    unlink("/tmp/tabbed_output");

}


    
    
    
    
