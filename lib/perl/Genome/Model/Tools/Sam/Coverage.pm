package Genome::Model::Tools::Sam::Coverage;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;


class Genome::Model::Tools::Sam::Coverage {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        aligned_reads_file => {
            is  => 'String',
            doc => 'The input sam/bam file.',
        },
        reference_file => {
            is  => 'String',
            doc => 'The reference file used in the production of the aligned reads file.',
        },
        return_output => {
            is => 'Integer',
            doc => 'Flag to allow the return of the coverage report string.',
            default_value => 1,
        }, 
    ],
    has_optional => [
        output_file => {
            is =>'String',
            doc => 'The output file containing metrics.  If no file is provided a temporary file is created.',
        },
        coverage_command => {
            is =>'String',
            doc => 'The coverage command tool path.',
        },

    ],
};


sub help_brief {
    'Calculate haploid coverage using a BAM file.';
}

sub help_detail {
    return <<EOS
    Calculate haploid coverage using a BAM file.  
EOS
}


sub execute {
    my $self = shift;
    my $reference = $self->reference_file;
    my $aligned_reads = $self->aligned_reads_file;

    my $coverage_cmd;
    if (defined $self->coverage_command) {
       $coverage_cmd = $self->coverage_command;
    } else { 
        #currently this only from r350wu1 while default is r320wu1, so give use_version => r350wu1 should work
        #when we set G::M::T::Sam default to r350wu1, use_version will not be necessary. See Coverage.t
        my $samtools_cmd = $self->samtools_path;
        $coverage_cmd = "$samtools_cmd mapcheck"; 
    }

    my $output_file;
    if (defined $self->output_file) {
        $output_file = $self->output_file;
    } else {
        $output_file = Genome::Sys->create_temp_file_path('mapcheck_coverage_results');
    }
 
    my $cmd = "$coverage_cmd $aligned_reads -f $reference > $output_file";

    $self->status_message("Mapcheck coverage command: ".$cmd);
    my $report_rv = Genome::Sys->shellcmd(cmd=>$cmd,output_files=>[$output_file],input_files=>[$aligned_reads,$reference]);
    
    my @output_text;
    my $return_text; 
    if ($self->return_output) {
        my $out_fh = Genome::Sys->open_file_for_reading($output_file);
        if ($out_fh) {
            @output_text = $out_fh->getlines;
        } 
        $out_fh->close;
        $return_text = join("",@output_text);
        return $return_text;
    } 

    return 1; 

}


1;
