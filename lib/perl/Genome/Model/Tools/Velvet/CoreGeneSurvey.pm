package Genome::Model::Tools::Velvet::CoreGeneSurvey;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';
use AMOS::AmosLib;

use Bio::Seq;
use Bio::Seq::Quality;
use Bio::SeqIO;

class Genome::Model::Tools::Velvet::CoreGeneSurvey {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
        assembly_directory => {
            is  => 'Text',
            doc => 'Assembly main assembly directory',
        },
        fasta_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'Assembly consensus fasta file',
        },
        percent_id => { 
            is => 'Number',
            doc => 'Acceptable percent identity',
        },
        fraction_of_length => { 
            is => 'Number',
            doc => 'Minimum fraction of the total length of a gene that needs to be present to be counted as covered',
        },
        domain => {
            is => 'Text',
            doc => 'Determines which core gene set/group file to use, either bacteria or archaea',
        },
        output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'Result output file',
        },
    ],
};

sub help_brief {
    'Tool to run core gene survey',
}

sub help_detail {

}

sub execute {
    my $self = shift;

    my @valid_domains = qw/ BACTERIA ARCHAEA UNKNOWN /;
    my $domain = uc $self->domain;
    if ( not grep {$domain eq $_ } @valid_domains ) {
        # don't die if not archaea or bacteria for apipe automated builds 
        $self->status_message('Skipping CoreGeneSurvey running.  Assembly domain must be archaea, bacteria or unknown (unknown domains default to bacteria) to run CoreGeneSurvey');
        return 1;
    }
    $domain = 'BACTERIA' if $domain eq 'UNKNOWN';

    my $input_fasta = ( $self->fasta_file ) ? $self->fasta_file : $self->contigs_bases_file;
    # TODO make contigs.bases if not present and not specified

    my $output_file = ( $self->output_file ) ? $self->output_file : $self->core_gene_survey_file;

    my $tool = Genome::Model::Tools::Bacterial::CoreGeneCoverage->create(
        fasta_file         => $input_fasta,
        survey_type        => 'assembly',
        cell_type          => $domain,
        percent_id         => $self->percent_id,
        fraction_of_length => $self->fraction_of_length,
        output_file        => $output_file,
        run_locally        => 1,
    );

    $self->error_message('Failed to create CoreGeneSurveyTool') and return
        if not $tool;

    $self->error_message('Failed to execute CoreGeneSurveyTool') and return
        if not $tool->execute;

    #zip output file
    
    return 1;
}

1;
