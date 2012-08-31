package Genome::Model::Tools::Hgmi::CoreGenes;

use strict;
use warnings;

use Genome;
use Carp 'confess';
use Cwd 'abs_path';
use File::Temp;

class Genome::Model::Tools::Hgmi::CoreGenes (
    is => 'Command',
    has => [
        cell_type => { 
            is => 'String',
            doc => 'Type of genome to check',
            valid_values => ['ARCHAEA', 'BACTERIA'],
        },
        sequence_set_id => { 
            is => 'String',
            doc => 'sequence set id of organism',
        },
    ],
    has_optional => [
        dev => { 
            is => 'Boolean',
            doc => 'development flag',
            default => 0,
        },
    ],
);

my %CORE_GENE_PARAMS = (
    'BACTERIA' => {
        percent_id => '30',
        fraction_of_length => '0.3',
    },
    'ARCHAEA' => {
        percent_id => '50',
        fraction_of_length => '0.7',
    },
);
    
sub help_brief { return 'Runs a core gene check on the predicted sequences' }

sub help_synopsis { return help_brief() }

sub help_detail {
    return 'Takes a bunch of predicted sequences (retrieved via sequence set id) ' . 
        'and makes sure the coverage of core genes meets some minimum percentage';
}

sub execute {
    my $self = shift;

    my ($protein_seq_fasta_fh, $protein_seq_fasta) = File::Temp->tempfile(
        "coregenes_XXXXXX",
        DIR => './', # FIXME Should not place files into CWD!
    );
    $protein_seq_fasta_fh->close;

    # Grabs all proteins from the sequence set that passed merging and writes their sequence to the supplied fasta
    my $protein_export_cmd = Genome::Model::Tools::Bacterial::ExportProteins->create(
        sequence_set_id => $self->sequence_set_id,
        dev => $self->dev,
        output_file => $protein_seq_fasta
    );
    confess 'Could not create protein export command!' unless $protein_export_cmd;
    confess 'Could not execute protein export for sequence set ' . $self->sequence_set_id unless $protein_export_cmd->execute;

    # Convert Coregene_results to an absolute path
    # TODO Still need to not write this to cwd...
    my $core_gene_output_file = abs_path('Coregene_results');
    
    my $core_gene_cmd = Genome::Model::Tools::Bacterial::CoreGeneCoverage->create(
        fasta_file => $protein_seq_fasta,
        percent_id => $CORE_GENE_PARAMS{$self->cell_type}{percent_id},
        fraction_of_length => $CORE_GENE_PARAMS{$self->cell_type}{fraction_of_length},
        cell_type => $self->cell_type,
        output_file => $core_gene_output_file, 
    );
    confess 'Could not create core gene command object!' unless $core_gene_cmd;
    confess 'Could not execute core gene check!' unless $core_gene_cmd->execute;

    confess "Core genes check failed, results can be found in $core_gene_output_file" unless $core_gene_cmd->_passed;
    $self->status_message("Core gene check passed, results can be found in $core_gene_output_file");
    return 1;
}

1;

