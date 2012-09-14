package Genome::Model::DifferentialExpression::Command::ConvergeTranscripts::Cuffcompare;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::Model::DifferentialExpression::Command::ConvergeTranscripts::Cuffcompare {
    is => ['Command::V2'],
    has => [
        build => { is => 'Genome::Model::Build', id_by => 'build_id', },
    ],
    has_input_output => {
        build_id => {},
    },
};

sub execute {
    my $self = shift;
    my $build = $self->build;
    my $reference_fasta_path = $build->reference_sequence_build->full_consensus_path('fa');
    my $annotation_gtf_path = $build->annotation_build->annotation_file('gtf',$build->reference_sequence_build->id);
    my $output_directory = $build->transcript_convergence_directory;
    
    # Setup as SoftwareResult
    my $cwd = getcwd();
    chdir($output_directory);

    # Start with these params
    #{
    #    include_contained => 1,
    #    generic_gtf_input => 1,
    #    generate_tracking_files => 0,
    #}
    
    my $transcript_convergence_params = eval {
        $self->build->transcript_convergence_params;
    };
    $transcript_convergence_params->{use_version} = $build->model->transcript_convergence_version;
    unless ($transcript_convergence_params->{input_gtf_paths}) {
        $transcript_convergence_params->{input_gtf_paths} = [$annotation_gtf_path];
    }
    $transcript_convergence_params->{reference_fasta_path} = $reference_fasta_path;
    $transcript_convergence_params->{reference_gtf_path} = $annotation_gtf_path;

    unless (Genome::Model::Tools::Cufflinks::Cuffcompare->execute($transcript_convergence_params)) {
        chdir($cwd);
        die('Failed to execute Cuffcompaer with params: '. Data::Dumper::Dumper($transcript_convergence_params));
    }
    chdir($cwd);
    return 1;
}

1;

