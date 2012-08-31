package Genome::Model::RnaSeq::Command::DetectFusions::TophatFusion;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::TophatFusion {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Base',
    doc => 'run the tophat-fusion transcript fusion detector',
};

sub help_synopsis {
    return <<EOS
 genome model rna-seq detect-fusions tophat-fusion index_dir/ f1.fastq f2.fast2 output_dir/

 genome model rna-seq detect-fusions tophat-fusion index_dir/ f1.fastq f2.fast2 output_dir/ --use-version 1.2.3 --params "-a -b -c"
EOS
}

sub help_detail {
    return <<EOS
Run the tophat-fusion gene fusion detector.

It is used by the RNASeq pipeline to perform fusion detection when the fusion detection strategy is set to something like:
 'tophat-fusion 1.2.3'

EOS
}

sub execute {
    my $self = shift;

    unless($self->_fetch_result('get_or_create')){
        die("Unable to create a software result for tophat-fusion");
    }

    return 1;
}

sub shortcut {
    my $self = shift;
    return $self->_fetch_result('get_with_lock');
}

sub _fetch_result {
    my $self = shift;
    my $method = shift;

    my $result = Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult->$method(
            test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
            version => $self->version,
            alignment_result => $self->alignment_result,
            detector_params => $self->detector_params,
            bowtie_version => $self->build->model->bowtie_version,
    );

    if ($result){
        $self->_link_build_to_result($result);
        return 1;
    }

    return 0;
}

1;
