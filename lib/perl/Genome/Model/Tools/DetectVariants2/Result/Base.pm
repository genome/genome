package Genome::Model::Tools::DetectVariants2::Result::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Base {
    is => ['Genome::SoftwareResult::Stageable'],
    is_abstract => 1,
    doc => 'This class represents the result of a detect-variants operation. This base class just unites the various result types',
};

sub path {
    my $self = shift;
    my ($str) = @_;

    return join('/', $self->output_dir, $str);
}

sub get_vcf_result {
    my $self = shift;
    my @result = Genome::Model::Tools::DetectVariants2::Result::Vcf->get(
        input_id => $self->id,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
    my $vcf_version = Genome::Model::Tools::Vcf->get_vcf_version;
    @result = grep{ $_->vcf_version eq $vcf_version } @result;
    unless(@result < 2){
        die $self->error_message("Found ".scalar(@result)." vcf results for vcf_version: ".$vcf_version . " and input_id: " . $self->id);
    }
    my $vcf_result = (@result == 1) ? $result[0] : undef;
    return $vcf_result;
}

1;
