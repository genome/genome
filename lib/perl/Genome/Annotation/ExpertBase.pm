package Genome::Annotation::ExpertBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::ExpertBase {
    is_abstract => 1,
};

sub dag {
    #   Must return a Genome::WorkflowBuilder::DAG
    # these usually just consist of a build_adaptor
    # followed by a single command, but could be
    # more complex.

    # DAG INPUTS:
    #   build_id
    #   input_result  (Genome::SoftwareResult that has a
    #                  'output_file_path' accessor that refers
    #                  to a .vcf or .vcf.gz file.)
    # DAG OUTPUTS:
    #   software_result (Same requirements as <input_result>)
    die "Abstract";
}

1;
