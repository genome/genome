package Genome::Model::Tools::DetectVariants2::Utilities;

use warnings;
use strict;

use Genome;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    final_result_for_variant_type
);

sub final_result_for_variant_type {
    my ($results, $variant_type) = @_;

    my @dv2_results = grep($_->class =~ /Genome::Model::Tools::DetectVariants2::Result/, @$results);
    @dv2_results = grep($_->class !~ /::Vcf/, @dv2_results);
    my @relevant_results = grep(scalar( @{[ glob($_->output_dir . '/' . $variant_type .'*') ]} ), @dv2_results);

    if(!@relevant_results) {
        return;
    }
    if(@relevant_results > 1) {
        my $error_string = sprintf('Found multiple results for variant type (%s), Found the following SoftwareResults: %s',
            $variant_type, join(", ", map {$_->id} @relevant_results));
        Carp::confess($error_string);
    }

    return $relevant_results[0];
}

1;
