package Genome::Model::Tools::DetectVariants2::Utilities;

use warnings;
use strict;

use Genome;
use Cwd qw(abs_path);
use File::Spec;
use List::Util qw(first);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    final_result_for_variant_type
    final_result_for_variants_directory
);

sub final_result_for_variant_type {
    my ($results, $variant_type) = @_;

    my @dv2_results = grep($_->class =~ /Genome::Model::Tools::DetectVariants2::Result/, @$results);
    @dv2_results = grep($_->class !~ /::(?:Vcf|LqUnion|Manual)/, @dv2_results);

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

sub final_result_for_variants_directory {
    my $dir = shift;
    my $variant_type = shift;

    unless($variant_type =~ /s$/) {
        $variant_type .= 's';
    }

    my $file = first { -e $_ } (File::Spec->join($dir, "$variant_type.hq"), File::Spec->join($dir, "$variant_type.hq.bed"));
    if($file) {
        my $abs_file = abs_path($file);
        my $alloc = Genome::Disk::Allocation->get_allocation_for_path($abs_file);
        if($alloc) {
            my $owner = $alloc->owner;
            if($owner and $owner->isa('Genome::SoftwareResult')) {
                return $owner;
            }
        }
    }

    return;
}

1;
