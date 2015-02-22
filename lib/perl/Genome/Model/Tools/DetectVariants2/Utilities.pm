package Genome::Model::Tools::DetectVariants2::Utilities;

use warnings;
use strict;

use Genome;
use Cwd qw(abs_path);
use File::Spec;
use List::Util qw(first);
use List::MoreUtils qw(any uniq);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    final_result_for_variant_type
    final_result_for_variants_directory
);

sub final_result_for_variant_type {
    my ($results, $variant_type) = @_;

    my @dv2_results = grep($_->class =~ /Genome::Model::Tools::DetectVariants2::Result/, uniq @$results);
    @dv2_results = grep($_->class !~ /::(?:Vcf|LqUnion|Manual)/, @dv2_results);

    my @relevant_results = grep(scalar( @{[ glob($_->output_dir . '/' . $variant_type .'*') ]} ), @dv2_results);

    my %results;
    $results{$_->id} = 1 for @relevant_results;
    my @final_results = grep { !any { $results{$_->id} } $_->descendents } @relevant_results;

    if(@final_results > 1) {
        my @combine_results = grep { $_->class =~ /::Combine/ } @final_results;
        @final_results = @combine_results if @combine_results;
    }

    if(!@final_results) {
        return;
    }
    if(@final_results > 1) {
        my $error_string = sprintf('Found multiple results for variant type (%s), Found the following SoftwareResults: %s',
            $variant_type, join(", ", map {$_->id} @final_results));
        Carp::confess($error_string);
    }

    return $final_results[0];
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
