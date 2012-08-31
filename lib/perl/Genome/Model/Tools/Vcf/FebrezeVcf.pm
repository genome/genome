package Genome::Model::Tools::Vcf::FebrezeVcf;

use strict;
use warnings;
use Genome;
use IO::File;
use Getopt::Long;
use FileHandle;

class Genome::Model::Tools::Vcf::FebrezeVcf {
    is => 'Command',
    has => [
    full_output_file => {
        is => 'Text',
        is_output => 1,
        is_optional => 0,
        doc => "Output vcf -- output is gzip format if the input vcf is gzipped format",
    },
    gtonly_output_file => {
        is => 'Text',
        is_output => 1,
        is_optional => 1,
        doc => "Output vcf with only the GT field -- output is gzip format if the input vcf is gzipped format",
    },
    vcf_file => {
        is => 'Text',
        is_optional => 0,
        doc => "Merged Multisample Vcf containing mutations from all samples",
    },
    remove_nonvariant_sites => {
        is => 'Boolean',
        default => 0,
        is_optional => 1,
        doc => "Remove sites that don't have variants now that it has been sanitized",
    },
    ],
};


sub help_brief {
    "Clean up a vcf, removing any sites that have a FT status that isn't passing or missing"
}


sub help_synopsis {
    <<'HELP';
Clean up a vcf, removing any sites that have a FT status that isn't passing or missing
HELP
}

sub help_detail {
    <<'HELP';
Clean up a vcf, removing any sites that have a FT status that isn't passing or missing
HELP
}

###############
sub execute {                               # replace with real execution logic.
    my $self = shift;

    my $vcf_file = $self->vcf_file;
    my $full_output_file = $self->full_output_file;
    my $gtonly_output_file;
    if ($self->gtonly_output_file) {
        $gtonly_output_file = $self->gtonly_output_file;
    }

    my $inFh_vcf;
    my $output_file1;
    my $output_file2;
    if(Genome::Sys->_file_type($vcf_file) eq 'gzip') {
        $inFh_vcf = Genome::Sys->open_gzip_file_for_reading($vcf_file);
        if ($self->gtonly_output_file) {
            $output_file1 = Genome::Sys->open_gzip_file_for_writing($gtonly_output_file) || die "can't open file\n";
        }
        $output_file2 = Genome::Sys->open_gzip_file_for_writing($full_output_file) || die "can't open file\n";
    }
    else {
        $inFh_vcf = Genome::Sys->open_file_for_reading($vcf_file);
        if ($self->gtonly_output_file) {
            $output_file1 = Genome::Sys->open_file_for_writing($gtonly_output_file) || die "can't open file\n";
        }
        $output_file2 = Genome::Sys->open_file_for_writing($full_output_file) || die "can't open file\n";
    }

    my $remove_nonvariant_sites = $self->remove_nonvariant_sites;
    if ($remove_nonvariant_sites) {
        die "This option (remove-nonvariant-sites) isn't currently supported. Blame science.\n";
    }

    my %filter_options;
    my @sample_ids;
    while(my $line = $inFh_vcf->getline) {
        chomp $line;
        if ($line =~ m/^##/) {
            if ($self->gtonly_output_file) {
                print $output_file1 "$line\n";
            }
            print $output_file2 "$line\n";
            next;
        }
        if ($line =~ m/^#/) {
            my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
            @sample_ids = @samples;
            if ($self->gtonly_output_file) {
                print $output_file1 "$line\n";
            }
            print $output_file2 "$line\n";
            next;
        }

        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
        unless ($line =~ m/PASS/gi && ($filter =~ m/PASS/gi || $filter =~ m/\./)) {
            next;
        }

        my $new_format = "GT";
        if ($self->gtonly_output_file) {
            print $output_file1 "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$new_format";
        }
        print $output_file2 "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format";

        my (@format_fields) = split(/:/, $format);
        my $gt_location; #genotype
        my $ft_location; #filter field
        my $count = 0;
        foreach my $format_info (@format_fields) {
            if ($format_info eq 'GT') {
                $gt_location = $count;
            }
            elsif ($format_info eq 'FT') {
                $ft_location = $count;
            }
            $count++;
        }

        my $matcher = "$chr\t$pos";
        $count = 0;
        foreach my $sample (@samples) {
            my (@sample_fields) = split(/:/, $sample);
            if ($sample eq ".") {
                if ($self->gtonly_output_file) {
                    print $output_file1 "\t./.";
                }
                print $output_file2 "\t$sample";
                $filter_options{"Missing"}++;
            }
            else {
                my $filter_info = $sample_fields[$ft_location];
                $filter_options{$filter_info}++;
                if ($filter_info eq "PASS" || $filter_info eq ".") {
                    print $output_file2 "\t$sample";
                    if ($self->gtonly_output_file) {
                        if ($sample_fields[$gt_location] =~ m/\d/) {
                            print $output_file1 "\t$sample_fields[$gt_location]";
                        }
                        else {
                            print $output_file1 "\t./.";
                        }
                    }
                }
                elsif ($filter_info eq "ForcedGenotype") {
                    if ($self->gtonly_output_file) {
                        print $output_file1 "\t./.";
                    }
                    print $output_file2 "\t.";
                }
                else {
                    if ($self->gtonly_output_file) {
                        print $output_file1 "\t./.";
                    }
                    print $output_file2 "\t.";
                }
            }
            $count++;
        }
        if ($self->gtonly_output_file) {
            print $output_file1 "\n";
        }
        print $output_file2 "\n";
    }
    if ($self->gtonly_output_file) {
        close($output_file1);
    }
    close($output_file2);

    return 1;
}

1;
