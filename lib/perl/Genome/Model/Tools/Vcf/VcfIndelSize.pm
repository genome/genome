package Genome::Model::Tools::Vcf::VcfIndelSize;

use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;
use List::Util qw(shuffle);

my %stats = ();

class Genome::Model::Tools::Vcf::VcfIndelSize {
    is => 'Command',
    has => [
        vcf_files => {
            is => 'Text',
            is_optional => 0,
            doc => "comma-seperated list of VCF or VCF.GZ files",
        },
        keep_only_passing => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => "Only count the indel if a position is labeled PASS (or \".\") in the file.",
        },
        verbose => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => "Extended output including stats of what lines were skipped.",
        },
        subsample_vcf_size => {
            is => 'Text',
            is_optional => 1,
            doc => "If selected, subsamples this many INS and DEL of each size(or all if there aren't that many in any size category",
        },
        subsample_output_file => {
            is => 'Text',
            is_optional => 1,
            doc => "If subsample size selected, outputs the subsample to this file.",
        },

    ],
};


sub help_brief {                            # keep this to just a few words <---
    "Count up indel size"
}


sub help_synopsis {
    <<'HELP';
Count indel sizes in a vcf. Also subsample indels based on size for random checking purposes.
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
    <<'HELP';
Count indel sizes in a vcf. Also subsample indels based on size for random checking purposes.
EXAMPLE: gmt vcf vcf-indel-size --vcf-files file1.vcf,file2.vcf.gz,file3.vcf --keep-only-passing
HELP
}

###############

sub execute {                               # replace with real execution logic.
    my $self = shift;
    my $vcf_files = $self->vcf_files;

    my @vcffiles = split(/,/,$vcf_files);
    if (@vcffiles < 1){
        die ("requires multiple VCF files to be input (comma-sep)")
    }

    my $subsample_output_file;
    if ($self->subsample_output_file) {
        $subsample_output_file = Genome::Sys->open_file_for_writing( $self->subsample_output_file ) || die "can't open file\n";
    }

    my %indel_size_hash;
#subsample_vcf_size
    foreach my $vcffile (@vcffiles) {
        my $inFh;
        if ($vcffile =~ m/gz/) {
            $inFh = Genome::Sys->open_gzip_file_for_reading( $vcffile ) || die "can't open file\n";
        }
        else {
            $inFh = Genome::Sys->open_file_for_reading( $vcffile ) || die "can't open file\n";
        }
        my @samplenames;
        while(my $line = $inFh->getline ) {
            chomp($line);
            if ($line =~ /^\#/){
                if ($line =~ /^\#CHROM/){
                    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samplestuff) = split("\t",$line);
                    @samplenames = @samplestuff;
                }
                next;
            } 

            my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samplestuff) = split("\t",$line);

            if ($self->keep_only_passing && $filter ne "PASS" && $filter ne "."){
#                print "Skipped line because it didn't pass filters: $line\n";
                $stats{$vcffile}{'Filtered'}{$filter}++;
                next;
            }

            if ($ref =~ m/,/ || $alt =~ m/,/) { #skip multiallelic alleles
#                print "Skipped line because it was multiallelic: $line\n";
                $stats{$vcffile}{'Multiallelic'}{$alt}++;
                next;
            }
            $ref = substr($ref, 1);
            $alt = substr($alt, 1);
            my $indel_length;
            my $indel_type;
            if (length($ref) == 0 && length($alt) == 0) { #skip SNVs
                $stats{$vcffile}{'SNV'}{"$ref\t$alt"}++;
                next;
            }
            elsif (length($ref) == 0) {
                $indel_length = length($alt);
                $indel_type = "INS";
            }
            elsif (length($alt) == 0) {
                $indel_length = length($ref);
                $indel_type = "DEL";
            }
            elsif (length($ref) > length($alt)) {
                $indel_length = length($ref) - length($alt);
                $indel_type = "DEL_AMBIGUOUS";
            }
            elsif (length($alt) > length($ref)) {
                $indel_length = length($alt) - length($ref);
                $indel_type = "INS_AMBIGUOUS";
            }
            else {
                print "Skipped line because it was ambiguous as an ins or del: $line\n";
                $stats{$vcffile}{'Ambiguous'}{"$ref\t$alt"}++;
            }
            $stats{$vcffile}{'Indel Length'}{$indel_length}++;

            if ($self->subsample_vcf_size) {
                if (length($ref) == 0) {
                    $ref = 0;
                }
                if (length($alt) == 0) {
                    $alt = 0;
                }
                my $annot_start;
                my $annot_stop;
                if ($indel_type =~ m/INS/) {
                    $annot_start = ($pos);
                    $annot_stop = ($pos+1);
                }
                elsif ($indel_type =~ m/DEL/) {
                    $annot_start = ($pos+1);
                    $annot_stop = ($pos+$indel_length);
                }

                my (@format_fields) = split(/:/, $format);
                my $gt_location; #genotype
                my $count = 0;
                foreach my $format_info (@format_fields) {
                    if ($format_info eq 'GT') {
                            $gt_location = $count;
                    }
                    $count++;
                }

                my @geno;
                unless (defined $gt_location) {
                    print "No GT field in your vcf, this makes indel genotype ambiguous\n";
                    foreach my $sample (@samplestuff){
                        push(@geno,"0/1");
                    }
                }
                else {
                    foreach my $sample (@samplestuff){
                        my (@sample_fields) = split(/:/, $sample);
                        push(@geno,$sample_fields[$gt_location]);
                    }
                }

                #if all samples are 1/1 then this isn't a variant site
                my $valid_variant = 0;
                foreach my $gt (@geno) {
                    if ($gt ne "1/1") {
                        $valid_variant = 1;
                    }
                }
                #dont subsample from places with missing genotypes
                foreach my $gt (@geno) {
                    unless ($gt =~ m/\//) {
                        $valid_variant = 0;
                    }
                }

                if ($valid_variant == 0) {
                    next;
                }

                $stats{$vcffile}{'Indel Length2'}{$indel_length}++;
                my $sample_genotypes = join(",",@geno);
                my $sample_names = join(",",@samplenames);
                my $position = "$chr\t$annot_start\t$annot_stop\t$ref\t$alt\t$indel_length\t$sample_names\t$sample_genotypes";
                $indel_size_hash{$vcffile}{$indel_length}{$indel_type}{$position}++;
            }
        }
        close($inFh);
    }

    if ($self->verbose) {
        foreach my $file (sort keys %stats) {
            print "File:$file\n";
            foreach my $category (sort keys %{$stats{$file}}) {
                if ($category =~ m/Length/) {
                    next;
                }
                print "Stat Category:$category\n";
                foreach my $type (sort keys %{$stats{$file}{$category}}) {
                    my $count = $stats{$file}{$category}{$type};
                    print "$type\t$count\n";
                }
            }
        }
    }

    foreach my $file (sort keys %stats) {
        print "File:$file\n";
        print "Indel size distribution:\n";
        print "Size\tCount\n";
        foreach my $type (sort { $a <=> $b } keys %{$stats{$file}{'Indel Length'}}) {
            my $count = $stats{$file}{'Indel Length'}{$type};
            print "$type\t$count\n";
        }
        if ($self->subsample_vcf_size) {
            print "Genotype Filtered Indel size distribution:\n";
            print "Size\tCount\n";
            foreach my $type (sort { $a <=> $b } keys %{$stats{$file}{'Indel Length2'}}) {
                my $count = $stats{$file}{'Indel Length2'}{$type};
                print "$type\t$count\n";
            }
        }
    }
    if ($self->subsample_vcf_size) {
        my $subsample_number = $self->subsample_vcf_size;
        foreach my $file (sort keys %indel_size_hash) {
            foreach my $indel_length (sort { $a <=> $b } keys %{$indel_size_hash{$file}}) {
                print "Length: $indel_length\n";
                foreach my $indel_type (sort keys %{$indel_size_hash{$file}{$indel_length}}) {
                    print "Type: $indel_type\n";
                    if ($indel_type =~ m/AMBIGUOUS/) { next;}
                    my @positions = keys %{$indel_size_hash{$file}{$indel_length}{$indel_type}}; #random hash order?
                    my $number_of_positions = @positions;
                    if ($number_of_positions < $subsample_number) {
                        foreach my $position (sort @positions) {
                            if ($self->subsample_output_file) {
                                print $subsample_output_file "$position\n";
                            }
                            print "$position\n";
                        }
                    }
                    else {
                        my $down_num = $subsample_number - 1;
                        my @downsample = shuffle @positions;
                        @downsample = @downsample[0..$down_num];
                        foreach my $position (@downsample) {
                            if ($self->subsample_output_file) {
                                print $subsample_output_file "$position\n";
                            }
                            print "$position\n";
                        }
                    }
                }
            }
        }
    }
    return 1;
}


