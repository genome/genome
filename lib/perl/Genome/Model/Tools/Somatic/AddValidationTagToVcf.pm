package Genome::Model::Tools::Somatic::AddValidationTagToVcf;

use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;
use List::MoreUtils qw(firstidx);
use List::MoreUtils qw(uniq);
use Data::Dumper;


class Genome::Model::Tools::Somatic::AddValidationTagToVcf {
    is => 'Command',
    has => [
    output_file => {
        is => 'Text',
        is_output => 1,
        doc => "List of mutations in Vcf format",
    },
    validation_statuses => {
        is => 'Text',
        doc => "File associating chromosome position ref var with a validation status. Tab-separated. No IUB codes.",
    },
    input_vcf => {
        is => 'Text',
        doc => "VCF of variants",
    },
    ],
};


sub help_brief {                            # keep this to just a few words <---
    "Add validation info to a Vcf File"
}


sub help_synopsis {
    <<'HELP';
Add validation info to a VCF File
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
    <<'HELP';
Parses the relevant files and creates a VCF containing all the SNVs. This includes those that fail filters (noted in the FILTER field).
HELP
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    #header lines look like this.
    my $self = shift;

    my $status_fh = IO::File->new($self->validation_statuses,"r");
    unless($status_fh) {
        $self->error_message("Unable to open " . $self->validation_statuses);
        return;
    }

    my %statuses;

    while(my $line = $status_fh->getline) {
        chomp $line;
        next if $line =~/^#/;
        my ($chr, $pos, $ref, $var, $status) = split /\t/, $line;
        $statuses{$chr}{$pos}{$ref}{$var} = $status;        
    }
    $status_fh->close;

    #open vcf file and do some stuff

    my $validation_format_header_line = q{##FORMAT=<ID=VVS,Number=1,Type=Integer,Description="Validation Status relative to non-adjacent reference normal 0=FP, 1=SNP, 2=Somatic, 3=Loss_Of_Heterozygosity, 4=Validation_Failed, 5=Not_Targeted">};
    
    my $vcf_fh = IO::File->new($self->input_vcf);
    unless($vcf_fh) {
        $self->error_message("Unable to open vcf file: ". $self->input_vcf);
        return;
    }

    my $vcf_ofh = IO::File->new($self->output_file,"w");
    unless($vcf_ofh) {
        $self->error_message("Unable to open vcf file for writing: ". $self->output_file);
        return;
    }

    while(my $vcf_line = $vcf_fh->getline) {
        chomp $vcf_line;
        if($vcf_line =~ /^##/) {
            #pass through
            print $vcf_ofh $vcf_line,"\n";
        }
        elsif($vcf_line =~ /^#CHROM/) {
            print $vcf_ofh $validation_format_header_line,"\n";
            print $vcf_ofh $vcf_line,"\n";
        }
        else {
            my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $normal, $tumor) = split /\t/, $vcf_line;
            $format .= ":VVS";
            $normal .= ":.";
            my @alleles = ($ref, split /,/, $alt);
            my ($tumor_genotype) = split /:/, $tumor;
            my @tumor_indices = split /\//, $tumor_genotype;
            if(exists($statuses{$chr}{$pos})) {
                #recode to vcf status and append the result. Not sure how this will work for multiple alt alleles if one is not real
                my $status = "";
                for my $ref (keys %{$statuses{$chr}{$pos}}) {
                    for my $var (keys %{$statuses{$chr}{$pos}{$ref}}) {
                        #unless($status) {
                            $status = $statuses{$chr}{$pos}{$ref}{$var};
                        #}
                        #else {
                        #    $self->error_message("Multiple variant alleles with validation results at $chr $pos. Not sure what to do.");
                        #    return;
                        #}
                    }
                }
                if($status =~ /^\s*Somatic\s*$/i) {
                    $tumor .= ":2";
                }
                elsif($status =~ /^\s*Reference\s*$/i) {
                    $tumor .= ":0";
                }
                elsif($status =~ /^\s*Germline\s*$/i) {
                    $tumor .= ":1";
                }
                elsif($status =~ /^\s*LOH\s*$/i) {
                    $tumor .= ":3";
                }
                else {
                    #assume validation failed
                    $self->debug_message("$status unrecognized as a validation status for $chr $pos. Reporting as a failed validation.");
                    $tumor .= ":4";
                }

            }
            else {
                $tumor .= ":5";
            }
            print $vcf_ofh join("\t",$chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $normal, $tumor),"\n";
        }
    }
    return 1;
}
