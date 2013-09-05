package Genome::Model::Tools::Vcf::VcfMakerSniperSomatic;

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
use Genome::Model::Tools::Vcf::Helpers qw/genGT order_chroms/;


class Genome::Model::Tools::Vcf::VcfMakerSniperSomatic {
    is => 'Genome::Model::Tools::Vcf::VcfMakerSomaticBase',
    has => [
        type => {
            is => 'Text',
            doc => "type of variant calls - one of \"snv\" or \"indel\"" ,
            is_optional => 0,
            is_input => 1,
        },
        genome_build => {
            is => 'Text',
            doc => "Reference genome build" ,
            is_optional => 1,
        },
    ],
};


sub help_brief {                            # keep this to just a few words <---
    "Generate Vcf File from sniper output"
}


sub help_synopsis {
    <<'HELP';
    Generate a VCF File from sniper output
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
    my $self = shift;

    my $output_file = $self->output_file;
    my $genome_build = $self->genome_build;
    my $chrom = $self->chrom;
    my $seq_center = $self->seq_center;
    my $skip_header = $self->skip_header;
    my $sniper_file = $self->input_file;
    my $sample_id = $self->sample_id;
    my $type = $self->type;
    my $dbsnp_file = $self->dbsnp_file;
    my $cp_score_to_qual = $self->cp_score_to_qual;


    if(($type ne "snv") && ($type ne "indel")){
        die("\"type\" parameter must be one of \"snv\" or \"indel\"");
    }

###########################################################################
# subs

    #get preceding base using sniper faidx
    sub getPrecedingBase{
        my ($chr,$pos) = @_;
        my $cmd = 'sniper faidx /gscmnt/sata921/info/medseq/cmiller/NCBI-human-build36/$chr.fa $chr:$pos-$pos | tail -n 1';
        return runGetPrecedingBase($chr, $pos, $cmd);
    }

###################################################################

#----------------------------------
    unless ($skip_header){
        print_header($genome_build, $sample_id, $output_file, $seq_center);
    }
    if($chrom) {
        my %sniper_hash = readVariantFile($sniper_file, $chrom, $type, $self->standard_chroms);
        ## add DBsnp labels, if --dbsnp is specified
        if ($dbsnp_file ne ""){
            addDbSnp($dbsnp_file, $chrom, \%sniper_hash)
        }
        print_body($output_file, \%sniper_hash, $cp_score_to_qual);
    }
    else {
        chomp(my @chroms = `cut -f 1 $sniper_file | sort -n | uniq`);

        my @complete_chrom_list = order_chroms(@chroms);
        for my $chrom (@complete_chrom_list) {
            $self->status_message("processing $chrom...");
            my %sniper_hash = readVariantFile($sniper_file, $chrom, $type, $self->standard_chroms);



            ## add DBsnp labels, if --dbsnp is specified
            if ($dbsnp_file ne ""){
                addDbSnp($dbsnp_file, $chrom, \%sniper_hash)
            }

            print_body($output_file, \%sniper_hash, $cp_score_to_qual);
        }
    }
    return 1;
}

1;
