package Genome::Model::Tools::Sam::Pileup;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Sam::Pileup {
    is  => 'Genome::Model::Tools::Sam',
    has_input => [
        bam_file => {
            is  => 'Text',
            doc => 'Input BAM File',
        },
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'Output Pileup File',
        },
        reference_sequence_path => {
            is => 'Text',
            doc => 'Path to the reference fa or fasta',
        },
        use_bgzip => {
            is => 'Boolean',
            doc => 'Set this to bgzip the output of pileup',
        },
        samtools_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'samtools version',
        },
        samtools_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'samtools params',
        },
    ],
    has_optional_input => [
        region_file => {
            is => 'Text',
            doc => '1-based file listing sites to get calls for',
        },
        view_limit_file => {
            is => 'Text',
            doc => 'true bed file used to limit samtools view input to pileup',
        },
    ],
    has_param => [
        lsf_resource => {
            default => "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>1000 && mem>16000] span[hosts=1] rusage[tmp=1000:mem=16000]' -M 1610612736",
        }
    ],
};

sub execute {
    my $self = shift;
    my $bam = $self->bam_file;
    my $output_file = $self->output_file;

    #check accessibility of inputs
    unless (-s $bam) {
        die $self->error_message("Could not locate BAM or BAM had no size: ".$bam);
    }
    my $refseq_path = $self->reference_sequence_path;
    unless(-s $refseq_path){
        die $self->error_message("Could not locate Reference Sequence at: ".$refseq_path);
    }

    #check to see if the region file is gzipped, if so, pipe zcat output to the -l
    my $rf = $self->region_file;
    my $region_file = "";
    my $view_region_file = "";

    #Create a temporary file and to dump region_limit file into, but in BED format
    my $vlf = Genome::Sys->create_temp_file_path;
    my $vlf_cmd = "bash -c \"zcat ".$rf." | awk \'BEGIN { OFS=\\\"\\t\\\"; }  { print \\\$1,\\\$2-1,\\\$2; }\' | bgzip -c > ".$vlf."\"";
    my $vlf_result = Genome::Sys->shellcmd( cmd => $vlf_cmd);
    unless($vlf_result){
        die $self->error_message("Could not convert region file into view limiting file..");
    }

    #set up the appropriate region/view files for the commmand
    if(defined($self->region_file)){
        if(Genome::Sys->file_is_gzipped($rf)){
            $region_file = "-l <(zcat ".$rf." | awk \'BEGIN { OFS=\\\"\\t\\\"; }  { print \\\$1,\\\$2,\\\$2; }\')";
            $view_region_file = "-L <(zcat ".$vlf.")";
        }elsif (defined($rf)){
            $region_file = "-l <(cat ".$rf." | awk \'BEGIN { OFS=\\\"\\t\\\"; }  { print \\\$1,\\\$2,\\\$2; }\')";
            $view_region_file = "-L $vlf ";
        }
    }

    #if bgzip is set, push the output through bgzip then to disk
    my $out = "> ".$output_file."\"";

    #put the command components together
    my $samtools = $self->path_for_samtools_version($self->samtools_version);
    my $params = $self->samtools_params || "";
    my $cmd = "bash -c \"samtools view -u $view_region_file $bam | $samtools pileup $params -c -f $refseq_path $region_file - $out";

    my $result = Genome::Sys->shellcmd( cmd => $cmd); #, input_files => [$refseq_path, $bam], output_files => [$output_file], skip_if_output_is_present => 1);
    unless($result){
        die $self->error_message("failed to execute cmd: ".$cmd);
    }

    #If we are using bgzip, place the output into bgzip format
    my $temp_sorted_output = Genome::Sys->create_temp_file_path;
    if($self->use_bgzip){
        my $sort = Genome::Model::Tools::Bed::ChromSort->create( input => $output_file, output => $temp_sorted_output);
        unless($sort->execute){
            die $self->error_message("Could not complete sorting of pileup output!");
        }
        my $bgzip_cmd = "cat $temp_sorted_output | bgzip -c > $output_file";
        my $bgzip_result = Genome::Sys->shellcmd( cmd => $bgzip_cmd );
        unless($bgzip_result){
            die $self->error_message("Could not successfully bgzip sorted pileup output!");
        }
    }

    return 1;
}


1;
