package Genome::Model::Tools::RefSeq::Prepare;

use warnings;
use strict;
use File::stat;

my $SAM_DEFAULT = Genome::Model::Tools::Sam->default_samtools_version;

class Genome::Model::Tools::RefSeq::Prepare {
    is => 'Command',
    has => [
            fasta_file => {
                         type => 'String',
                         is_optional => 0,
                         doc => "Fasta file for reference (will be copied over as all_sequences.fa)",
                     },
            reference_name => {
                    type => 'String',
                    is_optional => 0,
                    doc => 'New reference name (to be specified in processing profiles which use this reference)',
            },
            reference_base_path => {
                    type => 'String',
                    default_value=> Genome::Config::reference_sequence_directory() . '/',
                    doc => "Optional alternate base path for the reference.",
            },
            sam_version  => {
                    type => 'String',
                    default_value => $SAM_DEFAULT,
                    doc  => "samtools version to be used, default is $SAM_DEFAULT",
                    is_optional => 1,
            },
    ],
};

sub help_brief {
    "Prepares a fasta file to be used as a new refseq in processing profiles"
}

sub help_detail {
    "Copies a fasta file out to the reference path, and then schedules jobs which will " . 
    "create appropriate BWA, Maq, and Samtools index files."
}

sub execute {
    my $self = shift;
    
    my $path = $self->reference_base_path . "/" . $self->reference_name;

    print "New path for the reference will be $path\n";
    
    if (-d $path) {
        $self->error_message("This path already exists, can't create a new one.");
        return;
    }

    mkdir($path);

    if (!-s $self->fasta_file) {
        $self->error_message("Fasta file " . $self->fasta_file . " doesn't exist or is zero length!");
        return;
    }

    my $new_fasta_path = sprintf("%s/all_sequences.fa", $path); 
    print "Copying to $new_fasta_path\n";

    Genome::Sys->copy_file($self->fasta_file, $new_fasta_path);
    unless (-e $new_fasta_path) {
        $self->error_message("The file did not successfully copy!");
        return;
    }

    my $size = stat($new_fasta_path)->size;

    print "Size is $size\n";
    # bwa bwtsw doesn't work with indexes < 10MB.  give it an 1MB safety margin
    my $bwa_idx_alg = ($size < 11000000 ? "is" : "bwtsw");

    print "Will index with bwa $bwa_idx_alg\n";
   
    my $samtools_version = Genome::Model::Tools::Sam->path_for_samtools_version($self->sam_version);
    print "Using samtools version $samtools_version\n";
    my $bwa_version = Genome::Model::Tools::Bwa->create->bwa_path();
    print "Using bwa version $bwa_version\n";
    my $maq_version = Genome::Model::Tools::Maq->create->maq_path();

    print "Submitting requests to the queue to perform indexing.  You'll be sent an email from LSF with the output when they're done.\n";

    print "Submitting a BWA index request.\n";
    $self->_bsub_invoke(rusage=> ($bwa_idx_alg eq "bwtsw" ? "rusage[mem=4000]' -M 4000000" : ""),
                        job_name => 'bwa-idx',
                        cmd=> $bwa_version . " index -a $bwa_idx_alg $new_fasta_path");
    
    my $new_bfa_path = sprintf("%s/all_sequences.bfa", $path); 
    print "Submitting a Maq fasta2bfa request.\n";
    $self->_bsub_invoke(job_name => 'maq-fasta2bfa',
                        cmd=> $maq_version . " fasta2bfa $new_fasta_path $new_bfa_path");
    
    print "Submitting a Samtools faidx request.\n";
    $self->_bsub_invoke(job_name => 'samtools-faidx',
                        cmd=> $samtools_version . " faidx $new_fasta_path");

    return 1;
}

sub _bsub_invoke {
    my $self = shift;
    my %p = @_;

    my $rusage = delete $p{rusage} || "'";
    my $jobname = delete $p{job_name} || "ref-seq-prep";

    my $bsub64_params = "bsub -J $jobname -u " . Genome::Config->user_email . " -q apipe -R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1] $rusage"; 
    
    my $cmd = delete $p{cmd};
    
    my $bsub_cmd = $bsub64_params . " " . $cmd;
    system($bsub_cmd);

}

