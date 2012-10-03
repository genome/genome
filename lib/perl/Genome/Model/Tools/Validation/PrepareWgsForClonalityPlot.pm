package Genome::Model::Tools::Validation::PrepareWgsForClonalityPlot;

use strict;
use Genome;
use IO::File;
use warnings;


class Genome::Model::Tools::Validation::PrepareWgsForClonalityPlot{
    is => 'Command',
    has => [
	bam_file => {
	    is => 'String',
	    is_optional => 1,
	    doc => 'path to the bam file (to get readcounts)',
	},

	snv_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'File containing snvs in 1-based, 5-col format (chr, st, sp, var, ref). Indels will be skipped',
	},

        output_readcounts_file => {
            is => 'String',
	    is_optional => 1,
	    doc => 'output file containing nicely formatted readcounts',
        },

        output_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'output file suitable for input into clonality plot',
        },

        genome_build => {
            is => 'String',
            is_optional => 1,
	    doc => 'genome build (36, 37lite, or path to fasta file for genome)',
            default => '36',
        },
        
        readcounts_file => {
            is => 'String',
            is_optional => 1,
	    doc => 'path to a formatted readcounts file. If provided, a bam file is not necessary, as these readcounts will be used',
        },

        min_mapping_quality => {
            is => 'Number',
            is_optional => 1,
            default => 20,
            doc => 'mapping quality lower cutoff for bam-readcounts',
        },
        
        ]
};

sub help_brief {
    "Take wgs data with snvs called by sniper or samtools and create a varscan-like file suitable for use by the clonality plotting tool (gmt validation clonality-plot)"
}

sub help_detail {
    "Take wgs data with snvs called by sniper or samtools and create a varscan-like file suitable for use by the clonality plotting tool (gmt validation clonality-plot)"
}



sub execute {
    my $self = shift;
    my $bam_file = $self->bam_file;
    my $snv_file = $self->snv_file;
    my $output_file = $self->output_file;
    my $genome_build = $self->genome_build;
    my $output_readcounts_file = $self->output_readcounts_file;
    my $readcounts_file = $self->readcounts_file;
    my $min_mapping_quality = $self->min_mapping_quality;
    
    
    #one or the other must be given
    if(( defined($readcounts_file) && defined($bam_file) ) || 
       ( !(defined($readcounts_file)) && !(defined($bam_file)) )){
        die "must provide either bam or readcounts file (but not both)";
    }


    my $fasta;
    if ($genome_build eq "36") {
        my $reference_build_fasta_object= Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif ($genome_build eq "37lite") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    } elsif (-e $genome_build) {
        $fasta = $genome_build;
    } else {
        die "genome build must be 36, 37lite, or a path to your genome's fasta file";
    }


    #create temp directory for munging
    my $tempdir = Genome::Sys->create_temp_directory();
    unless($tempdir) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }


    #clean up the snvs file (remove indels)
    my $outFh = open(OUTFILE,">$tempdir/snvs") || die "can't open tmp snv file\n";
    my $inFh = IO::File->new( $snv_file ) || die "can't open snv-file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @fields = split("\t",$line);
        unless (($fields[3] eq "0") || 
            ($fields[3] eq "-") || 
            ($fields[4] eq "0") || 
            ($fields[4] eq "-") || 
            (length($fields[3]) > 1) || 
            (length($fields[4]) > 1)){
            print OUTFILE $line . "\n";
        }        
    }
    close($inFh);
    close(OUTFILE);
    
    if(defined($readcounts_file)){
        `cp $readcounts_file $tempdir/readcounts`;
    } else {
        #now run the readcounting
        my $cmd = "gmt analysis coverage bam-readcount --bam-file $bam_file --output-file $tempdir/readcounts --min-quality-score $min_mapping_quality --variant-file $tempdir/snvs --genome-build $fasta";
        my $return = Genome::Sys->shellcmd(
            cmd => "$cmd",
            );
        unless($return) {
            $self->error_message("Failed to execute: Returned $return");
            die $self->error_message;
        }
    }

   

    if (defined($output_readcounts_file)){
        `cp $tempdir/readcounts $output_readcounts_file`;
    }
    

    
    my %tumHash;

    #store readcounts
    $inFh = IO::File->new( "$tempdir/readcounts" ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @fields = split("\t",$line);
        $tumHash{$fields[0] . "|" . $fields[1]} = $line;
    }

    $outFh = open(OUTFILE,">$output_file") || die "can't open output file\n";
    #match them up with the snvs
    my $inFh2 = IO::File->new( $snv_file ) || die "can't open tmp snv file\n";
    while( my $line = $inFh2->getline )
    {
        chomp($line);
        my @fields = split("\t",$line);
        if(exists($tumHash{$fields[0] . "|" . $fields[1]})){
            my @tum = split("\t",$tumHash{$fields[0] . "|" . $fields[1]});           
            print OUTFILE join("\t",(@fields[0..1],@fields[3..4])) . "\t";
            print OUTFILE "0\t0\t0\t";
            print OUTFILE "NULL\t";
            print OUTFILE join("\t",@tum[4..6]) . "\t";
            print OUTFILE "NULL\tSomatic\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n";
        }
    }
    close($inFh);
    close(OUTFILE);

}
1;
