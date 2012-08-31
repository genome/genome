package Genome::Model::Tools::Somatic::DbsnpConcordance;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use IO::File;

my %insertions;
my %deletions;

class Genome::Model::Tools::Somatic::DbsnpConcordance {
    is => 'Command',
    has => [
        indel_bed_file => {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Indels to check for dbsnp, in bed format, chrom start stop ref var - with insertions having start == stop',
        },
        indel_annotation_file => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Indels to check for dbsnp, in annotation format, chrom start stop ref var',
        },
        _dbsnp_insertions => {
            type => 'String',
            is_optional => 1,
            default => '/gscmnt/ams1102/info/info/dbsnp130_indels/insertions_start_stop_adjusted_dbsnp130',
            doc => 'dbsnp insertion file',
        },
        _dbsnp_deletions => {
            type => 'String',
            is_optional => 1,
            default => '/gscmnt/ams1102/info/info/dbsnp130_indels/deletions_adjusted_dbsnp130',
            doc => 'dbsnp deletion file',
        },
        output_file => {
            type => 'String',
            is_optional => 0,
            is_output => 1,
            doc => 'Where to place the output',
        },
    ],
    has_param => [
         lsf_queue => {
             default_value => 'tcga',
         }, 
         lsf_resource => {
             default_value => "-M 6000000 -R 'select[type==LINUX64 && mem>16000] rusage[mem=16000]'",
         },
     ],
};

sub execute {
    my $self = shift;
    unless(defined($self->indel_bed_file) xor defined($self->indel_annotation_file)){
        $self->error_message("You must specify either indel_bed_file or indel_annotation_file, but not both.");
        die $self->error_message;
    }
    my $ifh = Genome::Sys->open_file_for_reading($self->_dbsnp_insertions);
    while (my $line = $ifh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
        next unless ($allele =~ m/-/);
        $allele = substr($allele, 2);
        $insertions{$chr}{$start}{$stop}{'allele'}=$allele;
        $insertions{$chr}{$start}{$stop}{'id'}=$id;
    }
    $ifh->close;
    my $dfh = Genome::Sys->open_file_for_reading($self->_dbsnp_deletions);
    while (my $line = $dfh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
        next unless ($allele =~ m/-/);
        $allele = substr($allele, 2);
        $deletions{$chr}{$start}{$stop}{'allele'}=$allele;
        $deletions{$chr}{$start}{$stop}{'id'}=$id;
    }
    $dfh->close;
    my $anno;
    my $fh;
    if(defined($self->indel_annotation_file)){
        $anno = 1;
        $fh = Genome::Sys->open_file_for_reading($self->indel_annotation_file);
    } else {
        $fh = Genome::Sys->open_file_for_reading($self->indel_bed_file);
    }
    my $output = Genome::Sys->open_file_for_writing($self->output_file);
    while (my $line = $fh->getline){
        chomp $line;
        my ($chr,$start,$stop,$ref,$var) = split /\t/, $line;
        if($ref =~ m/\//){
            ($ref,$var) = split "/", $ref;
        }
        my $refvar = "$ref/$var";
        if($anno){
            if($var ne '0'){
                $start++;
            } else {
                $start--;
            }
        }
        my $bed_line = join("\t",($chr,$start,$stop,$refvar));
        print $output $bed_line."\t".$self->dbsnp_lookup($bed_line)."\n";
    }
    $fh->close;
    $output->close;
}

sub dbsnp_lookup {
    my $self=shift;
    my $bed_line =shift;
    my $dbsnp_id="-";
    chomp $bed_line;
    my ($chr, $start, $stop, $refvar) = split "\t", $bed_line;
    my ($ref,$var) = split "/",$refvar;
    if($ref eq "0") {
        if(exists($insertions{$chr}{$start}{$stop}{'allele'})) {
            if ($var eq $insertions{$chr}{$start}{$stop}{'allele'}) {
                $dbsnp_id=$insertions{$chr}{$start}{$stop}{'id'};
            }
        }
    }
    else {        
        if(exists($deletions{$chr}{$start}{$stop}{'allele'})) {
            if ($ref eq $deletions{$chr}{$start}{$stop}{'allele'}) {
                $dbsnp_id=$deletions{$chr}{$start}{$stop}{'id'};
            }
        } 
    }
    return $dbsnp_id;
}
