package Genome::Model::Tools::Somatic::FilterPindel;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use File::Basename;

my %insertions;
my %deletions;



class Genome::Model::Tools::Somatic::FilterPindel{
    is => 'Command',
    has => [
        variant_bed_file =>{
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'indels in bed format',
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
        output_dir => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => "place to put results"
        },
    ],
};

sub help_brief {
    "Adjunct temporary pipeline to do tiering dbsnpfiltering and tier1 annotation",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic filter-pindel --output-dir=/gscmnt/sata921/info/medseq/testing_launch_pindel/MEL_calls/filtering_results --variant-bed-file=/gscmnt/sata921/info/medseq/testing_launch_pindel/MEL_calls/indels.hq.bed
EOS
}

sub help_detail {                           
    return <<EOS 
The temporary pipeline to run pindel removed dbsnp filtering, tiering, and annotation so this replaces that until DetectVariants2 is ready.  Supply a bed input and an output directory and it will filter out dbsnp sites, tier the remaining, and annotate them
EOS
}

sub execute {
    my $self = shift;

    my $ifh = Genome::Sys->open_file_for_reading($self->_dbsnp_insertions); #IO::File->new($self->_dbsnp_insertions);
    while (my $line = $ifh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
        next unless ($allele =~ m/-/);
        $allele = substr($allele, 2);
        $insertions{$chr}{$start}{$stop}{'allele'}=$allele;
        $insertions{$chr}{$start}{$stop}{'id'}=$id;
    }
    $ifh->close;
    my $dfh = Genome::Sys->open_file_for_reading($self->_dbsnp_deletions);#IO::File->new($self->_dbsnp_deletions);


    while (my $line = $dfh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $id, $allele, undef) = split /\t/, $line;
        next unless ($allele =~ m/-/);
        $allele = substr($allele, 2);
        $deletions{$chr}{$start}{$stop}{'allele'}=$allele;
        $deletions{$chr}{$start}{$stop}{'id'}=$id;
    }
    $dfh->close;
    Genome::Sys->validate_file_for_reading($self->variant_bed_file);
    my $output = $self->output_dir;
    Genome::Sys->create_directory($output);  #this seems to be a no op if it exists
 
    my $vfh = Genome::Sys->open_file_for_reading($self->variant_bed_file);
    my $output_file_name = $self->output_dir ."/" .  File::Basename::basename($self->variant_bed_file);
    my $novel_file_name = $output_file_name . ".novel";
    my $novel_ofh = Genome::Sys->open_file_for_writing($novel_file_name);
    my $dbsnp_ofh = Genome::Sys->open_file_for_writing($output_file_name . ".dbsnp");
    while(my $bed_line= $vfh->getline) {
        my $dbsnp = $self->dbsnp_lookup($bed_line);
        if($dbsnp eq '-') {
            $novel_ofh->print($bed_line);
        }
        else {
            chomp $bed_line;
            $dbsnp_ofh->print($bed_line . "\t$dbsnp\n");
        }
    }
    $novel_ofh->close;
    $dbsnp_ofh->close;
    my $tier_cmd ="gmt fast-tier fast-tier --variant-bed-file=$novel_file_name";
    Genome::Sys->shellcmd(
        cmd => $tier_cmd,
        input_files => [$novel_file_name],
        skip_if_output_is_present => 0,
    );
    my $tier1 = $novel_file_name . ".tier1";
    Genome::Sys->validate_file_for_reading($tier1);
    my $annotate_cmd = "gmt annotate transcript-variants --annotation=top --variant-bed-file=$tier1 --output=$tier1.annotated";
    Genome::Sys->shellcmd(
        cmd => $annotate_cmd,
        input_files => [$tier1],
        skip_if_output_is_present => 0,
    );


    return 1;
}
1;


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
;
