package Genome::Model::Tools::Snp::DbSnpConcordance;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Bio::DB::Fasta;

class Genome::Model::Tools::Snp::DbSnpConcordance {
    is => 'Command',
    has => [
    output_file =>
    {
        type => 'String',
        doc => "Output file containing the results of the command.",
    },
    snp_file => 
    { 
        type => 'String',
        doc => "maq0.6.8 cns2snp output",
    },
    dbsnp_file =>
    {
        type => 'String',
        doc => "input file of snp locations and calls from the intersection of the affy and illumina platforms",
    },
    exclude_y =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => "Don't consider SNPs present on the Y chromosome",
        default => 0,
    },        
    report_by_quality =>
    {
        type => 'Boolean',
        is_optional => 1,
        doc => "Report number of concordant positions for every possible quality threshold",
        default => 0,
    },        
    ]
};



sub execute {
    my $self=shift;

    #Check on the file names
    unless(-f $self->snp_file) {
        $self->error_message("Snps file is not a file: " . $self->snp_file);
        return;
    }
    unless(-f $self->dbsnp_file) {
        $self->error_message("dbSnp file is not a file: " . $self->dbsnp_file);
        return;
    }

    #Check and open filehandles
    my $snp_fh=IO::File->new($self->snp_file);
    unless($snp_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->snp_file );
        return;
    }
    my $dbsnp_fh=IO::File->new($self->dbsnp_file);
    unless($dbsnp_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->dbsnp_file );
        return;
    }

    my ($dbsnp_hash_ref,) = $self->create_dbsnp_hash($dbsnp_fh);
    close($dbsnp_fh);
    
    unless(defined($dbsnp_hash_ref)) {
        $self->error_message("Fatal error creating dbSNP hash");
        return;
    }

    #Grab metrics
    my ($total_snp_positions,$concordance_ref,) 
        = $self->calculate_metrics($snp_fh,$dbsnp_hash_ref,);
    close($snp_fh);

    $self->print_report($total_snp_positions,$concordance_ref,);
    return 1;
}

    


1;

sub help_brief {
    "Performs a by position comparison of a dbSNP file to a cns2snp output file";
}

sub help_detail {
    "This script performs a comparison of a maq cns2snp output file with a dbSNP file. The comparisons are made by position only as a dbSNP file does not include allele information at this time."
}


sub calculate_metrics {
    my ($self,$snp_fh,$dbsnp_hash_ref,) = @_;
    my $exclude_y = $self->exclude_y;
    
    my %concordance; #with an eye to doing multiple dbs/samples in the future. Do a hash
    my %total_snps;

    #no header in cns2snp
    while(my $line = $snp_fh->getline) {
        chomp $line;
        my ($chr,$pos,$ref,$call,$quality,@metrics) = split /\t/, $line; 

        if($exclude_y) {
            next if($chr eq 'Y'); #female patient these are BS
        }

        next if($ref eq ' ' || $ref eq '' || $ref eq 'N'); #skip 'SNPs' where the reference is N or non-existent
        $total_snps{$quality} += 1;
        
        if(exists($dbsnp_hash_ref->{$chr}{$pos})) {
            $concordance{'dbSNP'}{$quality} += 1;
        }
    }
    return (\%total_snps, \%concordance);

}

sub print_report {
    my ($self, $total_snps_ref, $concordance_ref) = @_; 
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file); 
    my $batch = $self->report_by_quality;
    foreach my $sample (sort keys %$concordance_ref) {
        # make sure we consider all positions in concordance ref OR total_snps_ref
        my $cref = $concordance_ref->{$sample};
        my %quality_levels = map { $_ => 1 } (keys %$cref, keys %$total_snps_ref);
        my $total_concordant = 0;
        my $total_snps = 0;
        foreach my $quality (sort {$b <=> $a} keys %quality_levels) {
            $total_concordant += $cref->{$quality} if exists $cref->{$quality};
            $total_snps += $total_snps_ref->{$quality} if exists $total_snps_ref->{$quality};
            if($batch) {
                print $output_fh "$quality\t$total_concordant\t$total_snps\n";
            }
        }
        print $output_fh "There were $total_snps SNVs with unambiguous reference positions\n";
        printf $output_fh "There were %d positions in %s for a concordance of %0.02f%%\n",$total_concordant,$sample,$total_concordant/$total_snps*100;
    }

    #produce something when there are no results i.e. Solexa testing...
    my $number_of_keys = scalar(keys %$concordance_ref);
    if ($number_of_keys == 0) {
        print $output_fh "There were 0 in the snp file.  No output was generated.\n";
    } 
 
    if (defined($output_fh) ) {
    	$output_fh->close;
    } 
}
#Functions to create hashes of data

#Create hash of dbSNP sites
#right now this is for a single column expressing dbSNP-129 status
sub create_dbsnp_hash {
    my ($self,$dbsnp_fh) = @_;

    my %dbsnp_at;
    #there is no header on this file
    #Format is tab separated
    #Chr\tPos\tPos\tSample1Present\t..SampleNPresent
    while(my $line = $dbsnp_fh->getline) {
        chomp $line;
        my ($chr, $pos, $pos2) = split /\t/, $line;
        $dbsnp_at{$chr}{$pos} = 1;
    }
    return \%dbsnp_at;
}

    
