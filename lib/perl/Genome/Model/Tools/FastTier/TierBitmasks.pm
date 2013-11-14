package Genome::Model::Tools::FastTier::TierBitmasks;

use strict;
use warnings;
use Bit::Vector;
use Genome;
use UR;
use IO::File;
use Data::Dumper;
use Sys::Hostname;

class Genome::Model::Tools::FastTier::TierBitmasks {
    is => 'Genome::SoftwareResult::Stageable',
    has => [
        reference_sequence_build => {
            type => 'Genome::Model::Build::ReferenceSequence',
            is_input => 1,
            doc => 'Reference sequence to use for tier mask creation, default is NCBI human build36',
        },
        annotation_structures => {
            type => 'Genome::Db::Ensembl::AnnotationStructures',
            is_input => 1,
            doc => 'which annotation structures to use',
        },
        ucsc_directory => { #TODO: UCSC files should be their own software result
            type => 'Text',
            is_input => 1,
            doc => 'The location of phastcons17,28, regulatory regions, etc',
        },
        species => {
            type => 'Text',
            is_input => 1,
        },
    ],
};

sub result_paths {
    return ("tier1.bitmask","tier2.bitmask","tier3.bitmask","tier4.bitmask");
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_staging_directory;

    $self->status_message('Create TierBitmasks');

    my $ref = $self->reference_sequence_build->full_consensus_path('fa');
    my $ref_list_fh = Genome::Sys->open_file_for_reading($self->reference_sequence_build->full_consensus_sam_index_path);
    my @chromosomes;

    while(my $line = $ref_list_fh->getline) {
        chomp $line;
        my ($chr) = split /\t/, $line;
        push @chromosomes, $chr; 
    }
    $ref_list_fh->close;

    my %genome;
    my $genome_size = 0;
    my $masked_genome_size = 0;
    for my $ref_chr (@chromosomes) {
        $self->status_message("Running samtools faidx on $ref_chr");
        unless(open(FAIDX,"samtools faidx $ref $ref_chr |")) {
            die "Couldn't pipe samtools faidx\n";
        }
        my $header = <FAIDX>;
        my $chr = $ref_chr;
        my $chr_length = 0;
        my $cur_nstart = 0;
        my $cur_nstop = -1;
        my @n_blocks;
        while(my $line = <FAIDX>) {
            chomp $line;
            my @sequence = split //, $line;  
            while(my $base = shift @sequence) {
                $chr_length++;
                if($base !~ /[ACTG]/i) {
                    if($cur_nstart) {
                        $cur_nstop = $chr_length;
                    }
                    else {
                        $cur_nstart = $chr_length;
                        $cur_nstop = $chr_length;
                    }
                }
                else {
                    if($cur_nstart) {
                        push @n_blocks, [$cur_nstart, $cur_nstop];
                        $cur_nstart = 0;
                        $cur_nstop = -1;
                    }
                }
            }
        }
        $genome{$chr} = Bit::Vector->new($chr_length);
        while(my $interval = shift @n_blocks) {
            $genome{$chr}->Interval_Fill($interval->[0]-1,$interval->[1]-1);
        }
        $genome_size += $chr_length;
        $masked_genome_size += $chr_length - $genome{$chr}->Norm;
        unless(close(FAIDX)) {
            die "Error reading from samtools faidx pipe\n";
        }
    }


    #This script calculates the number of the bases in the genome covered by each current tier definition (as of 5/26/2009)

    #Tier 1 
    #Tier 1 contains all coding alterations and alterations to rna genes. To calculate its coverage we will scan through the transcript table and add coding exon bases and rna transcript bases to the set. RNA are tracked separately because Tier2 also contains the bases in coding transcripts.

    #Get annotation model from genome model

    printf "Calculated genome size is %u\n", $genome_size;
    printf "Masked genome size is %u\n", $masked_genome_size;

    my $tier1_coding = $self->shadow_genome(\%genome);
    my $tier1_rna = $self->shadow_genome(\%genome);
    my $transcript_iterator;
    my $transcript;
    my @exons;
    #now iterate over all transcripts
    for my $chromosome_name(@chromosomes) {
        #Preload substructures
        my @substructures = Genome::TranscriptStructure->get(
            chrom_name => $chromosome_name,
            data_directory => $self->annotation_structures->output_dir
        );
        $transcript_iterator = Genome::Transcript->create_iterator(
            data_directory => $self->annotation_structures->output_dir,            
            chrom_name => $chromosome_name,            
            reference_build_id => $self->reference_sequence_build->id
        );
        $self->status_message("Parsing $chromosome_name\n");
        unless($transcript_iterator) {
            warn "No iterator because ", Genome::Transcript->error_message, " Skipping to next\n";
            next;
        }
        while( $transcript = $transcript_iterator->next) {
            if($transcript->source ne 'ccds' && $transcript->transcript_status ne 'unknown') {
                #then this transcript is considered for annotation
                @exons = $transcript->cds_exons;
                push @exons, $transcript->introns;
                push @exons, grep { $_->structure_type eq 'rna' } $transcript->ordered_sub_structures;
                my $type;
                for my $exon (@exons) {
                    $type = $exon->structure_type;
                    if($type eq 'rna') {
                        $self->add_substructure_to_set($tier1_rna, $exon, $chromosome_name);
                    } elsif ($type eq 'intron') {
                        $self->add_splice_sites_to_set($tier1_coding,$exon,$chromosome_name);
                    } else {
                        $self->add_substructure_to_set($tier1_coding, $exon, $chromosome_name);
                    }
                    $type = undef;
                }
            }
        }
        Genome::TranscriptStructure->unload;
        Genome::Transcript->unload;
    }

    my $tier1 = $self->union_genomes($tier1_coding, $tier1_rna); 
    $tier1 = $self->difference_genomes($tier1, \%genome);
    undef($tier1_rna); #no longer needed
    $self->write_genome_bitmask($self->temp_staging_directory."/tier1.bitmask", $tier1);
    printf "Tier1 encompasses %u bases. %f%% of the genome\n", $self->bases_covered($tier1), $self->bases_covered($tier1)/$masked_genome_size * 100;
    undef($tier1_coding);

    #Tier2
    #Tier2 contains all coding bases (silent mutations), as well as conserved bases via phastConsElements and non-repeat regulatory
    #Scan through phastConsElements
    my $tier2_conserved_set = $self->shadow_genome(\%genome);
    my @conserved_regions = glob($self->ucsc_directory."/conserved_regions/*");

    for my $conserved_region (@conserved_regions){
        my $fh = Genome::Sys->open_file_for_reading($conserved_region);

        while(my $line = $fh->getline) {
            chomp $line;
            my ($bin, $chr, $start, $end, $lod, $score) = split /\t/, $line;
            $chr =~ s/chr//g;
            $self->add_range_to_set($tier2_conserved_set, $chr, $start, $end) if $score >= 500;
        }
        $fh->close;
    }

    #my $tier2 = union_genomes($tier1_coding, $tier2_conserved_set); 
    #undef($tier2_conserved_set); 
    my $tier2 = $tier2_conserved_set;
    print STDERR "Calculated Tier2 conserved set\n";
    $self->status_message("Calculated Tier2 conserved set\n");
    #printf "Tier2 conserved set encompasses %u bases. %f%% of the genome\n", bases_covered($tier2), bases_covered($tier2)/$masked_genome_size * 100;

    #now do regulatory regions
    #
    #First determine repeats
    my $repeatmasker_regions = $self->shadow_genome(\%genome);

    my @repeatmasker_files = glob($self->ucsc_directory."/rmsk/*");
    for my $file (@repeatmasker_files) {
        my $fh = Genome::Sys->open_file_for_reading($file);
        while(my $line = $fh->getline) {
            chomp $line;
            my @fields = split /\s+/, $line;
            my ($chr, $start, $end) = @fields[5,6,7];
            $chr =~ s/chr//g;
            $self->add_range_to_set($repeatmasker_regions,$chr, $start, $end); 
        }
        $fh->close;
    }
    print STDERR "Calculated repeatmasker regions\n";
    $self->status_message("Calculated repeatmasker regions\n");

    #now take union and calculate the coverage


    #next build up the regulatory annotated regions
    my $regulatory_regions = $self->shadow_genome(\%genome);
    #my @files = qw| targetScanS.txt oreganno.txt tfbsConsSites.txt vistaEnhancers.txt eponine.txt firstEF.txt wgEncodeUcsdNgTaf1ValidH3K4me.txt wgEncodeUcsdNgTaf1ValidH3ac.txt wgEncodeUcsdNgTaf1ValidRnap.txt wgEncodeUcsdNgTaf1ValidTaf.txt polyaDb.txt polyaPredict.txt switchDbTss.txt encodeUViennaRnaz.txt laminB1.txt |;
    my @regulatory_regions = glob($self->ucsc_directory."/regulatory_regions/*");

    for my $file (@regulatory_regions) {
        my $fh = Genome::Sys->open_file_for_reading($file);
        while(my $line = $fh->getline) {
            chomp $line;
            my ($bin, $chr, $start, $end, ) = split /\t/, $line; #ignoring scores on these tables (may be bad)
            $chr =~ s/chr//g;
            $self->add_range_to_set($regulatory_regions,$chr, $start, $end); 
        }
        $fh->close;
    }
    print STDERR "Calculated Tier2 regulatory regions\n";
    $self->status_message("Calculated Tier2 regulatory regions\n");

    #no bins in this file
    my @cpg_islands = glob($self->ucsc_directory."/cpg_islands/*");
    for my $cpg_islands (@cpg_islands){
        my $fh = Genome::Sys->open_file_for_reading($cpg_islands);
        while(my $line = $fh->getline) {
            chomp $line;
            my ($bin, $chr, $start, $end, ) = split /\t/, $line; #ignoring scores on these tables (may be bad)
            $chr =~ s/chr//g;
            $self->add_range_to_set($regulatory_regions,$chr, $start, $end); 
        }
        $fh->close;
    }

    $self->in_place_difference_genomes($regulatory_regions, $repeatmasker_regions); 
    #in_place_difference_genomes($regulatory_regions, $tier2);
    print STDERR "Calculated Tier2 regulatory regions / repeatmasker\n";
    $self->status_message("Calculated Tier2 regulatory regions / repeatmasker\n");
    #printf "Tier2 regulatory set encompasses %u bases. %f%% of the genome\n", bases_covered($regulatory_regions), bases_covered($regulatory_regions)/$masked_genome_size * 100;

    $self->in_place_union_genomes($tier2, $regulatory_regions);
    print STDERR "Calculated Tier2 conserved U regulatory regions / repeatmasker\n";
    $self->status_message("Calculated Tier2 conserved U regulatory regions / repeatmasker\n");
    $self->in_place_difference_genomes($tier2, $tier1); #exclude things hitting Tier1
    print STDERR "Calculated (Tier2 conserved U regulatory regions / repeatmasker) / Tier1\n";
    $self->status_message("Calculated (Tier2 conserved U regulatory regions / repeatmasker) / Tier1\n");
    $self->in_place_difference_genomes($tier2, \%genome); #account for masking
    print STDERR "Calculated (Tier2 conserved U regulatory regions / repeatmasker) / Tier1 / masked genome\n";
    $self->status_message("Calculated (Tier2 conserved U regulatory regions / repeatmasker) / Tier1 / masked genome\n");
    printf "Tier2 encompasses %u bases. %f%% of the genome\n", $self->bases_covered($tier2), $self->bases_covered($tier2)/$masked_genome_size * 100;
    $self->write_genome_bitmask($self->temp_staging_directory."/tier2.bitmask", $tier2);
    #free up some mem?
    undef($regulatory_regions);

    #Tier3
    #Tier3 contains all remaining non-repeat bases
    $self->in_place_union_genomes($tier2, $tier1); #everything in tier1 and tier2
    my $tier3 = $self->complement_genome($tier2); #Everything not in tier1 or tier2
    $self->in_place_difference_genomes($tier3,$repeatmasker_regions); #should be everything not in tier1 or tier2 and not in regulatory regions
    $self->in_place_difference_genomes($tier3, \%genome);
    $self->write_genome_bitmask($self->temp_staging_directory."/tier3.bitmask", $tier3);
    #free up some mem?
    printf "Tier3 encompasses %u bases. %f%% of the genome\n", $self->bases_covered($tier3), $self->bases_covered($tier3)/$masked_genome_size * 100;
    my $tier4 = $self->union_genomes($tier2, $tier3);

    ($tier1,$tier2,$tier3) = (undef,undef,undef);

    $self->in_place_complement_genome($tier4);
    $self->in_place_difference_genomes($tier4, \%genome);
    printf "Tier4 encompasses %u bases. %f%% of the genome\n", $self->bases_covered($tier4), $self->bases_covered($tier4)/$masked_genome_size * 100;
    $self->write_genome_bitmask($self->temp_staging_directory."/tier4.bitmask", $tier4);

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->annotation_structures->add_user(label => 'tier bitmasks', user => $self);

    return $self;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("tier_bitmasks-%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

sub add_substructure_to_set {
    my $self = shift;
    my ($set, $exon, $chromosome) = @_;
    my ($start, $end) = ($exon->structure_start, $exon->structure_stop);
    $self->add_range_to_set($set, $chromosome,$start,$end);
}
sub add_splice_sites_to_set {
    my $self = shift;
    my ($set, $exon, $chromosome) = @_;
    my ($start, $end) = ($exon->structure_start, $exon->structure_stop);
    if(defined($start)){
        $self->add_range_to_set($set, $chromosome, $start, $start+1 );
    }
    if(defined($end)){
        $self->add_range_to_set($set, $chromosome, $end-1, $end );
    }
}



sub add_range_to_set {
    my $self = shift;
    my ($set, $chromosome, $start, $end) = @_;
    ($start, $end) = ($end, $start) if($start > $end);
    my $vector = $set->{$chromosome};
    return unless defined $vector;
    unless(defined $start && defined $end && $start >= 0 && $end <= $vector->Size-1) {
        warn "Invalid range $start $end\n";
        return;
    }
    $vector->Interval_Fill($start,$end); #assuming 0 based coordinates
}



sub write_genome_bitmask {
    my $self = shift;
    my ($filename,$genome_ref) = @_;
    unless($filename) {
        die("No filename of file to write to");
    }
    unless($genome_ref) {
        die("No bitmask to write to file");
    }
    #do some stuff to write this to a file without making it suck
    my $out_fh = IO::File->new($filename,">:raw");
    unless($out_fh) {
        die("Unable to write to " . $filename);
    }
    my $header_string = join("\t", map {$_ => $genome_ref->{$_}->Size()} sort keys %$genome_ref);
    my $write_string = pack 'N/a*', $header_string;
    my $write_result = syswrite($out_fh,$write_string);
    unless(defined $write_result && $write_result == length($write_string)) {
        die("Error writing the header");
    }
    for my $chr (sort keys %$genome_ref) {
        #first write the length in bytes 
        my $chr_write_string = $genome_ref->{$chr}->Block_Read();
        $write_result = syswrite $out_fh, pack("N",length($chr_write_string));
        unless(defined $write_result || $write_result != 4) {
            die("Error writing the length of chromosome $chr");
        }
        $write_result = syswrite $out_fh, $genome_ref->{$chr}->Block_Read();
        unless(defined $write_result || $write_result != length($chr_write_string)) {
            die("Error writing the header");
        }
    }
    $out_fh->close;
    return 1;
}

sub shadow_genome {
    my $self = shift;
    my $genome = shift;
    my %new;
    for my $chr (keys %$genome) {
        $new{$chr} = $genome->{$chr}->Shadow;
    }
    return \%new;
}

sub union_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    my %union;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        $union{$chr} = $genome1->{$chr}->Clone;
        if(exists($genome2->{$chr})) {
            $union{$chr}->Union($genome2->{$chr},$union{$chr});  
        }
    }
    for my $chr (keys %$genome2) {
        next unless defined $genome2->{$chr};
        if(!exists($genome1->{$chr})) {
            $union{$chr} = $genome2->{$chr}->Clone;  
        }
    }
    return \%union;
}

sub in_place_union_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        if(exists($genome2->{$chr})) {
            $genome1->{$chr}->Union($genome1->{$chr},$genome2->{$chr});  
        }
    }
    for my $chr (keys %$genome2) {
        next unless defined $genome2->{$chr};
        if(!exists($genome1->{$chr})) {
            $genome1->{$chr} = $genome2->{$chr}->Clone;  
        }
    }
}

sub difference_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    my %difference;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        $difference{$chr} = $genome1->{$chr}->Clone;
        if(exists($genome2->{$chr})) {
            $difference{$chr}->Difference($difference{$chr},$genome2->{$chr});  
        }
    }
    return \%difference;
}
sub in_place_difference_genomes {
    my $self = shift;
    my ($genome1, $genome2) = @_;
    for my $chr (keys %$genome1) {
        next unless defined $genome1->{$chr};
        if(exists($genome2->{$chr})) {
            $genome1->{$chr}->Difference($genome1->{$chr},$genome2->{$chr});  
        }
    }
}
sub complement_genome {
    my $self = shift;
    my ($genome) = @_;
    my %result;
    for my $chr (keys %$genome) {
        next unless defined $genome->{$chr};
        $result{$chr} = $genome->{$chr}->Clone;
        $result{$chr}->Complement($result{$chr});#in-place calc. Perhaps more mem efficient
    }
    return \%result;
}
sub in_place_complement_genome {
    my $self = shift;
    my ($genome) = @_;
    for my $chr (keys %$genome) {
        next unless defined $genome->{$chr};
        $genome->{$chr}->Complement($genome->{$chr});#in-place calc. Perhaps more mem efficient
    }
}


sub bases_covered {
    my $self = shift;
    my ($genome) = @_;
    my $total = 0;
    for my $chr (keys %$genome) {
        next unless defined $genome->{$chr};
        $total += $genome->{$chr}->Norm();
    }
    return $total;
}
