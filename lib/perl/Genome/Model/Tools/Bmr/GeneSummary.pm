package Genome::Model::Tools::Bmr::GeneSummary;

use strict;
use warnings;

use Genome;
use IO::File;
use Bit::Vector;
use Benchmark;

class Genome::Model::Tools::Bmr::GeneSummary {
    is => 'Genome::Command::Base',
    has_input => [
    refseq_build_name => {
        is => 'String',
        is_optional => 1,
        default => 'NCBI-human-build36',
        doc => 'The reference sequence build, used to gather and generate bitmask files for base masking and sample coverage.',
    },
    roi_bedfile => {
        type => 'String',
        is_optional => 0,
        doc => 'BED file used to limit background regions of interest when calculating background mutation rate',
    },
    mutation_maf_file => {
        type => 'String',
        is_optional => 0,
        doc => 'List of mutations used to calculate background mutation rate',
    },
    wiggle_file_dirs => {
        type => 'Csv',
        is_optional => 0,
        doc => 'directories containing wiggle files detailing genome-wide coverage of each sample in dataset (comma-delimited)',
    },
    class_summary_file => {
        type => 'String',
        is_optional => 0,
        doc => 'Background Mutation Rates for each class of mutation from this sample set, found using \'gmt bmr class-summary\'.',
    },
    output_file => {
        type => 'String',
        is_optional => 0,
        doc => 'File to contain results table.',
    },
    rejected_mutations => {
        type => 'String',
        is_optional => 1,
        doc => 'File to store mutations that did not fall within the ROIs, or whose gene names did not match any in the ROI list. Default operation is to print to STDOUT.',
    },
    ]
};

sub help_brief {
    "Calculate coverage and number of mutations for every gene in ROIs."
}

sub help_detail {
    return <<HELP;
This script calculates the per-gene coverage and number of mutations found within the regions
specified in the input ROI list. The input mutation list provides the number of non-synonomous
mutations (missense, nonsense, nonstop, splice-site) found in the sample set for each mutation
class. The coverage is found by intersecting the input wiggle files and the regions of interest,
and the bitmasks for each mutation category. Data output is of this format: [Gene  Mutation_Class
Coverage  #Mutations  BMR] for each gene, for each class of mutation.
HELP
}

sub execute {
    my $self = shift;
    my $t0 = Benchmark->new;

    #resolve refseq
    #my $ref_build_name = $self->refseq_build_name;
    #my ( $ref_model_name, $ref_build_version ) = $ref_build_name =~ /^(\S+)-build(\S*)$/;
    #my $ref_model = Genome::Model->get( name=>$ref_model_name );
    #my $ref_build = $ref_model->build_by_version( $ref_build_version );
    #my $ref_dir = $ref_build->data_directory;
    my $ref_dir = "/gscmnt/gc2106/info/medseq/ckandoth/refseq"; #This is much faster
    my $ref_index = $ref_dir . "/all_sequences.fa.fai";

    #WigToBitmask.pm contains some useful functions for handling bitmasks
    my $bitmasker = Genome::Model::Tools::Bmr::WigToBitmask->create(
        reference_index => $ref_index,
    );

    #Load bitmasks
    my $at_bitmask_file = $ref_dir . "/all_sequences.AT_bitmask";
    my $cpg_bitmask_file = $ref_dir . "/all_sequences.CpG_bitmask";
    my $cg_bitmask_file = $ref_dir . "/all_sequences.CG_bitmask";
    my $at_bitmask = $bitmasker->read_genome_bitmask( $at_bitmask_file );
    my $cpg_bitmask = $bitmasker->read_genome_bitmask( $cpg_bitmask_file );
    my $cg_bitmask = $bitmasker->read_genome_bitmask( $cg_bitmask_file );

    #Make sure bitmasks were loaded successfully
    unless ($at_bitmask) {
        $self->error_message("AT bitmask was not loaded.");
        return;
    }
    unless ($cpg_bitmask) {
        $self->error_message("CpG bitmask was not loaded.");
        return;
    }
    unless ($cg_bitmask) {
        $self->error_message("CG bitmask was not loaded.");
        return;
    }

    #load ROIs into a hash %ROIs -> chr -> gene -> start = stop;
    my %ROIs = ();
    my $roi_bedfile = $self->roi_bedfile;
    my $bed_fh = new IO::File $roi_bedfile,"r";
    while (my $line = $bed_fh->getline) {
        chomp $line;
        my ($chr,$start,$stop,$exon_id) = split /\t/,$line;
        #(my $gene = $exon_id) =~ s/^([^\.]+)\..+$/$1/;
        my $gene = $exon_id;
        if ($chr eq "M") { $chr = "MT"; } #for broad roi lists
        if (exists $ROIs{$chr}{$gene}{$start}) {
            next if $stop < $ROIs{$chr}{$gene}{$start};
        }
        $ROIs{$chr}{$gene}{$start} = $stop;
    }
    $bed_fh->close;

    #Create a new ROI bitmask
    my $roi_bitmask = $self->create_empty_genome_bitmask($ref_index);
    for my $chr (keys %ROIs) {
        for my $gene (keys %{$ROIs{$chr}}) {
            for my $start (keys %{$ROIs{$chr}{$gene}}) {
                my $stop = $ROIs{$chr}{$gene}{$start};
                $roi_bitmask->{$chr}->Interval_Fill($start,$stop);
            }
        }
    }

    #Parse wiggle file directories to obtain the path to each wiggle file
    my @wiggle_files = ();
    my @wiggle_dirs = split( /,\s*/, $self->wiggle_file_dirs );
    for my $wiggle_dir (@wiggle_dirs) {
        $wiggle_dir = (( $wiggle_dir =~ m/\/$/ ) ? $wiggle_dir : "$wiggle_dir/" );
        opendir(WIG_DIR, $wiggle_dir) or die "Cannot open directory $wiggle_dir $!\n";
        my @files = readdir(WIG_DIR);
        closedir(WIG_DIR);
        @files = grep { /\.wig$/ } @files;
        @files = map { $wiggle_dir . $_ } @files;
        push(@wiggle_files, @files);
    }

    #Loop through samples to build COVMUTS hash %COVMUTS -> gene -> class -> coverage,mutations
    my %COVMUTS = ();
    my @classes = qw(CG.transit CG.transver AT.transit AT.transver CpG.transit CpG.transver Indels);
    for my $wiggle_file (@wiggle_files) {
        #Load the bitmask with the coverage data for this sample
        $bitmasker->wig_file($wiggle_file);
        $bitmasker->output_file("$wiggle_file.bitmask");
        if ($bitmasker->is_executed) { #This shouldn't have been executed, but just in case
            $bitmasker->is_executed('0');
        }
        #If the bitmask was already created, then don't re-create it, just load it
        if( -e "$wiggle_file.bitmask" ) {
            $bitmasker->read_genome_bitmask("$wiggle_file.bitmask");
        }
        else {
            $bitmasker->execute;
        }

        unless ($bitmasker) {
            $self->error_message("Not seeing a bitmask object for file $wiggle_file.");
            return;
        }
        my $cov_bitmask = $bitmasker->bitmask;
        unless ($cov_bitmask) {
            $self->error_message("Not seeing a hashref to the bitmask.");
            return;
        }

        #find intersection of ROIs and sample's coverage
        for my $chr (keys %$cov_bitmask) {
            $cov_bitmask->{$chr}->And($cov_bitmask->{$chr},$roi_bitmask->{$chr});
        }

        #Calculate coverage of each class in each ROI in this sample
        for my $chr (keys %ROIs) {
            my $chr_length_test_vec = $cov_bitmask->{$chr}->Shadow();
            my $at_cov_vec = $cov_bitmask->{$chr}->Shadow();
            my $cg_cov_vec = $cov_bitmask->{$chr}->Shadow();
            my $cpg_cov_vec = $cov_bitmask->{$chr}->Shadow();
            $at_cov_vec->And($cov_bitmask->{$chr},$at_bitmask->{$chr});
            $cg_cov_vec->And($cov_bitmask->{$chr},$cg_bitmask->{$chr});
            $cpg_cov_vec->And($cov_bitmask->{$chr},$cpg_bitmask->{$chr});

            for my $gene (keys %{$ROIs{$chr}}) {
                unless (grep { /^$gene$/ } keys %COVMUTS) {
                    for my $class (@classes) {
                        $COVMUTS{$gene}{$class}{'coverage'} = 0;
                        $COVMUTS{$gene}{$class}{'mutations'} = 0;
                    }
                }

                for my $start (keys %{$ROIs{$chr}{$gene}}) {
                    my $stop = $ROIs{$chr}{$gene}{$start};

                    #Indels use the coverage in all the ROIs if this gene
                    my $bits = $self->count_interval($cov_bitmask->{$chr},$chr_length_test_vec,$start,$stop);
                    $COVMUTS{$gene}{'Indels'}{'coverage'} += $bits;

                    #AT
                    $bits = $self->count_interval($at_cov_vec,$chr_length_test_vec,$start,$stop);
                    $COVMUTS{$gene}{'AT.transit'}{'coverage'} += $bits;
                    $COVMUTS{$gene}{'AT.transver'}{'coverage'} += $bits;

                    #CG
                    $bits = $self->count_interval($cg_cov_vec,$chr_length_test_vec,$start,$stop);
                    $COVMUTS{$gene}{'CG.transit'}{'coverage'} += $bits;
                    $COVMUTS{$gene}{'CG.transver'}{'coverage'} += $bits;

                    #CpG
                    $bits = $self->count_interval($cpg_cov_vec,$chr_length_test_vec,$start,$stop);
                    $COVMUTS{$gene}{'CpG.transit'}{'coverage'} += $bits;
                    $COVMUTS{$gene}{'CpG.transver'}{'coverage'} += $bits;
                }#end, for my $start
            }#end, for my $gene
        }#end, for my $chr
        undef $cov_bitmask; #clean up any memory
    }#end, for my wiggle file

    #clean up the object memory
    $bitmasker->delete;
    undef $bitmasker;

    #Loop through mutations, assign them to a gene and class in %COVMUTS
    my $mutation_file = $self->mutation_maf_file;
    my $mut_fh = new IO::File $mutation_file,"r";

    #print rejected mutations to a file or to STDOUT
    my $rejects_file = $self->rejected_mutations;
    my $rejects_fh;
    if ($rejects_file) {
        $rejects_fh = IO::File->( $rejects_file, ">" );
    }
    else {
        open $rejects_fh, ">&STDOUT";
    }

    while (my $line = $mut_fh->getline) {
        next if (( $line =~ /^Hugo\_Symbol/ ) || ( $line =~ /^\#/ ));
        chomp $line;
        my @segs = split( /\t/, $line );
        my ($gene,$geneid,$center,$refbuild,$chr,$start,$stop,$strand,$mutation_class,$mutation_type,$ref,$var1,$var2) = @segs;
        my $inCluster = $segs[54];

        #fix broad chromosome name
        $chr =~ s/chr//;
        #Ignore Silent variant and those in Introns, RNA, UTRs, or Flanks
        next if ( $mutation_class =~ m/RNA|Intron|Silent|3'Flank|3'UTR|5'Flank|5'UTR/ );
        #Skip this variant is it's not within the ROIs this job is processing
        next if ($self->count_bits($roi_bitmask->{$chr},$start,$stop) == 0);

        #SNVs
        if ($mutation_type =~ m/SNP|DNP|ONP|TNP/) {
            #if this mutation is non-synonymous
            if ($mutation_class =~ m/Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/) {
                #and if this gene is listed in the ROI list since it is listed in the MAF and passed the bitmask filter
                if (grep { /^$gene$/ } keys %COVMUTS) {
                    #determine the classification for ref A's and T's
                    $ref = substr( $ref, 0, 1 ); #In case of DNPs or TNPs
                    $var1 = substr( $var1, 0, 1 ); #In case of DNPs or TNPs
                    $var2 = substr( $var2, 0, 1 ); #In case of DNPs or TNPs
                    if ($ref eq 'A') {
                        #is it a transition?
                        if ($var1 eq 'G' || $var2 eq 'G') {
                            $COVMUTS{$gene}{'AT.transit'}{'mutations'}++;
                        }
                        #else, it must be a transversion
                        elsif ($var1 =~ /C|T/ || $var2 =~ /C|T/) {
                            $COVMUTS{$gene}{'AT.transver'}{'mutations'}++;
                        }
                        #otherwise, classification is impossible - quit.
                        else {
                            $self->error_message("Cannot classify variant: $gene, chr$chr:$start-$stop, $var1, $var2");
                            return;
                        }
                    }#end, if ref = A
                    elsif ($ref eq 'T') {
                        #is it a transition?
                        if ($var1 eq 'C' || $var2 eq 'C') {
                            $COVMUTS{$gene}{'AT.transit'}{'mutations'}++;
                        }
                        #else, it must be a transversion
                        elsif ($var1 =~ /G|A/ || $var2 =~ /G|A/) {
                            $COVMUTS{$gene}{'AT.transver'}{'mutations'}++;
                        }
                        #otherwise, classification is impossible - quit.
                        else {
                            $self->error_message("Cannot classify variant: $gene, chr$chr:$start-$stop, $var1, $var2");
                            return;
                        }
                    }#end, if ref = T
                    #determine the classification for ref C's and G's
                    elsif ($ref eq 'C') {
                        #is it a transition?
                        if ($var1 eq 'T' || $var2 eq 'T') {
                            #is it inside a CpG island?
                            if ($self->count_bits($cpg_bitmask->{$chr},$start,$stop)) {
                                $COVMUTS{$gene}{'CpG.transit'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{$gene}{'CG.transit'}{'mutations'}++;
                            }
                        }
                        #if not a transition, is it a transversion?
                        elsif ($var1 =~ /G|A/ || $var2 =~ /G|A/) {
                            #is it inside a CpG island?
                            if ($self->count_bits($cpg_bitmask->{$chr},$start,$stop)) {
                                $COVMUTS{$gene}{'CpG.transver'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{$gene}{'CG.transver'}{'mutations'}++;
                            }
                        }
                        #otherwise, classification is impossible - quit.
                        else {
                            $self->error_message("Cannot classify variant: $gene, chr$chr:$start-$stop, $var1, $var2");
                            return;
                        }
                    }#end, if ref = C
                    elsif ($ref eq 'G') {
                        #is it a transition?
                        if ($var1 eq 'A' || $var2 eq 'A') {
                            #is it inside a CpG island?
                            if ($self->count_bits($cpg_bitmask->{$chr},$start,$stop)) {
                                $COVMUTS{$gene}{'CpG.transit'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{$gene}{'CG.transit'}{'mutations'}++;
                            }
                        }
                        #if not a transition, is it a transversion?
                        elsif ($var1 =~ /T|C/ || $var2 =~ /T|C/) {
                            #is it inside a CpG island?
                            if ($self->count_bits($cpg_bitmask->{$chr},$start,$stop)) {
                                $COVMUTS{$gene}{'CpG.transver'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{$gene}{'CG.transver'}{'mutations'}++;
                            }
                        }
                        #otherwise, classification is impossible - quit.
                        else {
                            $self->error_message("Cannot classify variant: $gene, chr$chr:$start-$stop, $var1, $var2");
                            return;
                        }
                    }#end, if ref = G
                    else {
                        warn("Ref DNA is weird: $ref, $gene, chr$chr:$start-$stop, $center\n");
                        next;
                    }
                }#end, if ROI and MAF genes match 
                #if the ROI list and MAF file do not match, quit.
                else {
                    print $rejects_fh ("Variant within ROI, but gene name unknown: $gene, chr$chr:$start-$stop, $center\n");
                    next;
                }
            }#end, if mutation is non-synonymous
            else {
                warn("Variant classification is weird: $gene, $mutation_class, chr$chr:$start-$stop, $center\n");
                next;
            }
        }#end, if mutation is a SNV
        #Indels
        elsif ($mutation_type =~ m/INS|DEL/ ) {
            #verify this gene is listed in the ROI list since it is listed in the MAF and passed the bitmask filter
            if (grep { /^$gene$/ } keys %COVMUTS) {
                $COVMUTS{$gene}{'Indels'}{'mutations'}++;
            }
            else {
                print $rejects_fh ("Variant within ROI, but gene name unknown: $gene, chr$chr:$start-$stop, $center\n");
                next;
            }
        }#end, if mutation is an indel
        else {
            warn("Variant type is weird: $mutation_type, $gene, chr$chr:$start-$stop, $center\n");
            next;
        }
    }#end, loop through MAF
    $mut_fh->close;

    #Parse class-summary file and load %BMR hash (%BMR -> class = bmr)
    my %BMR; # %class_bmr -> class = b.m.r.
    my $summary_file = $self->class_summary_file;
    my $summaryfh = new IO::File $summary_file,"r";
    while (my $line = $summaryfh->getline) {
        next if ($line =~ /Class/);
        chomp $line;
        my ($class,$bmr) = split /\t/,$line;
        $BMR{$class} = $bmr;
    }
    $summaryfh->close;

    #Print all results
    my $output_file = $self->output_file;
    my $out_fh = new IO::File $output_file,"w";
    print $out_fh "Gene\tClass\tBases_Covered\tNon_Syn_Mutations\tBMR\n";
    for my $gene (keys %COVMUTS) {
        for my $class (sort keys %{$COVMUTS{$gene}}) {
            print $out_fh "$gene\t$class\t";
            print $out_fh $COVMUTS{$gene}{$class}{'coverage'} . "\t";
            print $out_fh $COVMUTS{$gene}{$class}{'mutations'} . "\t";
            print $out_fh $BMR{$class} . "\n";
        }
    }
    my $t3 = Benchmark->new;
    print " Total Time: ", timestr(timediff($t3,$t0)), "\n";
    return 1;
}

sub create_empty_genome_bitmask {
    my $self = shift;
    my $ref_index_file = shift;
    my %genome;
    my $ref_fh = new IO::File $ref_index_file,"r";
    while (my $line = $ref_fh->getline) {
        chomp $line;
        my ($chr,$length) = split /\t/,$line;
        $genome{$chr} = Bit::Vector->new($length + 1); #adding 1 for 1-based coordinates
    }
    $ref_fh->close;
    return \%genome;
}

sub count_interval {
    my ($self,$cov_vec,$test_vec,$start,$stop) = @_;
    $test_vec->Empty();
    my $roi_length = $stop - $start + 1;
    $test_vec->Interval_Copy($cov_vec,0,$start,$roi_length);
    my $count = $test_vec->Norm();
    return $count;
}

sub count_bits {
    my ($self,$vector,$start,$stop) = @_;
    my $count = 0;
    for my $pos ($start..$stop) {
        if ($vector->bit_test($pos)) {
            $count++;
        }
    }
    return $count;
}

1;

#more specific hash structure
#%COVMUTS->gene->class->(coverage,#mutations)
#
#classes
#
#CG.C.transit.T   CG.G.transit.A
#CG.C.transver.A CG.G.transver.T
#CG.C.transver.G CG.G.transver.C
#
#CpG.C.transit.T   CpG.G.transit.A
#CpG.C.transver.A CpG.G.transver.T
#CpG.C.transver.G CpG.G.transver.C
#
#AT.A.transit.G   AT.T.transit.C
#AT.A.transver.T AT.T.transver.A
#AT.A.transver.C AT.T.transver.G
#
#and Indels
