package Genome::Model::Tools::Bmr::ClassSummary;

use strict;
use warnings;

use Genome;
use IO::File;
use Bit::Vector;
use Benchmark;

class Genome::Model::Tools::Bmr::ClassSummary {
    is => 'Genome::Command::Base',
    has_input => [
    refseq_build_name => {
        is => 'String',
        is_optional => 1,
        default => 'NCBI-human-build36',
        doc => 'The reference sequence build, used to gather and generate bitmask files for base masking and sample coverage',
    },
    roi_bedfile => {
        type => 'String',
        is_optional => 0,
        doc => 'BED file used to limit background regions of interest when calculating background mutation rate',
    },
    mutation_maf_file => {
        type => 'String',
        is_optional => 0,
        doc => 'List of mutations needed to calculate background mutation rate',
    },
    wiggle_file => {
        type => 'String',
        is_optional => 0,
        doc => 'A wiggle track format file detailing genome-wide coverage of a sample',
    },
    output_file => {
        type => 'String',
        is_optional => 0,
        doc => 'Output file containing BMR for 7 classes for this group of regions, mutations, and a wiggle file',
    },
    genes_to_exclude => {
        type => 'Csv',
        is_optional => 1,
        doc => 'Comma-delimited list of genes to exclude in BMR calculation',
    },
    rejected_mutations => {
        type => 'String',
        is_optional => 1,
        doc => 'File to store mutations that did not fall within the ROIs, or whose gene names did not match any in the ROI list. Default operation is to print to STDOUT.',
    },
    ]
};

sub help_brief {
    "Calculate the background mutation rate for a single wiggle file."
}

sub help_detail {
    return <<HELP;
This script calculates and prints the background mutation rate for transitions, transversions, and
indel classes given a single wiggle file - (results from many files to be combined afterwards;
results denote specific coverage regions), given a specific set of mutations, and given specific
regions of interest. The input mutation list provides the number of non-synonomous mutations
(missense, nonsense, nonstop, splice-site) found in the sample set for each mutation class. Then,
this number of mutations is divided by the total coverage across each base category (A&T, CpG
islands, C&G not in CpG islands) that could lead to that mutation type in each ROI. Data output is
of this format: [Mutation_Class  BMR  Coverage  #_of_Mutations].
HELP
}

sub execute {
    my $self = shift;
    my $t0 = Benchmark->new;

    #parse gene exclusion list
    my $gene_exclusion_list = $self->genes_to_exclude;
    my @genes_to_exclude = ();
    if (defined $gene_exclusion_list) {
        @genes_to_exclude = split ",", $gene_exclusion_list;
    }

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

    #load ROIs into a new ROI hash %ROIs -> chr -> gene -> start = stop;
    my %ROIs = ();
    my $roi_bedfile = $self->roi_bedfile;
    my $bed_fh = new IO::File $roi_bedfile,"r";
    while (my $line = $bed_fh->getline) {
        chomp $line;
        my ($chr,$start,$stop,$exon_id) = split /\t/,$line;
        #(my $gene = $exon_id) =~ s/^([^\.]+)\..+$/$1/;
        my $gene = $exon_id;
        next if (scalar grep { /^$gene$/ } @genes_to_exclude);
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

    #Initialize a hash for recording coverage and mutations: %COVMUTS -> class -> coverage,mutations
    my %COVMUTS = ();
    my @classes = qw(CG.transit CG.transver AT.transit AT.transver CpG.transit CpG.transver Indels);
    for my $class (@classes) {
        $COVMUTS{$class}{'coverage'} = 0;
        $COVMUTS{$class}{'mutations'} = 0;
    }

    #Create or load a bitmask with the coverage data of this sample
    my $t1 = Benchmark->new;
    my $wiggle_file = $self->wiggle_file;
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
    my $t2 = Benchmark->new;

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

    #find intersection of sample's coverage in ROIs and mutation class regions
    for my $chr (keys %$cov_bitmask) {
        my $tempvector = $cov_bitmask->{$chr}->Clone();

        #Indels use the coverage in all the ROIs
        my $bits_on = $tempvector->Norm();
        $COVMUTS{'Indels'}{'coverage'} += $bits_on;

        #AT
        $tempvector->And($cov_bitmask->{$chr},$at_bitmask->{$chr});
        $bits_on = $tempvector->Norm();
        $COVMUTS{'AT.transit'}{'coverage'} += $bits_on;
        $COVMUTS{'AT.transver'}{'coverage'} += $bits_on;

        #CG
        $tempvector->And($cov_bitmask->{$chr},$cg_bitmask->{$chr});
        $bits_on = $tempvector->Norm();
        $COVMUTS{'CG.transit'}{'coverage'} += $bits_on;
        $COVMUTS{'CG.transver'}{'coverage'} += $bits_on;

        #CpG
        $tempvector->And($cov_bitmask->{$chr},$cpg_bitmask->{$chr});
        $bits_on = $tempvector->Norm();
        $COVMUTS{'CpG.transit'}{'coverage'} += $bits_on;
        $COVMUTS{'CpG.transver'}{'coverage'} += $bits_on;
    }

    undef $cov_bitmask; #clean up any memory

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
        $rejects_fh = new IO::File $rejects_file,"w";
    }
    else {
        open $rejects_fh, ">&STDOUT";
    }

    while (my $line = $mut_fh->getline) {
        next if ( $line =~ m/^Hugo\_Symbol/ or $line =~ m/^#/ );
        chomp $line;
        my @segs = split( /\t/, $line );
        my ($gene,$geneid,$center,$refbuild,$chr,$start,$stop,$strand,$mutation_class,$mutation_type,$ref,$var1,$var2) = @segs;
        my $inCluster = $segs[54];

        #fix broad chromosome name
        $chr =~ s/chr//;
        #highly-mutated genes to ignore
        next if (scalar grep { /^$gene$/ } @genes_to_exclude);
        #Ignore Silent variant and those in Introns, RNA, UTRs, or Flanks
        next if ( $mutation_class =~ m/RNA|Intron|Silent|3'Flank|3'UTR|5'Flank|5'UTR/ );
        #make sure mutation is inside the ROIs
        if($self->count_bits($roi_bitmask->{$chr},$start,$stop) == 0) {
            print $rejects_fh ( "Variant does not fall within any ROI: $gene, chr$chr:$start-$stop, $center\n" );
            next;
        }

        #SNVs
        if ($mutation_type =~ m/SNP|DNP|ONP|TNP/) {
            #if this mutation is non-synonymous
            if ($mutation_class =~ m/Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/) {
                #and if this gene is listed in the ROI list since it is listed in the MAF and passed the bitmask filter
                if (grep { /^$gene$/ } keys %{$ROIs{$chr}}) {
                    #determine the classification for ref A's and T's
                    $ref = substr( $ref, 0, 1 ); #In case of DNPs or TNPs
                    $var1 = substr( $var1, 0, 1 ); #In case of DNPs or TNPs
                    $var2 = substr( $var2, 0, 1 ); #In case of DNPs or TNPs
                    if ($ref eq 'A') {
                        #is it a transition?
                        if ($var1 eq 'G' || $var2 eq 'G') {
                            $COVMUTS{'AT.transit'}{'mutations'}++;
                        }
                        #else, it must be a transversion
                        elsif ($var1 =~ /C|T/ || $var2 =~ /C|T/) {
                            $COVMUTS{'AT.transver'}{'mutations'}++;
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
                            $COVMUTS{'AT.transit'}{'mutations'}++;
                        }
                        #else, it must be a transversion
                        elsif ($var1 =~ /G|A/ || $var2 =~ /G|A/) {
                            $COVMUTS{'AT.transver'}{'mutations'}++;
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
                                $COVMUTS{'CpG.transit'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{'CG.transit'}{'mutations'}++;
                            }
                        }
                        #if not a transition, is it a transversion?
                        elsif ($var1 =~ /G|A/ || $var2 =~ /G|A/) {
                            #is it inside a CpG island?
                            if ($self->count_bits($cpg_bitmask->{$chr},$start,$stop)) {
                                $COVMUTS{'CpG.transver'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{'CG.transver'}{'mutations'}++;
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
                                $COVMUTS{'CpG.transit'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{'CG.transit'}{'mutations'}++;
                            }
                        }
                        #if not a transition, is it a transversion?
                        elsif ($var1 =~ /T|C/ || $var2 =~ /T|C/) {
                            #is it inside a CpG island?
                            if ($self->count_bits($cpg_bitmask->{$chr},$start,$stop)) {
                                $COVMUTS{'CpG.transver'}{'mutations'}++;
                            }
                            else {
                                $COVMUTS{'CG.transver'}{'mutations'}++;
                            }
                        }
                        #otherwise, classification is impossible - quit.
                        else {
                            $self->error_message("Cannot classify variant: $gene, chr$chr:$start-$stop, $var1, $var2");
                            return;
                        }
                    }#end, if ref = G
                    else {
                        print $rejects_fh ("Ref DNA is weird: $ref, $gene, chr$chr:$start-$stop, $center\n");
                        next;
                    }
                }#end, if ROI and MAF genes match 
                #if the ROI list and MAF file do not match, record this with a status message.
                else {
                    print $rejects_fh ("Variant within ROI, but gene name unknown: $gene, chr$chr:$start-$stop, $center\n");
                    next;
                }
            }#end, if mutation is missense
            else {
                print $rejects_fh ("Variant classification is weird: $mutation_class, $gene, chr$chr:$start-$stop, $center\n");
                next;
            }
        }#end, if mutation is a SNV
        #Indels
        elsif ($mutation_type =~ m/INS|DEL/ ) {
            #verify this gene is listed in the ROI list since it is listed in the MAF and passed the bitmask filter
            if (grep { /^$gene$/ } keys %{$ROIs{$chr}}) {
                $COVMUTS{'Indels'}{'mutations'}++;
            }
            else {
                print $rejects_fh ("Variant within ROI, but gene name unknown: $gene, chr$chr:$start-$stop, $center\n");
                next;
            }
        }#end, if mutation is an indel
        else {
            print $rejects_fh ("Variant type is weird: $mutation_type, $gene, chr$chr:$start-$stop, $center\n");
            next;
        }
    }#end, loop through MAF
    $mut_fh->close;

    #Loop through COVMUTS hash and tabulate group BMRs into %BMR hash (%BMR -> class = bmr)
    my %BMR;
    for my $class (keys %COVMUTS) {
        my $rate;
        if ($COVMUTS{$class}{'coverage'}) {
            $rate = $COVMUTS{$class}{'mutations'} / $COVMUTS{$class}{'coverage'};
        }
        else {
            $rate = 'No Coverage';
        }
        $BMR{$class} = $rate;
    }

    #Loop through %BMR to print summary file
    my $output_file = $self->output_file;
    my $out_fh = new IO::File $output_file,"w";
    print $out_fh "Class\tBMR\tCoverage(Bases)\tNon_Syn_Mutations\n";
    for my $class (sort keys %BMR) {
        print $out_fh "$class\t$BMR{$class}\t$COVMUTS{$class}{'coverage'}\t$COVMUTS{$class}{'mutations'}\n" or die "\nFailed on: $class $!\n";
    }
    $out_fh->close;
    my $t3 = Benchmark->new;
    print " Total Time: ", timestr(timediff($t3,$t0)), "\n";
    print "Load Wiggle: ", timestr(timediff($t2,$t1)), "\n";
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
