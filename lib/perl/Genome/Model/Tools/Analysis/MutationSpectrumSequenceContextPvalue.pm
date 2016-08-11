package Genome::Model::Tools::Analysis::MutationSpectrumSequenceContextPvalue;


use warnings;
use strict;
use Genome;
use FileHandle;
use IO::File;
use Genome::Info::IUB;
use Cwd ('getcwd','abs_path');

class Genome::Model::Tools::Analysis::MutationSpectrumSequenceContextPvalue {
    is => ['Command'],
    has => [
    roi_file => { 
        is  => 'String',
        is_input=>1, 
	is_optional => 0,
        doc => 'name of the input file with SNVs (standard 5 column annotation format) ',
    },
    window_size => {
	is  => 'String',
        is_input=>1, 
	is_optional => 1,
        default_value => '10',
        doc => 'number of bases before and after each position in ROI file to look',
    },
    pvalue_output_file => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 0,
        doc => 'The name of the file to write pvalue info to',
    },
    ref_seq => {
	is_input => 1,
	is_optional => 1,
	doc => 'full path to build36/37/n of reference sequence fasta file. Default: build37 fasta',
	example_values => ['/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa'],
    },
    hypothesized_population_proportion => {
        is_input => 1,
        is_optional => 0,
        doc => 'a fraction between 0 and 1 that the window size of Ts before Cs is compared against',
    },
    location_of_interest => {
        is_input => 1,
        is_optional => 1,
        default => 'before',
        doc => 'Use either "before" or "after", indicating which side of the mutated base you would like a p-value on',
    },
    ],
};

sub help_brief {
    "This tool creates a P-value for the pyrimidine bases before and after C->T transitions using an R function."
}

sub help_detail {                           
    return <<EOS 
    This tool creates a P-value for the pyrimidine bases before and after C->T transitions using an R function.
EOS
    }

sub execute {
    my $self = shift;

    # input params #
    my $ROI_file = $self->roi_file;
    my $window_size = $self->window_size;
    my $pvalue_output_file = abs_path($self->pvalue_output_file);
    my $ref_seq_fasta = $self->ref_seq;
    my $hprop = $self->hypothesized_population_proportion;
    my $location_of_interest = $self->location_of_interest;

    if($window_size !~ /^\d+$/) {
        print STDERR "--window_size $window_size is not a integer!\n";
        return 0;
    }

    my $mutation_context={
#        'A->C' => {},
#        'A->G' => {},
#        'A->T' => {},
#        'C->A' => {},
#        'C->G' => {},
        'C->T' => {},
    };

    # grab reference sequence around snps #
    my $joinxRefstatInput = makeinput4Joinx($ROI_file,$window_size);
    my ($fh2, $joinxOUT) = Genome::Sys->create_temp_file;
    my $cmd = "joinx1.6 ref-stats -b $joinxRefstatInput -f $ref_seq_fasta -r  | cut -f 1-4,8 > $joinxOUT ";
    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
    );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        return;
    }
    $fh2->close;

    # variables to record sample population stats (n and pbar) #
    my $sample_n = 0;
    my $Ts_before_or_after_Cs = 0;
    my $Cs_before_or_after_Cs = 0;

    # read through sites/sequences and record sample stats #
    my $joinxFH = IO::File->new($joinxOUT) or die "Unable to open the file $joinxOUT due to $!"; 
    while (my $line = $joinxFH->getline) {
        next if($line =~/\#/);
        chomp $line;
        my @list = split(/\t/,$line);
        my ($chr,$pos,$ref,$var) = split(/\_/,$list[3]);
        my $key = join("->",($ref,$var));

        my $rev_compl=0; #boolean to determine if we should reverse complement a seq
        if(!exists($mutation_context->{$key})) {
            $key = return_mutation_spectrum_category($ref,$var);
            if(!defined($key)) {
                $key = join("->",($ref,$var));
                print STDERR "Skipping category $key...\n";
                next;
            }
            $rev_compl= 1;
        }

        my $relative_pos = -1*$window_size;
        my $neg_one_pos = $window_size - 1;
        my $plus_one_pos = $window_size + 1;
        my $seq = $list[4]; #reference sequence fragment
        $seq = reverse_complement($seq) if($rev_compl);
        my @seq = split('',$seq);

        # run through sequence before the C position and see if the entire window contains Ts #
        my $all_Ts = 1;
        my $all_Cs = 1;

        if ($location_of_interest eq 'before') {
            # check seq for Ts #
            for my $i (0 .. $neg_one_pos) {
                if ($seq[$i] eq 'T') { next; }
                else { $all_Ts = 0; last; }
            }
            # check seq for Cs #
            for my $i (0 .. $neg_one_pos) {
                if ($seq[$i] eq 'C') { next; }
                else { $all_Cs = 0; last; }
            }
        }

        if ($location_of_interest eq 'after') {
            # check seq for Ts #
            for my $i ($plus_one_pos .. $#seq) {
                if ($seq[$i] eq 'T') { next; }
                else { $all_Ts = 0; last; }
            }
            # check seq for Cs #
            for my $i ($plus_one_pos .. $#seq) {
                if ($seq[$i] eq 'C') { next; }
                else { $all_Cs = 0; last; }
            }
        }

        $sample_n++;
        if ($all_Ts) { $Ts_before_or_after_Cs++; }
        if ($all_Cs) { $Cs_before_or_after_Cs++; }

    }
    $joinxFH->close;

    # find P-value for sample population stats #
    my $pval_cmd = qq{sequence_context_pvalue(p0="$hprop",ts_before_cs="$Ts_before_or_after_Cs",n="$sample_n",outfile="$pvalue_output_file",cs_before_cs="$Cs_before_or_after_Cs",before_or_after="$location_of_interest") };
    print STDERR "Running Command:\n$pval_cmd\n";
    my $call = Genome::Model::Tools::R::CallR->create(command=>$pval_cmd, library=> "MutationSpectrum.R");
    $call->execute;

    return 1;                               
}

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub return_mutation_spectrum_category {

    my $ref = shift;
    my $var = shift;

    my $mutation_spectrum = {
#        'A->C' => 'A->C',
#        'T->G' => 'A->C',
#        'A->G' => 'A->G',
#        'T->C' => 'A->G',
#        'A->T' => 'A->T',
#        'T->A' => 'A->T',
#        'C->A' => 'C->A',
#        'G->T' => 'C->A',
#        'C->G' => 'C->G',
#        'G->C' => 'C->G',
        'C->T' => 'C->T',
        'G->A' => 'C->T',
    };

    if (!exists($mutation_spectrum->{"${ref}->${var}"})) {
        return undef;
    } else {
        return $mutation_spectrum->{"${ref}->${var}"};
    }
}

sub makeinput4Joinx {

    my $file = shift;
    my $window_size = shift;

    open(ROI, $file) or die "Unable to open the file $file due to $!";
    my ($fh3, $tempfile) = Genome::Sys->create_temp_file;
    #open(OUT, "> ROI.out") or die "Can't write to ROI.out\n";
    while(<ROI>) {
        chomp;
        next if(/chromo/);
        my @list = split(/\t/,$_);
        my $key = join("_",@list[0,1,3,4]);
        my $win_start = $list[1]-$window_size-1; #zero based for the joinx tool
        my $win_end = $list[1]+$window_size;
        $fh3->print("$list[0]\t$win_start\t$win_end\t$key\n");
    }
    close ROI;
    $fh3->close;
    return $tempfile;
}

sub rem_white_space {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}
