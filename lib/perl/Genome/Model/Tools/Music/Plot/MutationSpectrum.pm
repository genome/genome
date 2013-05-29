package Genome::Model::Tools::Music::Plot::MutationSpectrum;

use warnings;
use strict;
use Genome;
use FileHandle;
use IO::File;
use Carp;
use Genome::Info::IUB;
use POSIX qw( WIFEXITED );
use Cwd ('getcwd','abs_path');

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Plot::MutationSpectrum {
    is => ['Command::V2'],
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
    plot_title => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => 'Mutation Spectrum Sequence Context',
        doc => 'The title of the plot',
    },
    output_file => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 0,
        #default_value => 'output.pdf',
        doc => 'The name of pdf file to save the plot to',
    },
    file4plot => {
        is_input => 1,
        is_optional => 1,
        is => 'String',
        doc => 'The name of the file to save the sequence context data to be plotted.  If not specified, file will be deleted after use',
    },
    ref_seq => {
        is_input => 1,
        is_optional => 1,
        doc => 'full path to build36/37/n of reference sequence fasta file. Default: build37 fasta',
        example_values => ['/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa'],
    },
    ],
};

sub help_brief { 
    "Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands."; 
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis mutation-spectrum-sequence-context --roi-file=SJMEL001003-0260.alltier.snv --output-file=SJMEL.pdf --plot-title="SJMEL001003-0260 --window-size=10"
EOS
}

sub help_detail {                           
    return <<EOS 
    This tool summarizes the mutation spectrum sequence context.  It produces a stacked barplot for each mutation cateogry showing the proportion of bases around the point of interest (position 0) +/- window_size basepairs.
EOS
}

sub execute {
    my $self = shift;
    $DB::single = 1;

    my $ROI_file = $self->roi_file;
    my $window_size = $self->window_size;
    my $plot_output_file = abs_path($self->output_file);
    my $ref_seq_fasta = $self->ref_seq;
    my $plot_title = $self->plot_title;

    my $plot_input_file;
    if($self->file4plot) {
	$plot_input_file = abs_path($self->file4plot);
    }else {
	my ($fh, $tempfile) = Genome::Sys->create_temp_file;
	$plot_input_file = abs_path($tempfile);
    }

    if($window_size !~ /^\d+$/) {
	print STDERR "--window_size $window_size is not a integer!\n";
	return 0;
    }

    my $mutation_context={ 'A->C' => {},
		           'A->G' => {},
		           'A->T' => {},
		           'C->A' => {},
		           'C->G' => {},
		           'C->T' => {},
                          };

    my $joinxRefstatInput = makeinput4Joinx($ROI_file,$window_size);
    my ($fh, $joinxOUT) = Genome::Sys->create_temp_file;
    my $cmd = "joinx1.6 ref-stats -b $joinxRefstatInput -f $ref_seq_fasta -r  | cut -f 1-4,8 > $joinxOUT ";
    my $return = Genome::Sys->shellcmd(
	cmd => "$cmd",
	);
    unless($return) {
	#$self->error_message("Failed to execute: Returned $return");
	die "$cmd failed to execute!";
    }
    $fh->close;

    my $joinxFH = IO::File->new($joinxOUT) or die "Unable to open the file $joinxOUT due to $!"; 
    #my $roiFH = IO::File->new($ROI_file) or die "Unable to open the file $ROI_file due to $!";
    #open (ROI, $ROI_file) or die "Unable to open the file $ROI_file due to $!";
    while(my $line = $joinxFH->getline) {
	next if($line =~/\#/);
	chomp $line;
	my @list = split(/\t/,$line);
	my ($chr,$pos,$ref,$var) = split(/\_/,$list[3]);
	my $key = join("->",($ref,$var));

	my $rev_compl=0; #boolean to determine if we should reverse complement a seq
	if(!exists($mutation_context->{$key})) {
	    $key = return_mutation_spectrum_category($ref,$var);
	    if(!defined($key)) {
		print STDERR "Warning, cannot classify mutation category for $key, skipping...\n";
		next;
	    }
	    $rev_compl= 1;
	}

	my $relative_pos = -1*$window_size;
	my $seq = $list[4]; #reference sequence fragment
	#my $seq = `samtools faidx $ref_seq_fasta $roi | grep -v '>' `;
	$seq = reverse_complement($seq) if($rev_compl);
	my @seq = split('',$seq);
	foreach my $base (@seq) {
	    if(!exists($mutation_context->{$key}->{$relative_pos})) {  #initialize the counts 
		$mutation_context->{$key}->{$relative_pos} = { 'A'=> 0,'T'=> 0,'C'=> 0,'G'=> 0 };
	    }
	    $mutation_context->{$key}->{$relative_pos}->{$base}++;
	    $relative_pos++;
	}
    }
    $joinxFH->close;
    make_file4plot($mutation_context, $plot_input_file);


    # process plot title
    $plot_title =~ s/ /_/g;
    print $plot_title."\n"; 

    my $plot_cmd = "R --slave --args < " . __FILE__ . ".R $plot_input_file $plot_title $plot_output_file";
    WIFEXITED( system $plot_cmd ) or croak "Couldn't run: $plot_cmd ($?)";

=mod
    my $plot_cmd;
    $plot_cmd = qq{ plot_mutation_spectrum_seq_context(input_file="$plot_input_file",output_file="$plot_output_file",plot_title="$plot_title") };
    my $call = Genome::Model::Tools::R::CallR->create(command=>$plot_cmd, library=> "MutationSpectrum.R");
    $call->execute;
    # Call R to create a plot for the user-specified genes
    my $plot_cmd = "R --slave --args < " . __FILE__ . ".R ";
    WIFEXITED( system $plot_cmd ) or croak "Couldn't run: $plot_cmd ($?)";
=cut

    return 1;
}


sub make_file4plot {
    my $data = shift;
    my $outputfile = shift;
    my $outFH = IO::File->new($outputfile, "w") or die "Can't write to $outputfile\n";
    #open(OUT, ">", $outputfile) or die "Can't write to $outputfile\n";
    foreach my $category (keys %$data) {
	foreach my $rel_pos (sort{$a<=>$b} keys %{ $data->{$category} } ) {
	    foreach my $base ( keys %{ $data->{$category}->{$rel_pos} }   ) {
		my $count = $data->{$category}->{$rel_pos}->{$base};
		$outFH->print("$category\t$rel_pos\t$base\t$count\n");
	    }
	}
    }
    $outFH->close;
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
    my $mutation_spectrum = { 'A->C' => 'A->C',
			      'T->G' => 'A->C',
                              'A->G' => 'A->G',
                              'T->C' => 'A->G',
                              'A->T' => 'A->T',
                              'T->A' => 'A->T',
			      'C->A' => 'C->A',
			      'G->T' => 'C->A',
			      'C->G' => 'C->G',
			      'G->C' => 'C->G',
			      'C->T' => 'C->T',
			      'G->A' => 'C->T',
    };
    if(!exists($mutation_spectrum->{"${ref}->${var}"})) {
	return undef;
    }else {
	return $mutation_spectrum->{"${ref}->${var}"};
    }

}

sub makeinput4Joinx {
    my $file = shift;
    my $window_size = shift;

    open(ROI, $file) or die "Unable to open the file $file due to $!";
    my ($fh, $tempfile) = Genome::Sys->create_temp_file;
    #open(OUT, "> ROI.out") or die "Can't write to ROI.out\n";
    while(<ROI>) {
	chomp;
	next if(/chromo/);
	my @list = split(/\t/,$_);
	my $key = join("_",@list[0,1,3,4]);
	my $win_start = $list[1]-$window_size-1; #zero based for the joinx tool
	my $win_end = $list[1]+$window_size;
	$fh->print("$list[0]\t$win_start\t$win_end\t$key\n");
    }
    close ROI;
    $fh->close;

    return $tempfile;
}

sub rem_white_space {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
