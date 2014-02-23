package Genome::Model::Tools::Analysis::MutationSpectrumSequenceContext;


use warnings;
use strict;
use Genome;
use Workflow;
use FileHandle;
use IO::File;
use Genome::Info::IUB;
use Cwd ('getcwd','abs_path');

class Genome::Model::Tools::Analysis::MutationSpectrumSequenceContext {
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
    proportiontest => {
	is  => 'String',
        is_input=>1,
        is_optional => 0,
        #default_value => 'output.pdf',
        doc => 'The name of the file to save the proportion result to',
    },
    file4plot => {
	is_input => 1,
        is_optional => 1,
	is => 'String',
	doc => 'The name of the file to save the sequence context data to be plotted.  If not specified, file will be deleted after use',
    },
    ref_seq => {
	is_input => 1,
	is_optional => 0,
	doc => 'specify full path to a ref seq fasta',
    },
    random_seed => {
	is_input => 1,
	is_optional => 1,
	doc => 'set the seed for random number generator (useful for generating consistent results for testing purpose.  Use a integer',
    },
    random_trials => {
	is_input => 1,
	is_optional => 1,
	doc => 'number of random sampling - useful for background mutation context calculation',
	default => '10000',
    },
    ],
};

sub help_brief {
    "Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis mutation-spectrum-sequence-context --roi-file=SJMEL001003-0260.alltier.snv --proportiontest SJMEL001003.prop.test --output-file=SJMEL.pdf --plot-title=SJMEL001003-0260 --window-size=10

 build36 refseq fasta  => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fasta",
 build37 refseq fasta  => "/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa",


EOS
}

sub help_detail {
    return <<EOS
    This tool summarizes the mutation spectrum sequence context.  It produces a stacked barplot for each mutation cateogry showing the proportion of bases around the point of interest (position 0) +/- window_size basepairs.  In addition, it also generates 2 barplots showing the background distribution of 4 bases, and compares the proportion of each base at each position with that of the random distribution.

EOS
    }

sub execute {
    my $self = shift;
    $DB::single = 1;


    my $ROI_file = $self->roi_file;
    my $window_size = $self->window_size;
    my $plot_output_file = abs_path($self->output_file);
    my $proportiontestFile = abs_path($self->proportiontest);
    my $random_trials = $self->random_trials;
    my $random_number_seed = $self->random_seed;
    my $ref_seq_fasta = $self->ref_seq;

    die "Ref Seq Fasta: $ref_seq_fasta cannot be found!\n" if(!-e $ref_seq_fasta);


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


    STDERR->print("step1: Calculating Mutation Context based on variants using window_size $window_size\n");
    my ($mutation_context4type,$mutation_context2type) =  generate_mutation_seq_context($ROI_file,$window_size,$ref_seq_fasta);
    STDERR->print("step2: Calculating Background Mutation Context based on $random_trials randomized samplings\n");
    my ($random_context4type,$random_context2type) = generate_random_seq_context($ref_seq_fasta,$window_size,$random_trials,$random_number_seed);
    make_file4plot_4type("${plot_input_file}.4type",$mutation_context4type,$random_context4type);
    make_file4plot_2type("${plot_input_file}.2type",$mutation_context2type,$random_context2type);


    my $call;
    my $R_cmd;
    STDERR->print("step3: Performing proportion test on Mutation vs Background\n");
    prepare_file4_proportion_test_4type($mutation_context4type,$random_context4type,"$proportiontestFile.4type");
    $R_cmd = qq{ compare_prop2populations(input_file="${proportiontestFile}.4type",output_file="${proportiontestFile}.4type") };
    $call = Genome::Model::Tools::R::CallR->create(command=>$R_cmd, library=> "MutationSpectrum.R");
    $call->execute;

    prepare_file4_proportion_test_2type($mutation_context2type,$random_context2type,"$proportiontestFile.2type");
    $R_cmd = qq{ compare_prop2populations(input_file="${proportiontestFile}.2type",output_file="${proportiontestFile}.2type") };
    $call = Genome::Model::Tools::R::CallR->create(command=>$R_cmd, library=> "MutationSpectrum.R");
    $call->execute;

    STDERR->print("step4: Plotting Mutation Context\n");
    $R_cmd = qq{ plot_mutation_spectrum_seq_contextV2(input4type="${plot_input_file}.4type",input2type="${plot_input_file}.2type",output_file="$plot_output_file",plot_title="$plot_title") };
    $call = Genome::Model::Tools::R::CallR->create(command=>$R_cmd, library=> "MutationSpectrum.R");
    $call->execute;


    return 1;
}





sub generate_mutation_seq_context {

    my $ROI_file = shift;          #input variant files (5 columns)
    my $window_size = shift;       #size of the context around each position sampled
    my $ref_fasta = shift;         #absolute path to ref seq fasta


    die "$ref_fasta not found!\n" if(!-e $ref_fasta);

    #data structures to hold the sequence context count
    my $mutation_context1={'A->C' => {},
		           'A->G' => {},
		           'A->T' => {},
		           'C->A' => {},
		           'C->G' => {},
		           'C->T' => {},
                         };
    my $mutation_context2={'A->C' => {},
		           'A->G' => {},
		           'A->T' => {},
		           'C->A' => {},
		           'C->G' => {},
		           'C->T' => {},
                         };


    my $nitrogen_base={ 'A'=>'purine',
                        'G'=>'purine',
                        'C'=>'pyrimidine',
                        'T'=>'pyrimidine',
    };


    my $joinxRefstatInput = makeinput4Joinx($ROI_file,$window_size);
    #my $joinxOUT = "joinx.output";
    my ($fh, $joinxOUT) = Genome::Sys->create_temp_file;
    my $cmd = "joinx1.6 ref-stats -b $joinxRefstatInput -f $ref_fasta -r  | cut -f 1-4,8 > $joinxOUT ";
    system($cmd);

    my $joinx_FH = IO::File->new($joinxOUT) or die "Failed to open the file $joinxOUT\n";
    #open (JOINX, $joinxOUT) or die "Can't read the file $joinxOUT due to $!";
    while(my $line = $joinx_FH->getline) {
	next if($line =~ /\#/);
	chomp $line;
	my @list = split(/\t/,$line);
	my ($chr,$pos,$ref,$var) = split(/\_/,$list[3]);

	my $key = join("->",($ref,$var));
	my $rev_compl=0; #boolean to determine if we should reverse complement a seq
	if(!exists($mutation_context1->{$key})) {
	    $key = return_mutation_spectrum_category($ref,$var);
	    if(!defined($key)) {
		STDERR->print("Warning, cannot classify mutation category for $list[3], skipping...\n");
		next;
	    }
	    $rev_compl= 1;
	}
	my $relative_pos = -1*$window_size;
	my $seq = $list[4]; #reference sequence fragment
	$seq = reverse_complement($seq) if($rev_compl);
	my @seq = split('',$seq);
	foreach my $base (@seq) {
	    $mutation_context1->{$key}->{$relative_pos}->{$base}++;
	    #instead of keeping count for all 4 bases, collapse down to purine,pyrimidines
	    my $pur_pyr = $nitrogen_base->{$base} || 'unknown';
	    $mutation_context2->{$key}->{$relative_pos}->{$pur_pyr}++;
	    $relative_pos++;
	}



    }
    $joinx_FH->close;

    #returns to 2 hashrefs (1 with all 4 bases, 1 with purine/pyrimidines)
    return ($mutation_context1,$mutation_context2);

}


sub generate_random_seq_context {

    my $ref_fasta = shift;                  #absolute path to ref seq fasta
    my $window_size = shift;                #size of the context around each position sampled
    my $number_trials = shift;              #number of random trials to perform
    my $random_number_seed = shift;

    srand($random_number_seed) if(defined($random_number_seed)); #user specified random seed

    die "$ref_fasta not found!\n" if(!-e $ref_fasta);

    my $chr = 2; #chr to randomly sample from
    my $max_pos = 240000000; #max position on chr to sample from  (0-240 million)


    my $mutation_context1={'C'=>{}, 'A'=>{} }; #keeps track of counts for each base of interest
    my $mutation_context2={'C'=>{}, 'A'=>{} }; #keeps track of counts for purines/pyrimdines
    my $base_count={'C'=>0, 'A'=>0 }; #keeps track how many trials done for each base of interest
    my $previously_sampled={}; #holds the random number already used to prevent sampling same base twice

    my $nitrogen_base={ 'A'=>'purine',
                        'G'=>'purine',
                        'C'=>'pyrimidine',
                        'T'=>'pyrimidine',
                       };

    while($base_count->{'C'} < $number_trials || $base_count->{'A'} < $number_trials) {
	my $random_number = int(rand($max_pos));
	next if(exists($previously_sampled->{$random_number})); #avoid randomly sample the same base twice
	$previously_sampled->{$random_number} = 1;
	my $win_start = $random_number - $window_size;
	my $win_end = $random_number + $window_size;
	next if($win_start < 0);
	my $ROI = "${chr}:$win_start-$win_end"; #region to grab seq from ref_seq
	my $seq_fragment = `samtools faidx $ref_fasta $ROI | grep -v '>' `;
	$seq_fragment =~ s/\n//g; #remove all newline in the sequence string
	my $base = substr($seq_fragment,$window_size,1);
	next if($seq_fragment =~ /[^ATCG]/); #skip if the sequence fragment contain at least 1 base that is not A,T,C,G

	if(exists($mutation_context1->{$base})) {
	    if($base_count->{$base} < $number_trials) {
		my $relative_pos = -1*$window_size;
		my @seq = split('',$seq_fragment);
		foreach (@seq) {
		    $mutation_context1->{$base}->{$relative_pos}->{$_}++;
                    #instead of keeping count for all 4 bases, collapse down to purine,pyrimidines
		    my $pur_pyr = $nitrogen_base->{$_} || 'unknown';
		    $mutation_context2->{$base}->{$relative_pos}->{$pur_pyr}++;
		    $relative_pos++;
		}
		$base_count->{$base}++;
	    }
	}

    }


    return ($mutation_context1,$mutation_context2);


}


sub prepare_file4_proportion_test_4type {

  my $contextA = shift;  #based on mutation
  my $contextB = shift;  #based on random
  my $output_file = shift;

  #my $tmp_file = "test.proportion.in";
  my $fh_out = IO::File->new($output_file,"w") or die "Unable to write to $output_file\n";
  my @mutation_class = keys %$contextA;
  foreach my $category(@mutation_class) {
    my ($Base) = $category =~ /([ATCG])->[ATCG]/;
    my $Tot_category_size = $contextA->{$category}->{0}->{$Base} || 0;
    my $Tot_background_size = $contextB->{$Base}->{0}->{$Base} || 0;  #usually number of random trials
    if($Tot_category_size) {
      foreach my $rel_pos (sort{$a<=>$b} keys %{ $contextA->{$category} } ) {
        foreach (qw(A T C G)) {
          my $count = $contextA->{$category}->{$rel_pos}->{$_} || 0;
          my $background_count = $contextB->{$Base}->{$rel_pos}->{$_} || 0;
          $fh_out->print("$category\t$rel_pos\t$_\t$count\t$background_count\t$Tot_category_size\t$Tot_background_size\n");
        }
      }
    }
  }
  $fh_out->close;
}


sub prepare_file4_proportion_test_2type {

  my $contextA = shift; #based on mutation
  my $contextB = shift; #based on random
  my $output_file = shift;

  my $nitrogen_base={ 'A'=>'purine',
    'G'=>'purine',
    'C'=>'pyrimidine',
    'T'=>'pyrimidine',
  };


  #my $tmp_file = "test.proportion.in";
  my $fh_out = IO::File->new($output_file,"w") or die "Unable to write to $output_file\n";
  my @mutation_class = keys %$contextA;
  foreach my $category(@mutation_class) {
    my ($Base) = $category =~ /([ATCG])->[ATCG]/;
    my $Tot_category_size = $contextA->{$category}->{0}->{ $nitrogen_base->{$Base} } || 0;
    my $Tot_background_size = $contextB->{$Base}->{0}->{ $nitrogen_base->{$Base} } || 0;
    if($Tot_category_size) {
      foreach my $rel_pos (sort{$a<=>$b} keys %{ $contextA->{$category} } ) {
        foreach (qw(purine pyrimidine)) {
          my $count = $contextA->{$category}->{$rel_pos}->{$_} || 0;
          my $background_count = $contextB->{$Base}->{$rel_pos}->{$_} || 0;
          #my $Tot_category_size = $contextA->{$category}->{0}->{$_} || 0;
          #my $Tot_background_size = $contextB->{$Base}->{0}->{$_} || 0;
          $fh_out->print("$category\t$rel_pos\t$_\t$count\t$background_count\t$Tot_category_size\t$Tot_background_size\n");
        }
      }
    }
  }
  $fh_out->close;
}



sub make_file4plot_4type {

    my $output_file = shift;     #pop off the 1st argument
    my @mutation_contexts = @_;  #list of hashref that hold mutation context info


    my $output_FH = IO::File->new($output_file,"w") or die "Can't write to $output_file\n";
    foreach my $hashref (@mutation_contexts) {
	next if(ref($hashref) ne 'HASH');
	foreach my $category (keys %$hashref) {
	    foreach my $rel_pos (sort{$a<=>$b} keys %{ $hashref->{$category} } ) {
		foreach my $base(qw(A T C G)) {
		    my $count = $hashref->{$category}->{$rel_pos}->{$base} || 0;
		    $output_FH->print("$category\t$rel_pos\t$base\t$count\n");
		}

	    }
	}
    }
    $output_FH->close;

}

sub make_file4plot_2type {

    my $output_file = shift;     #pop off the 1st argument
    my @mutation_contexts = @_;  #list of hashref that hold mutation context info


    my $output_FH = IO::File->new($output_file,"w") or die "Can't write to $output_file\n";
    foreach my $hashref (@mutation_contexts) {
	next if(ref($hashref) ne 'HASH');
	foreach my $category (keys %$hashref) {
	    foreach my $rel_pos (sort{$a<=>$b} keys %{ $hashref->{$category} } ) {
		foreach my $base(qw(purine pyrimidine)) {
		    my $count = $hashref->{$category}->{$rel_pos}->{$base} || 0;
		    $output_FH->print("$category\t$rel_pos\t$base\t$count\n");
		}

	    }
	}
    }
    $output_FH->close;

}

#---------------------------------------------------------


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

    #open(ROI, $file) or die "Unable to open the file $file due to $!";
    my $input_fh = IO::File->new($file) or die "Failed to open the file $file\n";
    my ($output_fh, $tempfile) = Genome::Sys->create_temp_file;
    #open(OUT, "> ROI.out") or die "Can't write to ROI.out\n";
    while(my $line = $input_fh->getline) {
	chomp $line;
	next if($line =~ /chromo/);
	my @list = split(/\t/,$line);
	my $key = join("_",@list[0,1,3,4]);
	my $win_start = $list[1]-$window_size-1; #zero based for the joinx tool
	my $win_end = $list[1]+$window_size;
	$output_fh->print("$list[0]\t$win_start\t$win_end\t$key\n");
    }
    $input_fh->close;
    $output_fh->close;

    return $tempfile;

}




sub rem_white_space {

    my $string = shift;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
