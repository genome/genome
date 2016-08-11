package Genome::Model::Tools::Analysis::MutationSpectrumByStrand;

#####################################################################################################################################
# MutationSpectrumByStrand - Given a somatic variation model ID or a file of SNVs, will plot mutation spectrum barplot
#
#			
#####################################################################################################################################

use warnings;
use strict;
use Math::Round;
use Genome;
use FileHandle;
use IO::File;
use Genome::Info::IUB;
use Cwd ('getcwd','abs_path');

class Genome::Model::Tools::Analysis::MutationSpectrumByStrand {
    is => ['Command'],
    has => [
    labels => { 
        is  => 'String',
        is_input=>1, 
	is_optional => 1,
	default_value => "Label",
        doc => 'explicitly specifify sample label (useful when trying to plot multiple samples',
    },
    annotation_file => {
	is  => 'String',
        is_input=>1, 
	is_optional => 0,
        doc => 'a five column, (chr,pos1,pos2,ref,var) input file ',
    },
    plot_title => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => 'Mutation Spectrums By Strand',
        doc => 'The title of the plot',
    },
    output_file => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 0,
        #default_value => 'output.pdf',
        doc => 'The name of pdf file to save the plot to',
    },
    mut_spec_file => {
	is_input => 1,
        is_optional => 1,
	is => 'String',
	doc => 'The name of the file to save the mutation spectrum data to be plotted.  If not specified, file will be deleted after use',
    },
    plot_width => {
	is => 'String',
	is_optional => 1,
	default => 6,
        doc => 'Set the width of the plot in inches',

    },
    number_row => {
	is => 'String',
	is_optional => 1,
	default => 2,
        doc => 'Only modify this if you are plotting multiple samples.',
    },
    plot_height => {
	is => 'String',
	is_optional => 1,
	default => 6,
        doc => 'Set the height of the plot in inches',

    },
    plot_graph => {
	is => 'Boolean',
	is_optional => 1,
	doc => 'set this flag if you want the results to be plotted, if set to false, be sure to specify --mut-spec-file',
	default => 1,
    },
    plot_this => {
	is_input => 1,
        is_optional => 1,
	is => 'String',
	doc => 'Plot the user-edited file and exit immediately. ',
    },

    ],
};

sub help_brief {
    "Given an annotation file, gives an output of transition/transversion for those on transcribed/untranscribed strand.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"

gmt analysis mutation-spectrum-by-strand  --annotation-file=SNV.file --output-file SJMELmutspec_strand.pdf --plot-title=SJMEL001002 --mut-spec-file=file4plot.txt
gmt analysis mutation-spectrum-by-strand  --annotation-file=x --output-file SJMELmutspec_strand.pdf --plot-title=SJMEL001002 --plot-this="User.edited.SNV"

EOS
}

sub help_detail {                           
    return <<EOS 
    This tool summarizes the mutation spectrum for a single sample.  It takes a list of SNVs, and only use ones that are mapped to to the genic and intrageneic region (UTR, exon,intron).  It classify each SNV as transcribed/untranscribed by comparing to the strand of the gene that it(SNV) maps to.  Current version only support build37 annotation, so do not forget to run lift-over if you have build36 coordinates. It outputs a barplot of counts for each mutation spectrum category with respect to the strand.

EOS
    }

sub execute {
    my $self = shift;
    $DB::single = 1;

    my $numberRow = $self->number_row;
    #my $group_id = $self->group_id;
    my $plot_input_file;
    my $plot_width = $self->plot_width;
    my $plot_height = $self->plot_height;
    my $plot_title = $self->plot_title;
    my $manual_label = $self->labels;
    my $plot_graph = $self->plot_graph;
    my $SNV_file = $self->annotation_file;

    if($self->mut_spec_file) {
	$plot_input_file = abs_path($self->mut_spec_file);
	unlink($plot_input_file) if(-e $plot_input_file); #remove existing file first
    }else {
	my ($fh, $tempfile) = Genome::Sys->create_temp_file;
	$plot_input_file = abs_path($tempfile);
    }
    my $plot_output_file = abs_path($self->output_file);

    



    if($self->plot_this) {
	$plot_input_file = abs_path($self->plot_this);
	my $plot_cmd;
	$plot_cmd = qq{ plot_mutation_spectrum_bystrand(inputFile="$plot_input_file",outputFile="$plot_output_file",plot_title="$plot_title",num_row=$numberRow,file_width=$plot_width,file_height=$plot_height) };

	my $call = Genome::Model::Tools::R::CallR->create(command=>$plot_cmd, library=> "MutationSpectrum.R");
	$call->execute;
	return 1;

    }

    # TODO: cleanup and normalize thies
    my $genefile = Genome::Sys->dbpath('ensembl','67_37l_v2') . '/exonIntron.blocks.bed';

    my $SNVs = intersectSNV2genes($SNV_file,$genefile);

    my $x=1;

    my $strand_rule = { 'C->T' => {'+'=>'untranscribe',
                                   '-'=>'transcribe',
		              },
                    'G->A' => {'+'=>'transcribe',
                               '-'=>'untranscribe',
		              },

		    'C->G' => {'+'=>'untranscribe',
                               '-'=>'transcribe',
		              },
		    'G->C' => {'+'=>'transcribe',
                               '-'=>'untranscribe',
		              },

		    'C->A' => {'+'=>'untranscribe',
                               '-'=>'transcribe',
		              },
                    'G->T' => {'+'=>'transcribe',
                               '-'=>'untranscribe',
		              },				   

		    'A->T' => {'+'=>'untranscribe',
                               '-'=>'transcribe',
			      },
                    'T->A' => {'+'=>'transcribe',
                               '-'=>'untranscribe',
		              },

                    'A->G' => {'+'=>'untranscribe',
                               '-'=>'transcribe',
			      },
                    'T->C' => {'+'=>'transcribe',
                               '-'=>'untranscribe',
		              },

                    'A->C' => {'+'=>'untranscribe',
                               '-'=>'transcribe',
			      },
                    'T->G' => {'+'=>'transcribe',
                               '-'=>'untranscribe',
		              },

                   };



    my $mutation_spectrum={};

    foreach my $SNV (keys %$SNVs) {
	my @strand = keys %{ $SNVs->{$SNV} };
	if(scalar(@strand) > 1) {
	    next;  #skip SNVs that intersect genes on both strands
	}
	my ($chr,$pos,$ref,$var) = split(/\_/,$SNV);
	my $k = join("->",($ref,$var));
	my $strand_type = $strand_rule->{$k}->{$strand[0]}; #returns 'transcribe' or 'untranscribe'
	if(!$strand_type) {
	    print "unknown strand info!\n";
	}
	my $category = return_mutation_spectrum_category($ref,$var);
	if(!defined($category)) {
	    print STDERR "unknown category!\n";
	}
	$mutation_spectrum->{$category}->{$strand_type}++;
	
	my $x=1;
    }


    my $MUTSPEC_FH = IO::File->new($plot_input_file,'w') or die "Unable to open the file $plot_input_file due to $!"; 

    foreach my $cat (keys %$mutation_spectrum) {
	foreach my $strand (keys %{ $mutation_spectrum->{$cat} }) {
	    print $MUTSPEC_FH "$cat\t$strand\t$mutation_spectrum->{$cat}->{$strand}\t$manual_label\n";
	}
	
    }
    $MUTSPEC_FH->close();


    if($plot_graph) {
	my $plot_cmd;
	$plot_cmd = qq{ plot_mutation_spectrum_bystrand(inputFile="$plot_input_file",outputFile="$plot_output_file",plot_title="$plot_title",num_row=$numberRow,file_width=$plot_width,file_height=$plot_height) };
	my $call = Genome::Model::Tools::R::CallR->create(command=>$plot_cmd, library=> "MutationSpectrum.R");
	$call->execute;
    }

    return 1;                               
}



sub intersectSNV2genes {

    my $SNV_file = shift;
    my $genefile = shift;


    my ($temp_FH, $tempfile) = Genome::Sys->create_temp_file;
    my $SNV_FH = IO::File->new($SNV_file,'r') or die "Unable to open the file $SNV_file due to $!"; 

    while(<$SNV_FH>) {
	chomp;
	next if(/chromo/);
	my @list = split(/\t/,$_);
	$list[1]--;
	my $str = join("\t",@list);
	print $temp_FH "$str\n";
    }
    $temp_FH->close;
    $SNV_FH->close;


    my ($temp_fh2, $tempfile2) = Genome::Sys->create_temp_file;
    my $cmd = "intersectBed -wao -a $tempfile -b $genefile | cut -f 1-9,11 > $tempfile2";

    system($cmd);

    #filter the SNV list to ones that intersect with gene(s)
    my $SNV_filtered = parse_IntersectBed_output($tempfile2);

    return $SNV_filtered;

    

}

sub parse_IntersectBed_output {

    my $file = shift;

    my $total=0;
    my $skip=0;
    my $intersect_gene=0;
    my $data={};
    
    open(FILE,$file) or die "Unable to open the file $file due to $!";
    while(my $line = <FILE>) {
	$total++;
	chomp $line;
	#my ($chr,$start,$stop,$ref_var,@rest) = split(/\t/,$_);
	my @list = split(/\t/,$line);
	if($list[5] eq '.' && $list[-1] eq '.') { #skip any SNV that does not intersect with any gene
	    $skip++;
	    next;
	}
	if($list[8] =~ /LOC\d+/ || $list[8] =~ /ENSG\d+/) { #skip any genes that start with LOC or ENSG
	    $skip++;
	    next;
	}

	$intersect_gene++;
	my $strand = $list[-1];
	my $key = join("_",@list[0,2,3,4]);
	$data->{$key}->{$strand}++;

    }
    close FILE;

    return $data;


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


sub parse_bed_file {

    my $file = shift;

    my $count = { 'A->C' => 0,
		  'A->G' => 0,
		  'A->T' => 0,
		  'C->A' => 0,
		  'C->G' => 0,
		  'C->T' => 0,
		  'G->A' => 0,
		  'G->C' => 0,
		  'G->T' => 0,
		  'T->A' => 0,
		  'T->C' => 0,
		  'T->G' => 0
    };

    open(BED,$file) or die "Unable to open the file $file due to $!";
    while(<BED>) {
	chomp;
	my ($chr,$start,$stop,$ref_var,@rest) = split(/\t/,$_);
	next if($chr =~ /GL/);	
	my ($ref,$var) = split(/\//,$ref_var);
	my @variants = Genome::Info::IUB::variant_alleles_for_iub($ref,$var);
	if(@variants>1) {
	    warn "more than 1 variant allele detected for '$_'\n";
	    next;
	}
	my $key = join("->",($ref,$variants[0]));
	$count->{$key}++;
    }
    close BED;

    return $count;


}



sub parse_grouping_file {

    my $file = shift;

    open(FILE, $file) or die "Can't open the file $file due to $!";
    
    my @files=();
    my @labels=();
    while(<FILE>) {
	chomp;
	next if(/^\s+$/); #remove empty lines
	my ($mut_spec_file,$label) = split(/\t/,$_);
	next if(!$mut_spec_file or !$label);
	$mut_spec_file = rem_white_space($mut_spec_file);
	$label =  rem_white_space($label);

	push(@files,$mut_spec_file);
	push(@labels,$label);
    }
    close FILE;

    return (\@files,\@labels);


}


sub rem_white_space {

    my $string = shift;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
