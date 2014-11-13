package Genome::Model::Tools::Analysis::SummarizeMutationSpectrum;

#####################################################################################################################################
# MutationSpectrum - Given a somatic variation model ID or modelgroup, will plot mutation spectrum barplot for each sample
#
#			
#####################################################################################################################################

use warnings;
use strict;
use Math::Round;
use Genome;
use Workflow;
use FileHandle;
use IO::File;
use Genome::Info::IUB;
use Cwd ('getcwd','abs_path');

class Genome::Model::Tools::Analysis::SummarizeMutationSpectrum {
    is => ['Command'],
    has => [
        input_file => {
            is => 'String',
            is_input => 1,
            is_optional => 1,
            doc => '.tsv file consiting of 2 entries per line: label and input file path. The input files should consist of 4-column, tab-separated bed files of type (chr start stop ref/var) or (chr start stop ref var)',
        },
        somatic_id => {
            is  => 'String',
            is_input=>1,
            is_optional => 1,
            doc => 'somatic variation id. If >1 id, comma-separate; ',
        },
        labels => {
            is  => 'String',
            is_input=>1,
            is_optional => 1,
            doc => 'explicit specifify sample label (normally ONLY used if user supplying a input file instead of using model id)',
        },
        group_id => {
            is  => 'String',
            is_input=>1,
            is_optional => 1,
            doc => 'model group id containing many somatic variation id.  Will supercede --somatic-id ',
        },
        plot_title => {
            is  => 'String',
            is_input=>1,
            is_optional => 1,
            default_value => 'Mutation Spectrums',
            doc => 'The title of the plot',
        },
        number_row => {
            is_input => 1,
            is_optional => 1,
            is => 'String',
            default => 'NULL',
            doc => 'for model_group/multiple samples, controls the number of rows the final plot will have (only valid if --make1plot is false)',
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
        make1plot => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'set this flag if you want multi-sample to be all plotted in 1 graph, default: samples are individually plotted',
            default => 0,
        },
        ymax => {
            is => 'float',
            is_optional => 1,
            default => 1.0,
            doc => 'Set the maximum Y axis: between 0 and 1',
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
            doc => 'Plot the user defined mutation spectrum file and exit.  ',
        },
        exclude_gl_contigs => {
            is_input => 1,
            is_optional => 1,
            is => 'Boolean',
            default => 1,
            doc => 'whether or not to exclude contigs beginning with the prefix GL.',
        },
    ],
};

sub help_brief {
    "Given lists of variants, produces an output of transition/transversion numbers and plots of these values";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"

gmt analysis summarize-mutation-spectrum --somatic-id=2880998721,2881022145 --mut-spec-file=mut_spec.file --output-file SJPHALL.pdf --number-row=2
gmt analysis summarize-mutation-spectrum --group-id=123456 --mut-spec-file=mut_spec.file --output-file SJPHALL.pdf --number-row=4
EOS
}

sub help_detail {
    return <<EOS
    This tool summarizes the mutation spectrum for a single/multiple models.  It produces a barplot for each mutation cateogry for each sample.  If user inputs multiple samples (i.e. via comma-separated somatic variation id or a model-group), the plot will contain miniplots for each sample.
EOS
}

sub execute {
    my $self = shift;
    $DB::single = 1;


    my $input_file = $self->input_file;
    my $somatic_id = $self->somatic_id;
    my $group_id = $self->group_id;
    my $plot_input_file;
    my $make1plot = $self->make1plot;
    my $ymax = $self->ymax;
    my $plot_title = $self->plot_title;
    my $numberRow = $self->number_row;
    my $manual_label = $self->labels;
    my $plot_graph = $self->plot_graph;

    if($self->mut_spec_file) {
        $plot_input_file = abs_path($self->mut_spec_file);
        unlink($plot_input_file) if(-e $plot_input_file); #remove existing file first
    }else {
        my ($fh, $tempfile) = Genome::Sys->create_temp_file;
        $plot_input_file = abs_path($tempfile);
    }
    my $plot_output_file = abs_path($self->output_file);

    if($self->plot_this) {
        my $user_mut_spec_file = abs_path($self->plot_this);
        my $plot_cmd;
        if(!$make1plot) {
            $plot_cmd = qq{ make_dodge_barplot_facet_sample(inputFile="$user_mut_spec_file",outputFile="$plot_output_file",plot_title="$plot_title",num_row=$numberRow,y_lim=c(0,$ymax)) };
        }
        else {
            $plot_cmd = qq{ make_dodge_barplot_sample(inputFile="$user_mut_spec_file",outputFile="$plot_output_file",plot_title="$plot_title",y_lim=c(0,$ymax)) };
        }
        my $call = Genome::Model::Tools::R::CallR->create(command=>$plot_cmd, library=> "MutationSpectrum.R");
        $call->execute;
        return 1;
    }


    if ($input_file){
        my $fh = IO::File->new($input_file, 'r');
        for my $line (<$fh>){
            chomp $line;
            my ($label, $file) = split("\t", $line);
            my $raw_count = $self->parse_bed_file($file);
            make_output($raw_count, $label, $plot_input_file);
        }
    }
    else{
        my @models=();
        if($group_id) {
            my $model_group = Genome::ModelGroup->get($group_id);
            @models = $model_group->models;
        } elsif($somatic_id) {
            my @modelIDs = split(/,/,$somatic_id);
            #@models = map{ Genome::Model->get($_) } @modelIDs;
            @models = Genome::Model->get(\@modelIDs);
        }

        foreach my $model(@models) {
            my ($input_file,$automatic_label) = make_input_file_from_model($model); #cat tier1,2,3 SNV bed file for each model
            my $raw_count = $self->parse_bed_file($input_file);
            unlink($input_file);
            if($manual_label) { #user specified a sample label
                make_output($raw_count,$manual_label,$plot_input_file);
            }
            else {
                make_output($raw_count,$automatic_label,$plot_input_file); #use automatically generated label for a sample
            }
        }
    }


    #my $input_plot_file = $out1;

    if($plot_graph) {
        my $plot_cmd;
        if(!$make1plot) {
            $plot_cmd = qq{ make_dodge_barplot_facet_sample(inputFile="$plot_input_file",outputFile="$plot_output_file",plot_title="$plot_title",num_row=$numberRow,y_lim=c(0,$ymax)) };
        }
        else {
            $plot_cmd = qq{ make_dodge_barplot_sample(inputFile="$plot_input_file",outputFile="$plot_output_file",plot_title="$plot_title",y_lim=c(0,$ymax)) };
        }
        my $call = Genome::Model::Tools::R::CallR->create(command=>$plot_cmd, library=> "MutationSpectrum.R");
        $call->execute;
    }

    return 1;
}


sub make_plots {



}



sub make_input_file_from_model {

    my $somatic_model = shift; #requires a somatic variation model object

    my $model_id = $somatic_model->id;
    die "Invalid somatic model for $model_id!  Aborting....\n" if(!defined($somatic_model));
    my $build=$somatic_model->last_succeeded_build;
    die "No successful build for model $model_id found! Aborting...\n"if(!defined($build));

    #Assign label in this order of preference if available: source common name -> sample name -> model id
    my $sample_label = $somatic_model->id;
    my $source_common_name = $build->tumor_build->model->subject->source_common_name;
    my $subject_name = $build->subject->name;
    $sample_label = $subject_name if ($subject_name);
    $sample_label = $source_common_name if ($source_common_name);

    #If available append the sample common name (e.g. tumor to the label)
    my $type = "_". $build->tumor_build->model->subject->common_name;
    $sample_label .= uc($type) if ($type);

    #find the tier 1,2,3 SNV bed file
    my $dir = $build->data_directory . "/effects";
    my $tier1 = "$dir/snvs.hq.novel.tier1.v2.bed";
    my $tier2 = "$dir/snvs.hq.novel.tier2.v2.bed";
    my $tier3 = "$dir/snvs.hq.novel.tier3.v2.bed";

    #these files should all exist for each build but just in case.....
    warn "$tier1 does not exist for $model_id." if(!-e $tier1);
    warn "$tier2 does not exist for $model_id." if(!-e $tier2);
    warn "$tier3 does not exist for $model_id." if(!-e $tier3);

    my ($fh,$temp_file) = Genome::Sys->create_temp_file;
    `cat $tier1 $tier2 $tier3 > $temp_file`;

    return ($temp_file,$sample_label);

}


sub parse_anno_file {

    my ($self, $file) = @_;

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


    open(ANNO,$file) or die "Can't open the file $file due to $!";
    while(<ANNO>) {
        chomp;
        my ($chr,$start,$stop,$ref,$var) = split(/\t/,$_);
        next if($chr =~ /^GL/ and $self->exclude_gl_contigs);
	my $key = join("->",($ref,$var));
	$count->{$key}++;
    }
    close ANNO;

    return $count;

}


sub parse_bed_file {

    my ($self, $file) = @_;

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

    my $fh = IO::File->new($file, 'r') or die "Unable to open $file due to $!";
    while(<$fh>) {
        chomp;
        my ($chr,$start,$stop,$ref_var,@rest) = split(/\t/,$_);
        next if($chr =~ /GL/);
        my ($ref,$var) = split(/\//,$ref_var);
        unless($var) {
          $var = $rest[0];
        }
        my @variants = Genome::Info::IUB::variant_alleles_for_iub($ref,$var);
        if(@variants>1) {
            #warn "more than 1 variant allele detected for '$_'\n";
            next;
        }
        my $key = join("->",($ref,$variants[0]));
        $count->{$key}++;
    }
    $fh->close;

    return $count;


}

sub make_output {

    my $raw_count = shift;
    my $label = shift;
    my $plotinput = shift;

    open (OUT, ">>", $plotinput) or die "Unable to write to $plotinput due to $!";

    my $transition=0;
    my $transversion=0;
    my $transition_percent=0;
    my $transversion_percent=0;
    my $total_SNV=0;
    my $total={};

    $total->{'A->C'}[0] = $raw_count->{'A->C'}+$raw_count->{'T->G'};
    $total->{'A->G'}[0] = $raw_count->{'A->G'}+$raw_count->{'T->C'};
    $total->{'A->T'}[0] = $raw_count->{'A->T'}+$raw_count->{'T->A'};
    $total->{'C->A'}[0] = $raw_count->{'C->A'}+$raw_count->{'G->T'};
    $total->{'C->G'}[0] = $raw_count->{'C->G'}+$raw_count->{'G->C'};
    $total->{'C->T'}[0] = $raw_count->{'C->T'}+$raw_count->{'G->A'};

    $transition = $total->{'A->G'}[0] + $total->{'C->T'}[0];
    $transversion = $total->{'A->C'}[0] + $total->{'A->T'}[0] + $total->{'C->A'}[0] + $total->{'C->G'}[0];
    $total_SNV = $transition+$transversion;
    if($total_SNV != 0) {
        $transition_percent = nearest(0.001,$transition/$total_SNV);
        $transversion_percent = nearest(0.001,$transversion/$total_SNV);
        #print OUT "BaseChg\tcount\tpercent\tlabel\n" if($print_header);
        foreach my $k (keys %$total) {
            $total->{$k}[1] = nearest(0.001,$total->{$k}[0]/$total_SNV);
            print OUT "$k\t$total->{$k}[0]\t$total->{$k}[1]\t$label\n";
        }
    } else {
        $transition_percent = 0;
        $transversion_percent = 0;
        #print OUT "BaseChg\tcount\tpercent\tlabel\n" if($print_header);
        foreach my $k (keys %$total) {
            $total->{$k}[1] = 0;
            print OUT "$k\t$total->{$k}[0]\t$total->{$k}[1]\t$label\n";
        }
    }

    print OUT "Transitions\t$transition\t$transition_percent\t$label\n";
    print OUT "Transversion\t$transversion\t$transversion_percent\t$label\n";
    close OUT;

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

1;
