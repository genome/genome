package Genome::Model::ClinSeq::Command::CreateWgsClonalityPlotByVariantSource;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::CreateWgsClonalityPlotByVariantSource {
    is => 'Command::V2',
    has_input => [
        build => { 
              is => 'Genome::Model::Build::ClinSeq',
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'clinseq build to get variant sources from',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
        callers => {
            is => 'String',
            doc => 'Names of the callers to graph',
            is_many => 1,
            valid_values => [qw(sniper strelka varscan)],
        },
    ],
    has_output => [
        _plots => {
              is => 'FilesystemPath',
              is_many => 1,
              is_optional =>1,
        },
        _source_file => {
            is => 'Hashref',
            is_optional => 1,
        },
    ],
    doc => 'plot clonality graphs for each of the sources of variants (i.e., which variant callers) for a ClinSeq build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq create-clonality-plot-by-variant-source --outdir=/tmp/  --callers strelka,sniper 128884819

EOS
}

sub help_detail {
    return <<EOS
Graph Clonality plots for each combination of source of variants (i.e., snv/indel caller) for ClinSeq build

(put more content here)
EOS
}

sub execute {
    my $self = shift;
    #grab the variant source file
    
    my $variant_source_file;
    eval {
        $variant_source_file = $self->_snv_variant_source_file($self->build, "wgs"); #this will die if wgs doesn't exist
    };
    if($@) {
        $self->error_message("Unable to find SNV source file for WGS data. This tool will not run.");
        $self->error_message($@->what);
        return 1; #this seems a little bad, but we don't want workflows to fail if no WGS data is available. Could be handled in the workflow, I suppose.
    }

    my $varscan_readcount_file = $self->_varscan_formatted_readcount_file($self->build);

    my $cnaseq_hmm_file = $self->_cnaseq_hmm_file($self->build);

    $self->_load_source_file_into_memory($variant_source_file);

    $DB::single = 1;
    my @caller_combinations = combination_of_callers($self->callers);
    while(my $callers = shift @caller_combinations) {
        $self->_generate_clonality_for_callers($varscan_readcount_file, $cnaseq_hmm_file, $callers);
    }
    return 1;
}



#most likely the number of calls is much smaller than we can fit in memory so just do that with a hash.
sub _load_source_file_into_memory {
    my ($self, $snv_variant_source_file) = @_;
    #  coord	chr	start	end	variant	score1	score2	callers	strelka	sniper	varscan	tier
    #  10:100062111-100062112	10	100062111	100062112	C/M	149	31	sniper,strelka,varscan	1	1	1	tier3
    #  10:102501519-102501520	10	102501519	102501520	C/G	33	34	sniper	0	1	0	tier2
    #  10:102766573-102766574	10	102766573	102766574	A/G	18	7	sniper	0	1	0	tier1

    my %source_file_line;
    my $fh = Genome::Sys->open_file_for_reading($snv_variant_source_file);
    my @header_fields = split /\t/, $fh->getline;   #parse header line
    while(my $line = $fh->getline) {
        chomp $line;
        my %cells;
        @cells{@header_fields} = split /\t/, $line;
        my ($ref, $genotype) = split /\//, $cells{variant};
        my $readcount_string = join("\t", @cells{qw(chr end)}, $ref, $genotype);
        my @callers = split /,/, $cells{callers};
        $source_file_line{$readcount_string} = \@callers;
    }
    $fh->close;

    $self->_source_file(\%source_file_line);
}

sub _create_readcount_file_for_callers {
    my ($self, $input_readcount_file, $callers, $output_file) = @_;

    unless(defined $self->_source_file) {
        die $self->error_message("Caller source file not loaded before trying to dump out readcount file");
    }

    my %caller_set = map {$_ => 1} @$callers;
    my %retain_line;
    for my $key (keys %{$self->_source_file}) {
        if( grep { exists($caller_set{$_}) } @{ $self->_source_file->{$key} } ) {
            $retain_line{$key} = 1;
        }
    }
    
    my $ofh = Genome::Sys->open_file_for_overwriting($output_file);
    my $fh = Genome::Sys->open_file_for_reading($input_readcount_file);
    while(my $line = $fh->getline) {
        my @fields = split /\t/, $line, 5;
        my $string = join("\t", @fields[0..3]);
        if(exists($retain_line{$string})) {
            print $ofh $line;
        }
    }
    $ofh->close;
    $fh->close;
}

sub _generate_clonality_for_callers {
    my ($self, $varscan_readcount_file, $cnaseq_hmm_file, $callers) = @_;
    my $file_suffix = join("_",@$callers);
    my $readcount_output_file = $self->outdir . "/varscan_formatted_readcount.$file_suffix";
    my $image_file = $self->outdir . "/${file_suffix}_calls.pdf";
    my $rscript_output = $self->outdir . "/${file_suffix}_calls.R";
    my $sample_id = $self->build->common_name;
    my $graph_title = "$sample_id " . join(" U ", @$callers) . " calls";

    $self->_create_readcount_file_for_callers($varscan_readcount_file, $callers, $readcount_output_file);

    my $retval = Genome::Model::Tools::Validation::ClonalityPlot->execute( cnvhmm_file => $cnaseq_hmm_file, output_image => $image_file, r_script_output_file => $rscript_output, varscan_file => $readcount_output_file, analysis_type => "wgs", sample_id => $graph_title); 
    unless($retval) {
        die $self->error_message("Error generating clonality plot.");
    }
}





        
#this is probably overkill below, but recursion is nice
sub combination_of_callers {
    my @array = @_;
    my @list = ([@array]);
    while( my $top_element = shift @array ) {
        push @list, [$top_element];
        for my $element (@array) {
            push @list, [$top_element, $element];
        }
    }
    return @list;
}





1;
