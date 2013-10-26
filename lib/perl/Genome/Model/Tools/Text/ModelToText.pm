package Genome::Model::Tools::Text::ModelToText;

use strict;
use Genome;
use warnings;

class Genome::Model::Tools::Text::ModelToText{
    is => 'Command',
    has => [

        # processing_profile_id => {
	#     is => 'String',
	#     is_optional => 0,
	#     doc => 'id of the processing profile to transform to text',
	# },

        model_id => {
            is => 'String',
            is_optional => 0,
            doc => "id of the model to describe",
        }
        ]
};


sub help_brief {
    "Turn a processing profile into human-readable text, including citations"
}

sub help_detail {
    "Turn a processing profile into human-readable text, including citations"
}

#########################################################################

sub execute {
    my $self = shift;

    my @output;
    ## really, this should all be written so that it parses the dv2 grammar properly
    ## but a few regexen do a decent enough job for now. (even if they are pretty brittle)
    ##


    my $model=Genome::Model->get($self->model_id);
    my $pp = $model->processing_profile;

    print STDERR $pp->id . "\n";

    #get a list of the 'slots' (alignment_strategy, snv_detection_strategy, etc)
    my @params = $pp->params_for_class;

    if($model->subclass_name =~ /SomaticVariation/){
        my $refpp = $model->tumor_model->processing_profile;
        my $text = gatherRefAlignCommands($model->tumor_model,$refpp);
        push(@output,$text);
    }

    #for each of these, grab the value
    for my $param (@params){

        #punting on this stuff for now
        next if $param =~ /refcov|tiering|loh|output_plot|dnp_proportion|varscan_validation_version|transcript_variant_annotator_version|filter_previously_discovered_variants|vcf_annotate_dbsnp_info_field_string/;

        my @values = $pp->$param;

        foreach my $value (@values) {
            if(defined($value)){
                if (Scalar::Util::blessed($value)) {
                    $value = $value->__display_name__;
                }

                $value = cleanDescription($param, $value, $model);
            }
        }


        my $value;
        if (@values != 0 and $values[0]) {
            $value = join(",",@values);
        }

        if(defined($value)){
            push(@output, $value);
        }
    }

    my %citations = hashCitations();

    my $ref = 1;
    my @refs;
    my %refhash;

    my @output2;
    foreach my $line (@output){
        my $newline;
        my @tmp = split("version",$line);
        
        #mark each tool with a reference
        while ($line =~ /[^filtered by] ([^\s]+) version ([^\s]+)/g) {            
            unless(defined($newline)){
                $newline = $line;
            }

            my $tool = $1;
            my $vers = $2;
            #reuse citation numbers if the same tool appears multiple times
            my $refins = "[$ref]";
            if(defined($refhash{$tool})){
                $refins = "[" . $refhash{$tool} . "]";
            } else {
                if(defined($citations{$tool})){
                    push(@refs,$citations{$tool});
                } else {
                    push(@refs,"$tool CITATION NEEDED");
                }
            }

            $newline =~ s/$tool version $vers/$tool version $vers $refins/;

            unless(defined($refhash{$tool})){
                $refhash{$tool} = $ref;
                $ref++;
            }

        }
        push(@output2,$newline);
    }

    print "\n" . join("\n\n", @output2) . "\n";
    
    my $i;
    for($i=0;$i<@refs;$i++){
        print "[" . ($i+1) . "] " . $refs[$i] . "\n";
    }

    return 1;
}

#----------------------------------------
sub cleanDescription{
    my ($param, $value, $model) = @_;

    #replace refseq build
    if($value =~ /reference_sequence_build/){
        my $refseq_name = $model->reference_sequence_build->name;
        $value =~ s/reference_sequence_build/reference sequence build $refseq_name/g;
    }
    
    #replace parentheses - probably should try to disambiguate these, but not going to bother for now
    $value =~ s/\(//g;
    $value =~ s/\)//g;


    if($param eq "alignment_strategy"){
        $value = clean_alignment_strategy($value)
    }

    #fix param descriptions
    $value =~ s/\[([^\]]+)\]/\(params: $1\)/g;

    if($param =~ 'detection_strategy'){
        $value = clean_detection_strategy($value)
    }

    $value = addPrefix($param,$value);

    return($value);

}

#----------------------------------------
sub clean_alignment_strategy{
    my ($value) = @_;
    #filtered
    $value =~ s/using ([^\s]+) ([^\s]+)/using $1 version $2/g;
    return($value)
}

#----------------------------------------
sub clean_detection_strategy{
    my ($value, %citations) = @_;

    $value = splitAndFix($value);
    $value =~ s/zzversion/version/g;
        

    $value =~ s/union/unioned with/;
    $value =~ s/intersect/intersected with/;
    $value =~ s/unique union/unioned with/;

    return($value);
}

#----------------------------------------
sub splitAndFix{
    my ($value) = @_;

    #strip whitespace
    $value =~ s/^\s+//;
    $value =~ s/\s+$//;
    #print STDERR "\nVALUE: " . $value . "\n";

    #finished
    if($value =~ /zzversion/){
        return($value);
    }

    #base case - add version
    unless($value =~ /union unique|union|intersect|filtered by|then/){
        my @vals = split(/ /,$value);
        #print STDERR "SPLIT: " . join("|",@vals) . "\n";

        my $v = join(" ",($vals[0],"zzversion",@vals[1..$#vals]));
        #print STDERR "BASE: " . $v . "\n";
        return($v)
    }

    for my $joiner ("union unique","union","intersect","filtered by","then"){
        if($value =~ /$joiner/){            
            my @vals = split($joiner,$value);
            my $newval;
            foreach my $val (@vals){
                if(defined($newval)){
                    $newval = $newval . " $joiner " . splitAndFix($val);
                } else {
                    $newval = splitAndFix($val);
                }
            }
            $value = $newval;
        }
    }
    return($value);
}


#----------------------------------------
sub addPrefix{
    my ($param,$value) = @_;

    if($param eq "alignment_strategy"){
        $value =~ s/^instrument_data aligned/Sequence data was aligned/g;
        $value = $value . ".";
    }

    if($param =~ /detection_strategy/){
        $param =~ s/_detection_strategy//;
        $value = "We detected " . $param . "s using " . $value . ".";
    }

    return $value;
}

#----------------------------------------
sub gatherRefAlignCommands{
    my ($model,$pp) = @_;
    my $text = "Sequence data was aligned to reference sequence build ";
    $text .= $model->reference_sequence_build->name;
    $text .= " using " . $pp->read_aligner_name;
    $text .= " version " . $pp->read_aligner_version;
    $text .= " (params: " . $pp->read_aligner_params . ")";
    $text .= " then merged using " . $pp->merger_name;
    $text .= " version " . $pp->merger_version;
    $text .= " then deduplicated using " . $pp->duplication_handler_name;
    $text .= " version " . $pp->duplication_handler_version;
    $text .= ".";
    return($text);
}

#------------------------------------------------
sub hashCitations{
    my %citations; 
    $citations{"bwa"} = "Li H. and Durbin R.Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. (2009)";
    $citations{"samtools"} = "Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079. (2009)";
    $citations{"pindel"} =  "Ye, K., Schulz, M. H., Long, Q., Apweiler,R. & Ning, Z. Pindel: a pattern growth approach to detect breakpoints of large deletions and medium sized insertions from paired-end short reads.Bioinformatics 25, 2865-2871. (2009)";
    $citations{"sniper"} = "Larson, D.E. et al. SomaticSniper: Identification of Somatic Point Mutations in Whole Genome Sequencing Data. Bioinformatics. (2011)";
    $citations{"varscan"} = "Koboldt D.C. et al. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Research, 22(3):568-76. (2012)";    
    $citations{"gatk"} = "McKenna, A. et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20, 1297-1303. (2010).";
    $citations{"breakdancer"} = "Chen, K. et al. BreakDancer: an algorithm for high-resolution mapping of genomic structural variation. Nat Methods 6, 677-681. (2009)";
    $citations{"squaredancer"} = "unpublished - available at https://github.com/genome/gms-core/blob/master/lib/perl/Genome/Model/Tools/Sv/SquareDancer.pl";
    $citations{"tigra-sv"} = "in preparation";
    $citations{"strelka"} = "Saunders, C.T., Wong, W., Swamy, S., et al. Strelka: Accurate somatic small-variant calling from sequenced tumor-normal sample pairs. Bioinformatics (2012).";
    $citations{"music"} = "Dees, N. D. et al. MuSiC: Identifying mutational significance in cancer genomes. Genome Res 22, 1589-1598. (2012)";
    return(%citations);
}
