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

## really, this should all be written so that it parses the dv2 grammar properly, but a few regexen do a decent enough job for now. (even if they are pretty brittle)

    
    my $model=Genome::Model->get($self->model_id);
    my $pp = $model->processing_profile;

    #get a list of the 'slots' (alignment_strategy, snv_detection_strategy, etc)
    my @params = $pp->params_for_class;
    
    my %text;

    #for each of these, grab the value
    for my $param (@params){
        
        #punting on this stuff for now
        next if $param =~ /refcov|tiering|loh|output_plot|dnp_proportion|varscan_validation_version/;

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
            print $value . "\n\n";
        }
        # printf(
        #     "%s: %s\n",
        #     $param,
        #     Term::ANSIColor::colored(( defined $value ? $value : '<NULL>'), 'red'),
        #     );
    }
    
    print "\n";
    
    return 1;
}


sub cleanDescription{
    my ($param, $value, $model) = @_;
    
    #replace refseq build
    my $refseq_name = $model->reference_sequence_build->name;
    $value =~ s/reference_sequence_build/reference sequence build $refseq_name/g;

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


sub clean_alignment_strategy{
    my ($value) = @_;
    #filtered
    $value =~ s/using ([^\s]+) ([^\s]+)/using $1 version $2/g;    
    return($value)
}

sub clean_detection_strategy{
    my ($value) = @_;
    

    return($value);
}


sub addPrefix{
    my ($param,$value) = @_;

    if($param eq "alignment_strategy"){
        $value =~ s/^instrument_data aligned/Sequence data was aligned/g;
        $value = $value . ".";
    }

    if($param =~ /detection_strategy/){        
        $param =~ s/_detection_strategy//;
        $value = "We detected " . $param . "s using " . $value . ".";
        
        $value =~ s/union/unioned with/;
        $value =~ s/intersect/intersected with/;
        $value =~ s/unique union/unioned with/;
    }

    return $value;
}


sub getCitation{
    return 1;
}
