package Genome::Model::Tools::Annotate::LookupConservationScore;

use strict;
use warnings;
use Genome;
use IO::File;
use Genome::Info::UCSCConservation;

class Genome::Model::Tools::Annotate::LookupConservationScore{
    is => 'Genome::Model::Tools::Annotate',
    has => [
        chromosome => {
            is => 'Text',
            is_input => 1,
            doc => "Name of the chromosome where the coordinates reside", 
            is_optional => 0,
        },
        coordinates => {
            is => 'Text',
            is_input => 1,
            doc => "Coordinates in the chromosome to check.  Should come in the form of a string or array of positions (ex: (1, 2, 3, 4, 5))", 
            is_optional => 0,
        },
        species => {
            is => 'Text',
            is_input => 1,
            doc => "Species used for scoring.  Species and version must be specified unless reference_transcripts is specified", 
            is_optional => 1,
        },
        version=> {
            is => 'Text',
            is_input => 1,
            doc => "Version of the species used for scoring.  Must be in a form like 54_36p.  Species and version must be specified unless reference_transcripts is specified", 
            is_optional => 1,
        },
        reference_transcripts =>{
            is => 'Text',
            is_input => 1,
            doc => "model name and version used for scoring (ex: NCBI-human.combined-annotation/54_36p).  reference_transcripts must be specified unless species and version are specified", 
            is_optional => 1,
        },
        conservation_scores_results => {
            is => 'SCALAR', 
            is_output => 1,
            doc => "This is populated with an array reference containing the results of the conservation score lookup post execute",
            is_optional => 1,
        }
    ],
};

sub help_synopsis {
    return <<EOS
my \$lookup = Genome::Model::Tools::Annotate::LookupConservationScore->execute(chromosome => 1, coordinates => '(1, 2, 3)', species => 'human', version => '54_36p_v2');
my \$results_array_ref = \$lookup->conservation_scores_results;
EOS
}

# For speed.  slightly unsafe, but there shouldn't ever be any problems, right?
sub __errors__() {
    return;
}

sub execute {

    my $self = shift; 

    my $results = 
        eval { lookup_conservation_score(
                   chromosome => $self->chromosome,
                   coordinates => $self->coordinates,
                   species => $self->species,
                   version => $self->version,
                   reference_transcripts => $self->reference_transcripts)
        };
    if ($@) {
        $self->error_message($@);
        return 0;
    } else {
        $self->conservation_scores_results($results);
        return 1;
    }
}


sub lookup_conservation_score {
    my %args = @_;
    my($chromosome, $coordinates, $species, $version, $reference_transcripts) = @args{'chromosome','coordinates','species','version','reference_transcripts'};

    unless(ref($coordinates) eq 'ARRAY'){
        Carp::croak("Malformed coordinates.  Coordinates should be specified as an array ref");
        return 0;
    }
    
    my ($model_name, $reference_version, $combined_annotation_build_version);
    if($reference_transcripts and $species and $version){
        Carp::croak("species and version OR reference_transcripts must be specified");
        return 0;
    } elsif($reference_transcripts){
        ($model_name, $combined_annotation_build_version) = split("/", $reference_transcripts); 
        $reference_version = $combined_annotation_build_version;
        $reference_version =~ m/[\d]+_([\d]+)[a-zA-Z]/ or die "Malformed reference_transcripts: $model_name/$reference_version";
        $version = $1;
    } elsif($species and $version){
        $model_name = "NCBI-" . $species . ".combined-annotation";
        if($version =~ m/[\d]+_([\d]+)[a-zA-Z]/){
            $version = $1; 
            $reference_version = $version;
        }else{
            Carp::croak("Malformed version: " . $version);
            return 0;
        }
    }else{
        Carp::croak("species and version OR reference_transcripts must be specified");
        return 0;
    }
    
    my $location;  
    
    if($combined_annotation_build_version){
        #check if the combined annotation build has UCSC conservation data.  use it if it exists
        my $combined_annotation_model = Genome::Model->get(name => $model_name);
        my $combined_annotation_build = $combined_annotation_model->build_by_version($combined_annotation_build_version);
        $location = $combined_annotation_build->ucsc_conservation_directory if (-d $combined_annotation_build->ucsc_conservation_directory);
    }

    unless($location){
        #combined annotation build doesn't have a UCSC conservation dir.  Try the old logic
        #TODO: Remove these restrictions as soon as possible
        unless($model_name eq 'NCBI-human.combined-annotation'){
            Carp::croak("Cannot look up conservation scores for model: $model_name.  This module currently supports the following models: NCBI-human.combined-annotation"); 
            return 0;
        }
        unless($version == 36 or $version == 37){
            Carp::croak("Cannot look up conservation scores for model $model_name version $version.  This module currently supports the following versions for $model_name: 36, 37"); 
            return 0;
        }

        my %location_index = Genome::Info::UCSCConservation->ucsc_conservation_directories;
        $location = $location_index{$version};
    }

    unless ($location){
        Carp::croak("Cannot find a UCSC conservation directory for model_name $model_name version $version"); 
        return 0;
    }

    #TODO: we're copying directly from MG/ConsScore.pm from this point on.  It probably sucks and
    #needs to be refactored
    # open file, 
    our %file_handles;
    my $file = $location . "/chr".$chromosome."-rec";
    unless (exists($file_handles{$file} ) ) {
        $file_handles{$file} = IO::File->new($file, 'r');
        unless ($file_handles{$file}) {
            die "can't open $file : $!";
        }
    }
    my $fh = $file_handles{$file};

    my @results;
    @$coordinates = sort { $a <=> $b } @$coordinates;
    foreach my $pos (@{$coordinates})
    {
        # start seeking until the first position
        my $seekpos = ($pos -1)*2;
        unless ($fh->seek($seekpos,0) ) {
            Carp::croak("seek to pos $seekpos of $file failed: $!");
        }

        my $tmpval;
        my $read_result = $fh->sysread($tmpval,2);
        if (! defined ($read_result) ) {
            Carp::croak("Reading from pos $seekpos of $file failed: $!");
        } elsif (!$read_result) {
            Carp::croak("Attempt to read from pos $seekpos of $file resulted in EOF");
        } elsif ($read_result != 2) {
            Carp::croak("Short read from pos $seekpos of $file.  Asked for 2 bytes, got $read_result");
        }

        # unpack value, divide by 1000 and push onto array.
        my $score = unpack("n",$tmpval)/1000;
        push(@results, [ $pos, $score ]);
    }
    # return array ref of array refs of positions, scores?
    #$fh->close;
    return \@results;
}

1;
