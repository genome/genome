package Genome::Model::Tools::Analysis::DumpIgvXmlMulti;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Cwd qw(abs_path);

class Genome::Model::Tools::Analysis::DumpIgvXmlMulti {
    is => 'Command',
    has => [
    bams => {
        type => 'String',
        is_optional=>0,
        doc => "Bam files (single bam or comma separated list)",
    },
    labels => {
        type => 'String',
        is_optional=>0,
        doc => "labels for the data (single label or comma separated list)",,
    },
    review_bed_file => {
        type => 'String',
        is_optional => 0,
        doc => "Bed file of sites to review",
    },
    genome_name => {
        type => 'String',
        is_optional => 0,
        doc => "String with which to name the file and tracks eg AML04",
    },
    review_description => {
        type => 'String',
        is_optional => 1,
        default => "Sites to review",
        doc => "description to be displayed in IGV of the sites to review eg Tier1 Indels",
    },
    output_file => {
        type => 'String',
        is_optional => 0,
        doc => 'Output XML file',
    },
    reference_name => {
        type => 'String',
        is_optional => 0,
        doc => 'the name of the reference (in IGV) that the bams are aligned to. E.g. "b37" for build 37 or "reference"for our internal build36',
    },
    roi => {
        is => 'Text',
        is_optional => 1,
        doc => 'A bed file with the regions of interest to load into the session',
    },
    ]
};

sub help_brief {
    "Makes an IGV session for review - takes more than just tumor/normal pairs - supports as many tracks as you want"
}

sub help_detail {
    <<'HELP';
This helps create files for manual reviewers
HELP
}

sub execute {
    my $self=shift;
    my $tumor_common_name = $self->genome_name;   
    my $output_file = $self->output_file;
    my $genome_name = $self->genome_name;
    my $review_bed_file = abs_path($self->review_bed_file);
    my $reference_name = $self->reference_name;

    my @bams = split(/\,/,$self->bams);
    my @labels = split(/\,/,$self->labels);

    unless (@bams == @labels){
        die("number of bams and number of labels must be equal");
    }



#---- header -----
    my $header = <<"XML";
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="$reference_name" locus="1:1-100" version="3">
    <Resources>
    <Resource path="$review_bed_file" relativePath="false"/>
XML

    my $i=0;    
    for($i=0;$i<@bams;$i++){
        my $path = abs_path($bams[$i]);
        $header .= <<"XML";
     <Resource path="$path" relativePath="false"/>
XML
}
    $header .= "    </Resources>\n";

#---- regions ----
my $regions;
    if ($self->roi) {
        $regions = "<Regions>\n";
        my $in = Genome::Sys->open_file_for_reading($self->roi);
        while (my $line = <$in>) {
            chomp $line;
            my @fields = split(/\t/, $line);
            $regions .= <<"XML";
            <Region chromosome="$fields[0]" end="$fields[2]" start="$fields[1]"/>
XML
        }
        $regions .= "</Regions>\n";
    }
#---- panels -----

    my $panels;

    for($i=0;$i<@bams;$i++){
        my $path = abs_path($bams[$i]);
        my $label = $labels[$i];
        my $cov = $path . "_coverage";

$panels .= <<"XML";        
   <Panel height="1000" name="Panel$i" width="1901">
        <Track color="200,200,200" expand="false" fontSize="9" height="40" id="$cov" name="$label Coverage" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track expand="true" fontSize="9" height="1000" id="$path" name="$label" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
   </Panel>
XML
    }
#---- panels -----
my $features = <<"XML"; 
    <Panel height="60" name="FeaturePanel" width="1901">
        <Track color="0,0,0" expand="false" fontSize="9" height="14" id="Reference" name="Reference" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track color="0,0,178" expand="true" featureVisibilityWindow="-1" fontSize="9" height="45" id="$review_bed_file" name="$review_bed_file" renderer="BASIC_FEATURE" showDataRange="true" visible="true" windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
XML


#----- do output ----
    open(OUTFILE,">$output_file") || die("could not open output file $output_file");
    print OUTFILE $header;
    print OUTFILE $panels;
    print OUTFILE $features;
    if ($regions) {
        print OUTFILE $regions;
    }
    
    print OUTFILE "</Session>\n";
    close(OUTFILE);
    
return 1;
    
}


1;
