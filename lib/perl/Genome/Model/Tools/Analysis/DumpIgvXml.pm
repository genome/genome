package Genome::Model::Tools::Analysis::DumpIgvXml;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Cwd qw(abs_path);

class Genome::Model::Tools::Analysis::DumpIgvXml {
    is => 'Command',
    has => [
    tumor_bam => {
        type => 'String',
        is_optional=>0,
        doc => "Bam file for the tumor data",
    },
    normal_bam => {
        type => 'String',
        is_optional=>0,
        doc => "Bam file for the normal data",
    },
    relapse_bam => {
        type => 'String',
        is_optional => 1,
        doc => "Bam file for any relapse data",
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
    output_dir => {
        type => 'String',
        is_optional => 0,
        doc => 'directory to dump session files. Will be named by GSC common name',
    },
    reference_name => {
        type => 'String',
        is_optional => 0,
        default => 'reference',
        doc => 'the name of the reference (in IGV) that the bams are aligned to. E.g. b37 for build 37 or reference for our internal build36',
    },
    ]
};


sub execute {
    my $self=shift;
    my $tumor_common_name = $self->genome_name;   
    my $output_dir = $self->output_dir;
    my $ofh = IO::File->new("$output_dir/$tumor_common_name.xml","w");
    unless($ofh) {
        $self->error_message("Unable to open $output_dir/$tumor_common_name.xml for writing");
        return;
    }

    if($self->relapse_bam) {
        print $ofh $self->generate_xml($self->reference_name,"$tumor_common_name Tumor",abs_path($self->tumor_bam), "$tumor_common_name Normal", abs_path($self->normal_bam), abs_path($self->review_bed_file), $self->review_description, "$tumor_common_name Relapse", abs_path($self->relapse_bam));
    }
    else {
        print $ofh $self->generate_xml($self->reference_name,"$tumor_common_name Tumor",abs_path($self->tumor_bam), "$tumor_common_name Normal", abs_path($self->normal_bam), abs_path($self->review_bed_file), $self->review_description);
    }
    $ofh->close;
    return 1;

}


1;

sub help_brief {
    "Makes an IGV session for review"
}

sub help_detail {
    <<'HELP';
This helps create files for manual reviewers
HELP
}

sub generate_xml {
    my ($self, $reference_name, $tumor_name, $tumor_bam, $normal_name, $normal_bam, $review_bed_file, $review_bed_name,$relapse_name, $relapse_bam) = @_;
    my $relapse_resource = "";
    if($relapse_bam) {
        $relapse_resource = qq{\n        <Resource path="$relapse_bam" relativePath="false"/>};
    }

    my $relapse_panel = "";
    if($relapse_bam) {
        $relapse_panel = <<"RELAPSE_PANEL";
    <Panel height="2040" name="Panel128767388888" width="1901">
        <Track color="200,200,200" expand="false" fontSize="9" height="40" id="${relapse_bam}_coverage" name="$relapse_name Coverage" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track expand="true" fontSize="9" height="2000" id="$relapse_bam" name="$relapse_name" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
RELAPSE_PANEL
    }   
    my $xml = <<"XML";
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="$reference_name" locus="1:1-100" version="3">
    <Resources>
        <Resource path="$tumor_bam" relativePath="false"/>
        <Resource path="$review_bed_file" relativePath="false"/>
        <Resource path="$normal_bam" relativePath="false"/>$relapse_resource
    </Resources>
$relapse_panel   <Panel height="2040" name="Panel1287673856180" width="1901">
        <Track color="200,200,200" expand="false" fontSize="9" height="40" id="${tumor_bam}_coverage" name="$tumor_name Coverage" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track expand="true" fontSize="9" height="2000" id="$tumor_bam" name="$tumor_name" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
    <Panel height="2040" name="Panel1287673878899" width="1901">
        <Track color="200,200,200" expand="false" fontSize="9" height="40" id="${normal_bam}_coverage" name="$normal_name Coverage" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track expand="true" fontSize="9" height="2000" id="$normal_bam" name="$normal_name" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
    <Panel height="60" name="FeaturePanel" width="1901">
        <Track color="0,0,0" expand="false" fontSize="9" height="14" id="Reference" name="Reference" showDataRange="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track color="0,0,178" expand="true" featureVisibilityWindow="-1" fontSize="9" height="45" id="$review_bed_file" name="$review_bed_name" renderer="BASIC_FEATURE" showDataRange="true" visible="true" windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
    </Panel>
</Session>
XML

}
