package Genome::File::IgvXml;

use strict;
use warnings;
use Genome;
use Carp;

use Genome::Utility::List;

our $VERSION = 4; #cribbing code from Malachi's G::M::ClinSeq::Comand::DumpIgvXml and its version is 4

class Genome::File::IgvXml {
    is => 'Genome::File::Base',
    has_optional => [
        _panel_count => {
            is => 'Number',
            default => 0,
            doc => 'The number of panels in the file. Used for layout.',
        },
        _resources => {
            is => 'Array',
            default => [],
            doc => 'IGV resource declarations',
        },
        _data_tracks => {
            is => 'Array',
            default => [],
            doc => 'Main features for IGV session',
        },
        _feature_tracks => {
            is => 'Array',
            default => [],
            doc => 'Extra features at the bottom (e.g. Genes, BED files etc)',
        },
        base_url => {
            is => 'String',
            default => $ENV{GENOME_SYS_SERVICES_FILES_URL},
            doc => 'base url for all file paths. Will be prepended.',
        },
        font_size => {
            is => 'Number',
            default => 11,
            doc => 'font size for all tracks',
        },
        color_option => {
            is => 'String',
            default => 'UNEXPECTED_PAIR',
            doc => 'Default coloring option for BAM tracks',
        },
        genome => {
            is => 'String',
            valid_values => [ qw( b37 ) ],
            is_optional => 0,
            doc => 'The name of the reference sequence in IGV',
        },
        ],

};

sub add_resource {
    my ($self, $resource) = @_;
    chomp $resource;
    push @{$self->_resources}, "    <Resource path=\"$resource\"/>";
}

sub add_data_track {
    my ($self, $xml) = @_;
    chomp $xml;
    push @{$self->_data_tracks}, $xml;
}

sub add_feature_track {
    my ($self, $xml) = @_;
    chomp $xml;
    push @{$self->_feature_tracks}, $xml;
}

sub add_bam_track {
    my ($self,%options) = @_;
    #check that minimum parameters have been passed in
    unless(defined $options{file}) {
        croak "No name passed to add_bam_track\n";
    }
    $options{name} ||= $options{file};
    my $coverage_track_name = "$options{name} Coverage";
    my $read_track_name = "$options{name} Reads";
    my $resource_file_coverage_url = Genome::Utility::List::join_with_single_slash($self->base_url, $options{file} . "_coverage");
    my $resource_file_url = Genome::Utility::List::join_with_single_slash($self->base_url, $options{file});
    my $font_size = $options{font_size} || $self->font_size;
    my $color_option = $options{color_option} || $self->color_option;
    my $max = defined $options{max_depth} ? $options{max_depth} : 100;
    my $min = defined $options{min_depth} ? $options{min_depth} : 0;

    my $xml=<<XML;
    <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;9062.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" id="$resource_file_coverage_url" name="$coverage_track_name" showDataRange="true" visible="true">
      <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="$max" minimum="$min" type="LINEAR"/>
    </Track>
    <Track altColor="0,0,178" color="0,0,178" colorOption="$color_option" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="$font_size" id="$resource_file_url" name="$read_track_name" showDataRange="true" sortByTag="" visible="true"/>
XML
    $self->add_resource($resource_file_url);
    $self->add_data_track($self->make_panel($xml));
}

sub add_bed_track {
    my ($self,%options) = @_;
    #check that minimum parameters have been passed in
    unless(defined $options{file}) {
        croak "No name passed to add_bam_track\n";
    }
    $options{name} ||= $options{file};
    my $resource_file_url = Genome::Utility::List::join_with_single_slash($self->base_url, $options{file});
    my $font_size = $options{font_size} || $self->font_size;
    my $bed_track_name = $options{name};
    my $xml=<<XML;
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" height="45" id="$resource_file_url" name="$bed_track_name" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
XML
    $self->add_resource($resource_file_url);
    $self->add_feature_track($xml);
}

sub add_junction_track {
    my ($self,%options) = @_;
    #check that minimum parameters have been passed in
    unless(defined $options{file}) {
        croak "No name passed to add_bam_track\n";
    }
    $options{name} ||= $options{file};
    my $resource_file_url = Genome::Utility::List::join_with_single_slash($self->base_url, $options{file});
    my $font_size = $options{font_size} || $self->font_size;
    my $read_track_name = $options{name};

    my $xml=<<XML;
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" height="60" id="$resource_file_url" name="$read_track_name" showDataRange="true" visible="true" windowFunction="count"/>
XML

    $self->add_resource($resource_file_url);
    $self->add_feature_track($xml);
}

sub make_gene_track {
    my ($self) = @_;
    my $gene_track_name = $self->genome . " Genes";
    my $font_size = $self->font_size;
    my $xml = <<XML;
    <Track altColor="0,0,178" color="0,0,178" colorScale="ContinuousColorScale;0.0;160.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" height="35" id="$gene_track_name" name="Gene" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
XML
    return $xml;
}

sub make_panel {
    my ($self, $xml, $panel_name) = @_;
    $self->_panel_count($self->_panel_count + 1);
    my $name = $panel_name;
    $name ||= "Panel " . $self->_panel_count;
    chomp $xml;
    my $panel_xml =<<XML;
  <Panel name="$name">
$xml
  </Panel>
XML
    return $panel_xml;
}

sub generate_resource_xml {
    my ($self) = @_;

    unless( @{$self->_resources} > 0) {
        $self->error_message("No resources found for creation of a resource XML section for IGV");
    }

    my $xml = join("\n",("  <Resources>", sort @{$self->_resources}, "  </Resources>"));
    return $xml."\n";
}

sub generate_feature_panel_xml {
    my ($self) = @_;

    my $font_size = $self->font_size;
    my $reference_track = <<REFXML;
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>
REFXML
    my $gene_track = $self->make_gene_track();

    my $xml = $self->make_panel(join(q{},$reference_track,$gene_track,@{$self->_feature_tracks}), "FeaturePanel");
    return $xml;
}

sub generate_panel_layout_xml{
    my ($self, %args) = @_;
    my $panel_count = $self->_panel_count;

    my $xml;

    #Every session should have at least 2 panels (the feature panel + at least one data panel) otherwise there may be an input problem and/or the igv session would 
    #If for some reason only the feature panel was present, simply return and empty string for the panel layout as a layout would be uneccessary
    if ($panel_count == 0) {
        $self->error_message("No panels defined for creation panel layout XML");
        return;
    } elsif($panel_count == 1) {
        $xml = '';
    } else {
        #If an IGV session has say 6 panels (including the feature panel), it will have 6 entries in this list
        #The first entry always seems to be near 0.  The last entry is 1 - (the desired relative space used by the feature panel)
        # <PanelLayout dividerFractions="0.007352941176470588,0.32107843137254904,0.6017156862745098,0.8762254901960784"/>
        my @divs;
        my $current_divider = 0;
        $current_divider = sprintf("%.3f", $current_divider);
        push(@divs, $current_divider);

        my $extra_feature_space = 0.1;
        my $space_left = 1 - $extra_feature_space;
        my $spacer = $space_left/$panel_count;
        for (my $i = 1; $i < $panel_count; $i++){
            $current_divider += $spacer;
            $current_divider = sprintf("%.3f", $current_divider);
            push(@divs, $current_divider);
        }
        my $div_string = join(",", @divs);
        $xml = "  <PanelLayout dividerFractions=\"$div_string\"/>";
    }
    return $xml;
}

sub generate_data_track_xml {
    my ($self) = @_;
    return join(q{},@{$self->_data_tracks});
}

sub xml {
    my ($self,$locus) = @_;
    my $genome_build = $self->genome;
    #Put all the XML blocks together
    my $starting_locus = "12:25398182-25398361";
    if($locus) {
        if($locus =~ /^.+:\d+-\d+$/) {
            $starting_locus = $locus;
        }
        else {
            die $self->error_message("Malformed locus specification $locus");
        }
    }
    my $xml_header = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
    my $session_begin = "<Session genome=\"$genome_build\" locus=\"$starting_locus\" version=\"$VERSION\">";
    my $session_end = "</Session>";

    #this is kind of dirty, but it lets us get a nicely formatted file and allows greater flexibility in how people use it.
    my @components = ($self->generate_resource_xml, $self->generate_data_track_xml, $self->generate_feature_panel_xml, $self->generate_panel_layout_xml);
    map { chomp } @components;
    return join("\n",$xml_header,$session_begin,@components,$session_end) . "\n";
}

1;

=pod

=head1 Name

Genome::File::IgvXML

=head1 Synopsis

Genome::File::IgvXML contains code for generating IGV XML sessions

=head1 Usage

my $igv = Genome::File::IgvXml->create(id => $filename);

$igv->add_bam_track(file => $bam, name => $sample)

my $igv_fh = $igv->open("w");

print $igv_fh $igv->xml($locus);

$igv_fh->close;

=head1 Methods

=head2 add_bam_track

=over

=item I<Synopsis>

adds a single bam file to the XML file

=item I<Arguments>

a hash of parameters. The following are supported (but only file is required):

=over

=item * file

the file path

=item * name

the label for the bam's track

=item * font_size

size of the font for this track's labels in points

=item * color_option

color option for reads e.g. 'UNEXPECTED_PAIR'

=item * max_depth

maximum depth of the coverage track associated with this BAM. Defaults to 100.

=item * min_depth

minimum depth of the coverage track associated with this BAM. Defaults to 0.

=back

=back

=head2 add_bed_track 

=over

=item I<Synopsis>

adds a BED file to the session

=item I<Arguments>

a hash of parameters. The following are supported (but only file is required):

=over

=item * file

the file path

=item * name

the label for the BED's track

=item * font_size

size of the font for this track's labels in points

=back

=back

=head2 add_junction_track 

=over

=item I<Synopsis>

Add a file of splice junctions

=item I<Arguments>

a hash of parameters. The following are supported (but only file is required):

=over

=item * file

the file path

=item * name

the label for the junction file's track

=item * font_size

size of the font for this track's labels in points

=back

=back

=head2 xml

=over

=item I<Synopsis>

Generate the XML representation of the file as a string

=item I<Arguments>

the starting locus position as a string. By default it is C<12:25398182-2539361>

=item I<Returns>

Text of the XML file

=back

=head1 See Also

L<Genome::File::Base>

=head1 Author(s)

B<David Larson> I<dlarson@genome.wustl.edu>, B<Malachi Griffith> I<mgriffit@genome.wustl.edu>

=cut
