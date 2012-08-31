
package Genome::Model::Tools::Sam::ListReadGroups;

use strict;
use warnings FATAL => 'all';

use Genome;
use File::Path;

class Genome::Model::Tools::Sam::ListReadGroups {
    is  => 'Genome::Model::Tools::Sam',
    has_input => [
        input => {
            is  => 'String',
            doc => 'Input SAM/BAM file to list read groups from.',
        },
        silence_output=> {
            is => 'Boolean',
            doc => 'Don\'t print the output to stdout',
            default_value => 0,
        }
    ],
    has_output => [
        read_groups => {
            is => 'Array',
            doc => 'List of read groups for use programmatically',
            is_optional=>1,
            is_many=>1
        }
    ]
};

sub help_brief {
    'Tool to break apart BAMs into smaller BAMs based on their component read groups'
}

sub help_detail {
    return <<EOS
Tool to break apart BAMs into smaller BAMs based on their component read groups
EOS
}

sub execute {
    my $self = shift;

    my $fd;
    unless ($fd = IO::File->new(sprintf("%s view -H %s %s", $self->samtools_path, ($self->input =~ m/\.sam$/ ? "-S" : ""), $self->input) . "|")) {
        $self->error_message("can't open the file " . $self->input);
        return;
    }

    my @read_groups;

    while (my $line = <$fd>) {
        if (substr($line,0,3) eq '@RG') {
            my $rg_info;
            unless ($rg_info = $self->get_read_group_from_sam_header($line)) {
                $self->error_message("failed parsing read group... " . $self->error_message);
                return;
            }
            push @read_groups, $rg_info;
        }
    }

    my @rg_ids = map {$_->{id}} @read_groups;

    unless ($self->silence_output) {
        print "Read Groups detected from BAM file: \n" . join "\n", @rg_ids;
    }
    $self->read_groups(\@rg_ids);

    return 1;
}

sub process_line {
my $self = shift;

my $line = shift;
my $rg_handles = shift;

my ($rg_id) =  $line =~ m/RG:Z:(.*?)(\t|\s+)/;
unless (defined $rg_id) {
    $self->error_message("can't get read group ID");
    return;
}

my $handle = $rg_handles->{$rg_id};

print $handle $line;
}


sub get_read_group_from_sam_header {
my $self = shift;
my $line = shift;

my ($id, $platform_unit, $library_name);

($id) = $line =~ m/ID:(.*?)(\t|\s+)/;
unless (defined $id)  {
    $self->error_message("failed to parse read group id from SAM line: $line");
    return;
}

($platform_unit) = $line =~ m/PU:(.*?)(\t|\s+)/;
if ($line =~ m/PU:/ && !defined $platform_unit)  {
    $self->error_message("failed to parse platform unit from SAM line: $line");
    return;
}

($library_name) = $line =~ m/LB:(.*?)(\t|\s+)/;
if ($line =~ m/LB:/ && !defined $library_name) {
    $self->error_message("failed to parse library name from SAM line: $line");
    return;
}

return {id=>$id, platform_unit=>$platform_unit, library_name=>$library_name};
}


1;
__END__

