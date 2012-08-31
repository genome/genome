
package Genome::Model::Tools::Sam::BamBuster;

use strict;
use warnings FATAL => 'all';

use Genome;
use File::Path;

class Genome::Model::Tools::Sam::BamBuster {
    is  => 'Genome::Model::Tools::Sam',
    has_input => [
        input => {
            is  => 'String',
            doc => 'Input SAM/BAM file to extract reads from. Required.',
        },
        output_directory => {
            is          => 'String',
            doc         => 'Output directory where data will be dumped.  In here will be subdirectories for each library',
        },
    ],
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
    unless ($fd = IO::File->new(sprintf("%s view -h %s", $self->samtools_path, $self->input) . "|")) {
        $self->error_message("can't open the file " . $self->input);
        return;
    }

    my $SAM_HEADER;

    my $line;

    my %rg_handles;

    while ($line = <$fd>) {
        if (substr($line,0,1) eq '@') {
            $SAM_HEADER .= $line;
            if (substr($line,0,3) eq '@RG') {
                my $rg_info;
                unless ($rg_info = $self->get_read_group_from_sam_header($line)) {
                    $self->error_message("failed parsing read group... " . $self->error_message);
                    return;
                }
                #print Data::Dumper::Dumper($rg_info);
                $self->status_message(sprintf("Opening bam file for writing for RG ID %s library %s platform unit %s", $rg_info->{id}, $rg_info->{library_name}, $rg_info->{platform_unit}));
                print $self->output_directory, "\n";
                my $path = sprintf("%s/%s", $self->output_directory, $rg_info->{library_name});
                unless (-d $path) {
                   mkpath($path); 
                }
                unless (-d $path && -w $path) {
                    $self->error_message("$path is not a valid directory or is not writable by you, " . Genome::Sys->username . ". Please check permissions.");
                    return;
                }
                my $handle = IO::File->new(sprintf("|%s view -S -b -o %s/%s.%s.bam -", $self->samtools_path, $path, $rg_info->{id}, $rg_info->{platform_unit}));
                $rg_handles{$rg_info->{id}} = $handle;
            }
        } else {
            last;
        }
    }

    # $line has the first line of real sam in it so we have to make sure to take care of it too before looping over the rest

    for my $handle (values %rg_handles) {
        print $handle $SAM_HEADER;
    }

    unless ($self->process_line($line, \%rg_handles)) {
        $self->error_message("failed processing $line: " . $self->error_message);
        return;
    }

    while ($line = <$fd>) {
        unless ($self->process_line($line, \%rg_handles)) {
            $self->error_message("failed processing $line: " . $self->error_message);
            return;
        }
    }
    
    $fd->close;

    for my $handle (values %rg_handles) {
        $handle->close;
    }

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
    
    my ($id) = $line =~ m/ID:(.*?)\t/;
    unless (defined $id)  {
        $self->error_message("failed to parse read group id from SAM line: $line");
        return;
    }
    
    my ($platform_unit) = $line =~ m/PU:(.*?)\t/;
    unless (defined $platform_unit)  {
        $self->error_message("failed to parse platform unit from SAM line: $line");
        return;
    }
    
    my ($library_name) = $line =~ m/LB:(.*?)\t/;
    unless (defined $library_name)  {
        $self->error_message("failed to parse library name from SAM line: $line");
        return;
    }

    return {id=>$id, platform_unit=>$platform_unit, library_name=>$library_name};
}


1;
__END__

