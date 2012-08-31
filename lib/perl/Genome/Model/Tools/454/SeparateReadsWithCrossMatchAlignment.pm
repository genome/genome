package Genome::Model::Tools::454::SeparateReadsWithCrossMatchAlignment;

use strict;
use warnings;

use Genome;

use File::Basename;

class Genome::Model::Tools::454::SeparateReadsWithCrossMatchAlignment {
    is => ['Genome::Model::Tools::454'],
    has => [
            sff_file => {
                            is_input => 1,
                            is => 'String',
                            doc => 'The sff file to divide into smaller sff files',
                     },
            cross_match_file => {
                                 is_input => 1,
                                 is => 'String',
                                 doc => 'The output file path',
                             },
        ],
};

sub help_brief {
    "separate reads from an sff file based on the results of a cross_match alignment"
}

sub help_detail {
    return <<EOS
parses the results of a cross_match alignment and divides the reads(query) based on the primer(hit) in the alignment
the cross_match alignment is assumed to have been ran with the following parameters:
'-tags -masklevel 0 -gap1_only -minmatch 6 -minscore 6 -minmargin 1'
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('This genome-model tool '. $self->command_name .' will only run on 64-bit');
        return;
    }
    unless (-s $self->cross_match_file) {
        $self->error_message('cross_match output file '. $self->cross_match_file .' does not exist or is zero size.');
        return;
    }
    unless (-s $self->sff_file) {
        $self->error_message('sff file '. $self->sff_file .' does not exist or is zero size.');
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $cross_match_file_basename = basename($self->cross_match_file);

    my $out_sff_file_root_name = $self->sff_file;
    $out_sff_file_root_name =~ s/\.sff$//;

    my $reader = Genome::Model::Tools::Crossmatch::Reader->new(
                                                    input => $self->cross_match_file,
                                                );
    unless ($reader) {
        $self->error_message('Failed to create cross_match reader');
        return;
    }

    my %query_hit;
    my %hit_query_ref;
    my %hit_file;
    my %open_fhs;
    while(my $hit = $reader->next) {
        my $query_name = $hit->{query_name};
        my $subject_name = $hit->{subject_name};
        if (defined $query_hit{$query_name}) {
            $self->error_message('Expecting one hit per query'. $query_name .' and found more in '. $self->cross_match_file);
            return;
        }
        $query_hit{$query_name} = $subject_name;
        push @{$hit_query_ref{$subject_name}}, $query_name;
        unless (defined $hit_file{$subject_name}) {
            my $out_file = $self->_tmp_dir .'/'. $cross_match_file_basename .'.'. $subject_name .'.reads';
            my $fh = IO::File->new($out_file,'w');
            $hit_file{$subject_name} = $out_file;
            $open_fhs{$out_file} = $fh
        }
        my $hit_fh = $open_fhs{$hit_file{$subject_name}};
        unless ($hit_fh) {
            $self->error_message('For some reason no filehanlde exists for hit '. $subject_name);
            return;
        }
        print $hit_fh $query_name ."\n";
    }

    for my $open_fh (values %open_fhs) {
        $open_fh->close;
    }

    for my $hit (keys %hit_query_ref) {
        $self->status_message($hit ."\t". scalar(@{$hit_query_ref{$hit}}));

        my $reads_file = $hit_file{$hit};
        my $out_sff_file = $out_sff_file_root_name .'.'. $hit .'.sff';
        $self->status_message('Writing '. scalar(@{$hit_query_ref{$hit}}) .' reads to file '. $out_sff_file);
        my $separate_reads_sfffile = Genome::Model::Tools::454::Sfffile->create(
                                                                                in_sff_files => [$self->sff_file],
                                                                                out_sff_file => $out_sff_file,
                                                                                params => '-i '. $reads_file,
										version => $self->version,
										version_subdirectory => $self->version_subdirectory,
                                                                            );
        unless ($separate_reads_sfffile) {
            $self->error_message('Failed to create sffile genome-model tool');
            return;
        }
        unless ($separate_reads_sfffile->execute) {
            $self->error_message('Failed to execute command '. $separate_reads_sfffile->command_name);
            return;
        }
    }

    return 1;
}

1;


