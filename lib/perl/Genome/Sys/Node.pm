use strict;
use warnings;
use Genome;

package Genome::Sys::Node;

class Genome::Sys::Node {
    id_by => [
        id => { is => 'Text', doc => 'the GMS system ID of the GMS in question' },
    ],
    has => [
        hostname => { is => 'Text' },
    ],
    has_optional => [
        id_rsa_pub    => { is => 'Text' },
        desc          => { is => 'Text' },
        ftp_detail    => { is => 'Number' },
        http_detail   => { is => 'Number' },
        ssh_detail    => { is => 'Number' },
        nfs_detail    => { is => 'Number' },
    ],
    has_calculated => [
        mount_point       => { is => 'FilesystemPath',
                              calculate_from => ['id'],
                              calculate => q|"/opt/gms/$id"|,
                              doc => 'the mount point for the system, when attached (/opt/gms/$ID)',
                            },

        is_current        => { is => 'Boolean',
                              calculate_from => ['id'],
                              calculate => q|$id eq $ENV{GENOME_SYS_ID}|,
                              doc => 'true for the current system',
                            },

        is_attached       => { is => 'Boolean', 
                              calculate_from => ['mount_point'],
                              calculate => q|-e $mount_point and not -e '$mount_point/NOT_MOUNTED'|,
                              doc => 'true when the given system is attached to the current system' 
                            },

        attachment_method => {
                              is => 'Text',
                              valid_values => ['nfs','ssh','ftp','http'],
                              doc => 'the protocol used to attach the system to the current one', 
                            },
    ],
    data_source => { 
        #uri => "file:$tmpdir/\$rank.dat[$name\t$serial]" }
        is => 'UR::DataSource::Filesystem',
        path  => $ENV{GENOME_HOME} . '/known-systems/$id.tsv',
        columns => ['hostname','id_rsa_pub','desc','ftp_detail','http_detail','ssh_detail','nfs_detail'],
        delimiter => "\t",
    },
};

sub attach {
    my $self = shift;
    my $protocol_override = shift;

    my @protocols_to_try;
    if ($protocol_override) {
        @protocols_to_try = ($protocol_override);
    }
    else {
        @protocols_to_try = @{ $self->__meta__->property('attachment_method')->valid_values };
    }
    $self->status_message("protocols to test @protocols_to_try");

    my $hostname = $self->hostname;

    for my $protocol (@protocols_to_try) {
        my $method = "_attach_$protocol";
        unless ($self->can($method)) {
            warn "no support for $protocol yet...\n";
            next;
        }
        $self->$method();
    }

    return 1;
}

sub mount_point_for_protocol {
    my $self = shift;
    my $mount_point_symlink = $self->mount_point;

}
sub _attach_ftp {
    my $self = shift;
    my $cmd = '';
    die "failed to ftp mount!";
}

1;

