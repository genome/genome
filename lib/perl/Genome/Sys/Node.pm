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
                              valid_values => [ __PACKAGE__->_supported_protocols() ],
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

sub _supported_protocols {
    return ('nfs','ssh','ftp','http');
}

sub attach {
    my $self = shift;
    my $protocol_override = shift;

    my @protocols_to_try;
    if ($protocol_override) {
        @protocols_to_try = ($protocol_override);
    }
    else {
        @protocols_to_try = $self->_supported_protocols;
    }
    $self->debug_message("protocols to test @protocols_to_try");

    my $is_already_attached_via = $self->attached_via;

    for my $protocol (@protocols_to_try) {
        my $method = "_attach_$protocol";
        unless ($self->can($method)) {
            $self->debug_message("no support for $protocol yet...\n");
            next;
        }

        $self->$method();
        
        my $underlying_mount_point = $self->_underlying_mount_point_for_protocol($protocol);
        my $mount_point_symlink = $self->mount_point;

        if (-e $mount_point_symlink) {
            unlink $mount_point_symlink;
            Genome::Sys->create_symlink($underlying_mount_point, $mount_point_symlink);  
        }

        if ($self->is_attached) {
            $self->status_message("attached " . $self->id . " via " . $protocol);
            last;
        }
        else {
            $self->error_message("error attaching " . $self->id . " via " . $protocol);
        }
    }

    return 1;
}

sub detach {
    my $self = shift;
    my @protocols = @_;
   
    unless (@protocols) {
        @protocols = $self->_supported_protocols();
    }

    my @errors;
    for my $protocol (@protocols) {
        my $underlying_mount_point = $self->_underlying_mount_point_for_protocol($protocol);
        unless (-e $underlying_mount_point) {
            next;
        }
        my $method = "_detach_$protocol";
        eval { $self->$method; };
        if ($@) {
            push @errors, $@;
        }
    }
    
    if (@errors) {
        die join("\n",@errors);
    }
    
    return 1;
}

sub attached_via {
    my $self = shift;
    my $mount_point_symlink = $self->mount_point;
    if (-l $mount_point_symlink) {
        # this can be true even if -e is false
        my $path = readlink($mount_point_symlink);
        my ($protocol) = ($path =~ /.([^\.]+$)/);
        return $protocol;
    }
    elsif (-e $mount_point_symlink) {
        die "expected $mount_point_symlink to be a symlink!";
    }
    elsif (not -e $mount_point_symlink) {
        return;
    }
}

sub _protocol_for_underlying_mount_point {
    my $self = shift;
    my $mount_point = shift;
    $mount_point ||= readlink($self->mount_point);
    unless ($mount_point) {
        die "no mount point given, and no mount point behind default symlink!";
    }
    my ($protocol) = ($mount_point =~ /.([^\.]+$)/);
    return $protocol;
}

sub _underlying_mount_point_for_protocol {
    my $self = shift;
    my $protocol = shift;
    die "no protocol specified!" unless $protocol;
    my $mount_point_symlink = $self->mount_point;
    my $underlying_mount_point = $mount_point_symlink;
    $underlying_mount_point =~ s|/opt/gms/|/opt/gms/.|;
    $underlying_mount_point .= '.' . $protocol;
    return $underlying_mount_point;
}

sub _attach_ftp {
    my $self = shift;
    my $hostname = $self->hostname;
    my $ftp_detail = $self->ftp_detail;
    my $underlying_mount_point = $self->_underlying_mount_point_for_protocol('ftp');
    my $cmd = "curlftpfs 'ftp://$hostname/$ftp_detail' '$underlying_mount_point' -o tcp_nodelay,kernel_cache,direct_io";
    Genome::Sys->shellcmd(cmd => $cmd);    

}

sub _detach_ftp {
    my $self = shift;
    my $underlying_mount_point = $self->_underlying_mount_point_for_protocol('ftp');
    my $cmd = "fusermount -u '$underlying_mount_point'";
    Genome::Sys->shellcmd(cmd => $cmd);
}

1;

