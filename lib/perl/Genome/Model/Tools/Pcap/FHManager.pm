package Genome::Model::Tools::Pcap::FHManager;

use Carp;
my %file_handles;
sub new { return __PACKAGE__ }

sub get_fh
{
    my ($pointless,$file_name) = @_;
    Carp::confess "Filename is not defined" unless $file_name;
    Carp::confess "$file_name does not exist" unless (-e $file_name);
    if(exists $file_handles{$file_name}&& defined $file_handles{$file_name})
    {
        return $file_handles{$file_name};
    }
    else
    {
        my $fh = IO::File->new($file_name);
        Carp::confess "There was an error opening $file_name" unless $fh;
        #print "Opening $file_name...\n";
        $file_handles{$file_name} = $fh;
        return $file_handles{$file_name};
    }

}

1;
