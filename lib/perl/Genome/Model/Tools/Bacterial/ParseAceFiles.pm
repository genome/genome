package Genome::Model::Tools::Bacterial::ParseAceFiles;

use strict;
use warnings;
use Genome;
use IPC::Run;
use IO::Dir;
use Carp;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        acedb_path => { is => 'String',
                        doc => "path to acedb location",
                      },
        acefiles_path => { is => 'String',
                           doc => "path to acefiles",
                         },
    ],
    has_optional => [
        tace_path => { is => 'String',
                       doc => "path to tace executable",
                       default => '/gsc/scripts/bin/tace',
                     },

        dry_run => { is => 'Boolean',
                     doc => "don't actually parse into acedb",
                     default => 0,
                 },
    ],
);

sub help_brief
{
    "tool for parsing a directory of acefiles";

}

sub help_synopsis
{
    my $self = shift;
    return <<"EOS"
parse a directory of ace files
EOS
}

sub help_detail
{
    my $self = shift;
    return <<"EOS"
parse a directory of ace files
EOS

}


sub execute
{
    my $self = shift;
    # find ace files and fasta files to parse in
    my @acefiles;
    my $d = IO::Dir->new($self->acefiles_path);
    while(my $entry = $d->read) {
        next if $entry =~ /^\.\.?$/;
        next if $entry =~ /\.gff$/;
        next if $entry =~ /\.txt$/;
        push(@acefiles, $self->acedb_path."/".$entry);
    }
    $d->close;
    $self->status_message("found ".scalar(@acefiles)." acefiles to parse");
    #interactively  start up tace

    my @command = ($self->tace_path, $self->acedb_path);
    
    my $input_stream;

    # set up stream of commands to send to tace. 
    $input_stream = join("\n", map { "parse ". $_ } @acefiles) ."\nsave\nquit\n\n";

    my ($stdout,$stderr);
    if($self->dry_run) {
        $self->status_message("dry-run turned on - not running commands");
        print $input_stream,"\n";
        return 1;
    }
    my $rv = IPC::Run::run( \@command,
        '<',
        \$input_stream,
        '>',
        \$stdout,
        '2>',
        \$stderr, );

    # check return value. 
    unless($rv) {
        $self->error_message("problem parsing ace files, error: $stderr");
        croak;
    }




    return 1;
}


1;
