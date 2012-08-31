
package Genome::Model::Tools::Galaxy::Update;

use strict;
use warnings;
use Genome;
use File::Copy;

class Genome::Model::Tools::Galaxy::Update {
    is  => 'Command::V2',
    has => [
        path => {
            is  => 'Text',
            is_optional => 1,
            shell_args_position => 1,
            doc => 'Galaxy installation directory to use.  Defaults to $HOME/galaxy.'
        },
        pull => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Update the Galaxy software itself, not just gmt tool configuration.',
            default => 0
        }
    ],
    doc => 'update the Galaxy configuration for all of GMT',
};

sub sub_command_sort_position { 2 } 

sub execute {
    my $self = shift;

    my $path = $self->path;
    if (!defined($path)) {
        $path = $ENV{HOME} . "/galaxy/";
    }

    # look for key files to make sure path is galaxy directory
    my @key_files = (".hg", "run_galaxy_listener.sh", "run.sh");
    foreach my $k (@key_files) {
        my $file_path = $path . "/" . $k;
        unless (-e $file_path) {
            $self->error_message("Does not appear to be valid galaxy folder");
            die();
        }
    }

    if ($self->pull) {
        $self->status_message("updating the Galaxy software itself...");
        chdir($path);
        Genome::Sys->shellcmd(cmd => 'hg pull -u');
    }

    # TODO: hit the whole GMT tree
    my @xml_files = ();
    for my $ns ('Genome::Model::Tools::Music') {
        eval "use $ns";
        if ($@) {
            die $self->error_message("Failed to use $ns: $@");
        }
        my $c = $ns;
        $c =~ s/::/\//g;
        $c .= '.pm';
        my $path = $INC{$c};
        $path =~ s/.pm//;        
        $self->status_message("Checking $path for galaxy configuration...");
        push @xml_files, `find $path | grep .galaxy.xml`;
        chomp @xml_files;
        $self->status_message("Found " . scalar(@xml_files) . " tools with galaxy configuration.");
        $self->status_message("Found galaxy configuration for: @xml_files\n");
    }

    Genome::Sys->create_directory("$path/tools/gmt");
    
    my @tool_names;
    foreach my $xml_file (@xml_files) {
        my $fn = $xml_file;
        $fn =~ s|^.*Genome/Model/Tools/||;
        $fn =~ s|/|-|g;
        $fn = lc($fn);
        $fn =~ s/.pm.galaxy.xml$/.xml/;
        push @tool_names,$fn;
        $self->status_message("copying to $path/tools/gmt/$fn from $xml_file");
        copy($xml_file, "$path/tools/gmt/$fn");
    }

    # handle tool_conf.xml rewrite
    # read tool_conf.xml
    open(tool_conf_ifh, '<', "$path/tool_conf.xml");
    print "$path/tool_conf.xml" . "\n";
    my $tool_xml = '';
    while (<tool_conf_ifh>) {
        $tool_xml .= $_;
    }
    close(tool_conf_ifh);
    
    # Either replace existing Genome section or put it right after toolbox tag
    my $new_genome_section = '<section name="GMT" id="gmt">' . "\n";
    foreach my $tool (@tool_names) {
        $new_genome_section .= '    <tool file="gmt/' . $tool . '" />' . "\n"; 
    }
    $new_genome_section .= "</section>\n";
    unless ($tool_xml =~ s/^\s+<section name="GMT" id="gmt">.*?<\/section>\n/$new_genome_section/ms) {
        $tool_xml =~ s/^<toolbox>\n/<toolbox>\n$new_genome_section/ms;
    }
     
    # write tool_conf.xml
    $self->status_message("replacing $path/tool_conf.xml...");
    copy("$path/tool_conf.xml","$path/tool_conf.xml.bak");
    open(tool_conf_ofh, '>', "$path/tool_conf.xml");
    print tool_conf_ofh $tool_xml;
    close(tool_conf_ofh);

    $self->status_message("galaxy can now be restarted!\n");
    return 1;
}


