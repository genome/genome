package Genome::Model::Tools::Wiki::DocumentCommands;

use strict;
use warnings;

class Genome::Model::Tools::Wiki::DocumentCommands {
    is => ['Genome::Model::Tools::Wiki'],
    has => [
        base_modules => { is => 'String', is_optional => 1, is_many => 1, shell_args_position => 1,
            default_value => ['Genome::Command','Workflow::Command','UR::Namespace::Command'],
            doc => 'Generate documentation for each subclass of these modules' },            
        svn_revision => { is => 'Number', doc => 'SVN revision to link to for command', is_optional => 1 },
        _wiki_basepath => { is => 'Text', is_constant => 1, default => 'Analysis Pipeline/Command Reference/' },
        _svn_basepath => { is => 'Text', is_constant => 1, default => 'http://svn/scm/viewvc/gscpan/perl_modules/trunk/' },
    ],
    doc => "Uploads documentation to the wiki for commands."
};

sub help_brief {
    __PACKAGE__->__meta__->doc
}

sub help_synopsis{
    return <<"EOS"
gmt wiki document-commands --svn-revision 50000
gmt wiki document-commands --svn-revision 50000 Genome::Model::Tools::Wiki
EOS
}

sub help_detail {
    return join("\n", 
        'This tool generates documentation for each command and uploads it to the',
        'wiki.  Additionally a link to a revision in Subversion is added if',
        'the svn-revision option is specified. One or more commands can be provided',
        'on the command line.  If none are given, this tool will document',
        '`genome`, `gmt`, `ur`, and `workflow` by default.  This tool should be run',
        'from within a directory of the namespace to document',
        '(e.g. /gsc/scripts/opt/genome/current/pipeline/lib/perl/Genome/)');
}

sub execute {
    my $self = shift;
    
    my @base_modules = $self->base_modules;
        
    my @modules = map( $self->get_all_subclasses($_), @base_modules);
    push @modules, @base_modules; #Do base_modules last in case child is mistaken about name
    
    my $edit_summary = 'gmt wiki document-commands: command auto-documentation upload';
    if(Genome::Sys->username) {
        $edit_summary .= ' by ' . Genome::Sys->username;
    }
    
    if($self->svn_revision) {
        $edit_summary .= ' from deploy of revision ' . $self->svn_revision;
    }

    for my $module (@modules) {
        my $wikified_pod;
        eval {
            $wikified_pod = $self->get_wikified_pod($module);

            if($module->is_sub_command_delegator) {
                $wikified_pod .= $self->get_sub_command_wikilinks_section($module);
            }

            $wikified_pod .= $self->get_module_link_section($module);
        };

        if($@) {
            $self->warning_message("Error generating wikified POD for $module:\n" . $@);
            next;
        }

        my $wiki_page_name = $self->get_wikilink_for_module($module);

        App::DBI->no_commit(0); #Some includes of some tools set this value, which prevents posting
        my $result = $self->post_page_to_wiki(
            page         => $wiki_page_name,
            page_text    => $wikified_pod,
            edit_summary => $edit_summary,
        );
        App::DBI->no_commit(1); #Just to be paranoid, turn this back on

        unless ($result) {
            $self->error_message("Couldn't post to wiki [[$wiki_page_name]] for $module. Upload aborted.");
            return;
        }
    }
}

sub get_all_subclasses {
    my $self = shift;
    my $module = shift;

    my @submodules;
    eval {
        @submodules = $module->sub_command_classes;
    };
    
    if($@) {
        $self->warning_message("Error getting subclasses for module $module: " . $@);
    }

    return unless @submodules and $submodules[0]; #Sometimes sub_command_classes returns 0 instead of the empty list

    #Preprocess children (mitigates modules overwriting "gmt")
    return map($self->get_all_subclasses($_), @submodules), @submodules;
}

#TODO There are some minor formatting issues with the generated wiki
#either cleanup or eventually move away from using pod2wiki
sub get_wikified_pod {
    my $self = shift;
    my $module = shift;
    
    my $pod = $module->help_usage_command_pod;
        
    my $parser = Pod::Simple::Wiki->new('mediawiki');
    $parser->no_whining(1);
    $parser->no_errata_section(1);
    $parser->complain_stderr(0);

    my $output_wiki_text;
    $parser->output_string( \$output_wiki_text );
    $parser->parse_string_document($pod);
    
    return $output_wiki_text;
}

sub get_wikilink_for_module {
    my $self = shift;
    my $module = shift;

    my $command_name = $module->command_name;
    
    if($command_name eq 'gmt' and $module ne 'Genome::Model::Tools') {
        warn("Module $module thinks its command name is 'gmt'. Check inheritance.");
    }
    
    $command_name =~ s! !/!g; #Use slashes to allow automatic wiki breadcrumbs to work

    my $wiki_page_name = $self->_wiki_basepath . $command_name;

    return $wiki_page_name;
}

sub get_sub_command_wikilinks_section {
    my $self = shift;
    my $module = shift;

    my @sub_command_classes = $module->sub_command_classes;
    return '' unless @sub_command_classes;

    my $links = "\n\n=Subcommand Links=\n";
    for my $sub_module (@sub_command_classes) {
        $links .= "*[[" . $self->get_wikilink_for_module($sub_module) . '|' . $sub_module->command_name . "]]\n";
    }

    return $links;
}

sub get_module_link_section {
    my $self = shift;
    my $module = shift;
    
    if($module =~ m!^UR!) {
        #Possibly get appropriate github link?
        return ''; #UR is not kept in Subversion
    }

    my $body = '';
    
    if($module =~ m!^Genome!) {
        $body .= '* ' . Genome::Model::Tools::Wiki::DocumentClasses->create()->class_to_wikilink($module);
    }
    
    if($self->svn_revision) {
        my $module_path = $module;
        $module_path =~ s!::!/!g;
        $module_path .= '.pm';
        my $url = $self->_svn_basepath . $module_path . '?view=markup&pathrev=' . $self->svn_revision;
    

        $body .= "\n* [$url $module in SVN (rev " . $self->svn_revision . ")]";
    }
    
    return $self->create_section('=See Also=', $body);
}
