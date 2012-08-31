package Genome::Model::Tools::Wiki::DocumentClasses;

use strict;
use warnings;

class Genome::Model::Tools::Wiki::DocumentClasses {
    is => ['UR::Namespace::Command::RunsOnModulesInTree', 'Genome::Model::Tools::Wiki'],
    has => [
        svn_revision => { is => 'Number', doc => 'SVN revision to link to for module', is_optional => 1 },
        _wiki_basepath => { is => 'Text', is_constant => 1, default => 'Analysis Pipeline Group/Class reference/' },
        _svn_basepath => { is => 'Text', is_constant => 1, default => 'http://svn/scm/viewvc/gscpan/perl_modules/trunk/' },
    ],
    doc => "Uploads documentation to the wiki for classes in the Genome namespace.",
};

sub help_brief {
    __PACKAGE__->__meta__->doc
}

sub help_synopsis{
    return <<"EOS"
gmt wiki document-classes --svn-revision 50000
gmt wiki document-classes --svn-revision 50000 Genome::Model::Tools::Wiki
EOS
}

sub help_detail {
    return join("\n", 
        'This tool generates documentation for each class and uploads it to the',
        'wiki. Additionally a link to a revision in Subversion is added if',
        'the svn-revision option is specified.  If classes are specified on the',
        'command line they will be used.  Otherwise the entire namespace will be',
        'documented.  This tool must be run from within a directory of the ',
        'namespace to document. (e.g. /gsc/scripts/opt/genome/current/pipeline/lib/perl/Genome/)');
}

#See also UR::Object::ModuleWriter->resolve_class_description_perl
sub for_each_class_object {
    my $self = shift;
    my $class_to_doc = shift;
    
    my $doc = '';
    eval {        
        #Module Name
        my $heading_name = ($class_to_doc->is_abstract ? 'Abstract Class' : 'Class');
        $doc .= $self->create_section($heading_name, $class_to_doc->class_name);
        
        #Inheritance
        my $parents = $class_to_doc->is;
        if($parents and scalar @$parents) {
            $doc .= $self->create_section('==Inherits From==', join(', ', map($self->class_to_wikilink($_), @$parents) ));
        }
        
        #Database Information
        if($class_to_doc->table_name) {
            $doc .= $self->create_section('==Database Table==', $class_to_doc->table_name);
        }
        
        #Class Documentation
        if($class_to_doc->doc) {
            $doc .= $self->create_section('Description', $class_to_doc->doc);
        }
        
        #Properties
        $doc .= $self->create_section('Properties', '');
        
        my %id_property_names = map { $_->property_name => 1 } $class_to_doc->direct_id_token_metas;
        
        my @direct_properties = $class_to_doc->direct_property_metas;
        
        my %properties_by_type;
        for my $property (@direct_properties) {
            my $location = $property->is_specified_in_module_header; 
            next unless $location;
            
            my ($header_type) = $location =~ m/::(\w+)$/;
            
            if($id_property_names{$property->property_name}) {
                $header_type = 'id_by';
            }
            
            push @{ $properties_by_type{$header_type} }, $property;
        }
        
        my @order = ('id_by', 'has', 'has_many', 'has_optional', 'has_many_optional');
        
        for my $type (@order) {
            next unless ($properties_by_type{$type});
            
            my $section_text = '';
            my @properties = sort @{ $properties_by_type{$type} };
            for my $property (@properties) {
                $section_text .= $self->property_to_wikitext($property);
            }
            
            $doc .= $self->create_section('=' . $type . '=', $section_text);
        }
        
        if($self->svn_revision) {
            my $module = $class_to_doc->class_name;
            my $module_path = $module;
            $module_path =~ s!::!/!g;
            $module_path .= '.pm';
            my $url = $self->_svn_basepath . $module_path . '?view=markup&pathrev=' . $self->svn_revision;
        
            my $svn_link = "\n* [$url $module in SVN (rev " . $self->svn_revision . ")]";
            
            $doc .= $self->create_section('See Also', $svn_link);
        }        
    };
    
    if($@) {
        $self->warning_message('Error generating documentation for ' . $class_to_doc->class_name . ': ' . $@);
        return -1; #Non-fatal error--just skip to next module
    }
    
    my $edit_summary = 'gmt wiki document-classes: class auto-documentation upload';
    if(Genome::Sys->username) {
        $edit_summary .= ' by ' . Genome::Sys->username;
    }
    
    if($self->svn_revision) {
        $edit_summary .= ' from deploy of revision ' . $self->svn_revision;
    }
    
    my $class_name = $class_to_doc->class_name;
    my $wiki_page_name = $self->class_to_wikipath($class_name);
    
    App::DBI->no_commit(0); #Some classes (typically test modules) turn on no-commit, which prevents posting
    my $result = $self->post_page_to_wiki(
        page         => $wiki_page_name,
        page_text    => $doc,
        edit_summary => $edit_summary,
    );
    App::DBI->no_commit(1); #Just to be paranoid, turn this back on

    unless ($result) {
        $self->error_message("Couldn't post to wiki [[$wiki_page_name]] for $class_name. Upload aborted.");
        return;
    }
    
    return 1;
}

sub class_to_wikilink {
    my $self = shift;
    my $class_name = shift;
    my $section = shift;
    
    my $display_part = $section || $class_name;
    
    return $display_part unless($class_name =~ /^Genome/); #Only link to other classes within our namespace
    
    my $section_path = $section? '#' . $section : '';
    
    return "[[" . $self->class_to_wikipath($class_name) . $section_path . '|' . $display_part . ']]';
}

sub class_to_wikipath {
    my $self = shift;
    my $class_name = shift;
    
    my $class_path = $class_name;
    $class_path =~ s!::!/!g;
    
    return $self->_wiki_basepath . $class_path;
}

sub property_to_wikitext {
    my $self = shift;
    my $property = shift;
    
    my $body = '';
    
    my $data_type_text = $property->data_type || 'undefined';
    $body .= $self->create_item('type', $self->class_to_wikilink($data_type_text));
    
    if($property->is_calculated) {
        my $calculate_from = $property->calculate_from;
        if($calculate_from) {
            #TODO This should link to the relevant superclass if calculating from an inherited property
            my @values = map($self->class_to_wikilink($property->class_meta->class_name, $_), @$calculate_from);
            $body .= $self->create_item('calculate_from', @values);
        }

        my $calculate_found = 0;
        foreach my $calculate_type ( qw( calculate calculate_sql calculate_perl calculate_js ) ) {
            if ($property->$calculate_type) {
                $calculate_found = 1;
                
                #Prepend a space to each line (Mediawiki's pre-format notation)
                my $calculate_display = "\n " . join("\n ", split("\n", $property->$calculate_type));
                
                $body .= $self->create_item($calculate_type, $property->$calculate_type);
            }
        }

        $body .= $self->create_item('is_calculated', '1') unless ($calculate_found);
    } 
    elsif ($property->column_name) {   
        if (uc($property->column_name) ne uc($property->property_name)) {
            $body .= $self->create_item('column_name', $property->column_name);
        }
    }

    if (defined($property->default_value)) {
        my $value = $property->default_value;
        $body .= $self->create_item('default_value', '<nowiki>' . $property->default_value . '</nowiki>');
    }
    
    if (defined($property->implied_by)) { 
        $body .= $self->create_item('implied_by', $self->class_to_wikilink( $property->class_meta->class_name, $property->implied_by ));
    }

    if (my @id_by = $property->id_by_property_links) {
        $body .= $self->create_item('id_by', map( $self->class_to_wikilink( $data_type_text, $_->property_name), @id_by));
    }

    if ($property->via) {
        $body .= $self->create_item('via', $self->class_to_wikilink( $property->class_meta->class_name, $property->via ));

        if ($property->to and $property->to ne $property->property_name) {
            #TODO This won't be correct when a subclass of the via property's data_type is implicitly expected
            my $via_property = $property->via_property_meta;
            my $to_class = $via_property? $property->via_property_meta->data_type : '';
            $body .= $self->create_item('to', $self->class_to_wikilink($to_class, $property->to) );
        }
    }
    if ($property->reverse_as) {
        $body .= $self->create_item('reverse_as', $self->class_to_wikilink($property->data_type, $property->reverse_as) );
    }

    if ($property->valid_values) {
        my @values = map('<nowiki>' . $_ . '</nowiki>', @{$property->valid_values});
        $body .= $self->create_item('valid_values', @values);
    }
    
    return $self->create_section('==' . $property->property_name . '==', $body);
}

1;
