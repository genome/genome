package Genome::Model::Tools::DetectVariants::Dispatcher;

use strict;
use warnings;

use Genome;
use Parse::RecDescent;

class Genome::Model::Tools::DetectVariants::Dispatcher {
    is => ['Genome::Model::Tools::DetectVariants::Somatic'],
    has_optional => [
        snv_detector_name => {
            is => "String",
            doc => 'The variant detector(s) to use for finding SNPs',
        },
        snv_detector_version => {
            is => "Version",
            doc => "Version"
        },
        indel_detector_name => {
            is => "String",
            doc => 'The variant detector(s) to use for finding indels',
        },
        indel_detector_version => {
            is => "Version",
            doc => "Version"
        },
        sv_detector_name => {
            is => "String",
            doc => 'The variant detector(s) to use for finding svs',
        },
        sv_detector_version => {
            is => "Version",
            doc => "Version"
        },
        control_aligned_reads_input => {
            doc => 'Location of the control aligned reads file to which the input aligned reads file should be compared (if using a detector that needs one)'
        },
    ],
    has_constant_optional => [
        version => {}, #We need separate versions for the dispatcher
    ],
    has_constant => [
        variant_types => {
            is => 'ARRAY',
            value => [('snv', 'indel', 'sv')],
        },
        #These can't be turned off--just pass no detector name to skip
        detect_snvs => { value => 1 },
        detect_indels => { value => 1 },
        detect_svs => { value => 1 },
    ],
    doc => 'This tool is used to handle delegating variant detection to one or more specified tools and combining the results',
};

sub help_brief {
    "A dispatcher for variant detectors.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants dispatcher ...
EOS
} #TODO Fill in this synopsis with a few examples

sub help_detail {
    return <<EOS 
A variant detector(s) specified under snv-detector, indel-detector, or sv-detector must have a corresponding module under `gmt detect-variants`.

In a future version, the parameter may be a boolean expression such as "(samtools && var-scan) || maq". Not in this one, however.
EOS
}

sub _should_skip_execution {
    my $self = shift;
    
    for my $variant_type (@{ $self->variant_types }) {
        my $name_property = $variant_type . '_detector_name';
        
        return if defined $self->$name_property;
    }
    
    $self->status_message('No variant detectors specified.');
    return $self->SUPER::_should_skip_execution;
}

sub _detect_variants {
    my $self = shift;

    my $detector_trees = {};
    my $detectors_to_run = {};
    
    
    #Build up detection request
    for my $variant_type (@{ $self->variant_types }) {
        my $name_property = $variant_type . '_detector_name';
        
        if($self->$name_property) {
            my $detector_tree = $self->parse_detector_string($self->$name_property);
            
            $detector_trees->{$variant_type} = $detector_tree;
            
            my $version_property = $variant_type . '_detector_version';
            my $param_property = $variant_type . '_params';
            
            my $versions = [split(':', ($self->$version_property || ''))];
            my $params = [split(':', ($self->$param_property || ''))];
            
            unless(scalar @$versions eq scalar @$params) {
                $self->error_message('Inconsistent number of versions and parameters passed!');
                die($self->error_message);
            }
            
            $self->build_detector_list($detector_tree, $detectors_to_run, $variant_type, $versions, $params);
        }
    }    
    
    #Run individual detectors
    for my $detector (keys %$detectors_to_run) {
        my $detector_module = $self->detector_class($detector);
                
        unless($self->is_valid_detector($detector_module)) {
            $self->error_message('Did not find expected module (' . $detector_module . ') for named detector ' . $detector);
            die($self->error_message);
        }   
        
        for my $version (keys %{ $detectors_to_run->{$detector} }) {
            my %run_parameters = ();
            
            for my $variant_type (@{ $self->variant_types }) {
                next unless($detector_module->can('detect_' . $variant_type . 's'));
                
                if(exists $detectors_to_run->{$detector}{$version}{$variant_type}) {
                    my $param_value = shift @{ $detectors_to_run->{$detector}{$version}{$variant_type} };
                    my $param_name = $variant_type . '_params';
                    
                    $run_parameters{$param_name} = $param_value;
                    $run_parameters{'detect_' . $variant_type . 's'} = 1;
                } else {
                    $run_parameters{'detect_' . $variant_type . 's'} = 0;
                }
            }
            
            #Make sure at least one thing to detect is turned on (This also pokes at params, but we don't set any unless something is 1)
            unless(grep($_ eq 1, values %run_parameters)) {
                last RUN_DETECTORS;
            }
            
            #TODO What to do about MAQ?  Check for a .map file and either complain or generate one (within the tool)?
            #Similarly for FASTA versus binary FASTA...  (Really our tools should accept a common format and MAQ should figure this out)
            $run_parameters{version} = $version;
            $run_parameters{reference_sequence_input} = $self->reference_sequence_input;
            $run_parameters{aligned_reads_input} = $self->aligned_reads_input;
            
            #TODO We should really use a temporary directory and copy the results to the working subdirectory after they're completed
            #TODO When we enable the boolean expressions, the directory name will need to take into account parameters as well
            my $command_output_directory = $self->calculate_detector_output_directory($detector, $version, '');
            Genome::Sys->create_directory($command_output_directory);
            $run_parameters{output_directory} = $command_output_directory;
            
            
            #TODO Make a workflow out of all these? (This requires the individual tools to be responsible for temporary file handling, etc.
            #Or perhaps we just make a dispatcher tool-runner wrapper that does it)
            my $detector_command = $detector_module->create(%run_parameters);
            unless($detector_command) {
                $self->error_message('Failed to instantiate command for ' . $detector. '!');
                die($self->error_message);
            }
            
            unless($detector_command->execute) {
                $self->error_message('Failed to execute ' . $detector . ': ' . $detector_command->error_message);
                die($self->error_message);
            }
        }
    }
    
    #Combine results
    for my $variant_type (@{ $self->variant_types }) {
        my $name_property = $variant_type . '_detector_name';
        
        if($self->$name_property) {
            my $detector_tree = $detector_trees->{$variant_type};
        
            my $version_property = $variant_type . '_detector_version';
            my $param_property = $variant_type . '_params';
            
            my $versions = [split(':', $self->$version_property)];
            my $params = [split(':', $self->$param_property)];
            
            my $final_files = $self->combine_results($detector_tree, $variant_type, $versions, $params);
            
            my $output_file_property = '_' . $variant_type . '_staging_output';
            my $filtered_output_file_property = '_filtered' . $output_file_property;
            
            Genome::Sys->copy_file($final_files->[0], $self->$output_file_property);
            Genome::Sys->copy_file($final_files->[1], $self->$filtered_output_file_property);
        }
    }
    
    return 1;
}

sub _verify_inputs {
    my $self = shift;
    
    #Skip the somatic checks since we might not be running a somatic detector.  (If we are the checks will be performed then.)
    return $self->Genome::Model::Tools::DetectVariants::_verify_inputs;
}

sub calculate_detector_output_directory {
    my $self = shift;
    my ($detector, $version, $param_list) = @_;
    
    my $subdirectory = join('-', $detector, $version, $param_list);
    
    return $self->output_directory . '/' . Genome::Utility::Text::sanitize_string_for_filesystem($subdirectory);
}

sub parse_detector_string {
    my $self = shift;
    my $string = shift;

    return unless $string;
    
    my $parser = $self->parser;
    
    my $result = $parser->startrule($string);
    
    unless($result) {
        $self->error_message('Failed to interpret detector string from ' . $string);
        die($self->error_message);
    }
    
    return $result;
}

sub is_valid_detector {
    my $self = shift;
    my $detector_class = shift;
    
    return if $detector_class eq $self->class; #Don't nest the dispatcher!
    
    my $detector_class_base = 'Genome::Model::Tools::DetectVariants';
    return $detector_class->isa($detector_class_base);
}

sub detector_class {
    my $self = shift;
    my $detector = shift;
    
    $detector = join ("", map { ucfirst(lc($_)) } split(/-/,$detector) ); #Convert things like "foo-bar" to "FooBar"
    
    my $detector_class_base = 'Genome::Model::Tools::DetectVariants';
    my $detector_class = join('::', ($detector_class_base, $detector));
    
    return $detector_class;
}

sub parser {
    my $self = shift;
    
    my $grammar = q{
        startrule: intersection
                    { $item[1]; }
        | union
                    { $item[1]; }
        | single
                    { $item[1]; }
        
        parenthetical: "(" startrule ")"
                    { $item[2]; }
                    
        intersection: single "&&" startrule
                    { $return = { intersect => [$item[1], $item[3] ] }; }
        
        union: single "||" startrule
                    { $return = { union => [ $item[1], $item[3] ] }; }
        
        single: parenthetical
                    { $item[1]; }
        | "somatic" name
                    { my @name = %{$item[2]}; $return = { ('somatic ' . $name[0]) => $name[1] }; }
        | name
                    { $item[1]; }
        
        name: /[\w\.-]+/ <defer: 'This defer is an ugly but convenient counter, since it returns the number of queued events.'; >
                    { $return = { $item[1] => ($item[2] - 1) }; }
    };

    my $parser = Parse::RecDescent->new($grammar)
        or die('Failed to construct parser from grammar.');
        
    return $parser;
}

#The parser returns a data structure like { union => [{ intersection => [{ samtools => 1}, {'var-scan' => 2}]}, maq => 3] },
#which we then take to build a hash like:
# { samtools => { r453 => {snv => ['']}}, 'var-scan' => { '' => { snv => ['']} }, maq => { '0.7.1' => { snv => ['a'] } }}
sub build_detector_list {
    my $self = shift;
    my ($detector_tree, $detector_list, $detector_type, $versions, $params) = @_;
    
    #Just recursively keep looking for detectors
    my $branch_case = sub {
        my $self = shift;
        my ($combination, $subtrees, $branch_case, $leaf_case, $detector_list, $detector_type, $versions, $params) = @_;
        
        for my $subtree (@$subtrees) {
            $self->walk_tree($subtree, $branch_case, $leaf_case, $detector_list, $detector_type, $versions, $params);
        }
        
        return $detector_list;
    };
    
    #We found the detector we're looking for
    my $leaf_case = sub {
        my $self = shift;
        my ($detector_name, $index, $branch_case, $leaf_case, $detector_list, $detector_type, $versions, $params) = @_;
        
        my $version = $versions->[$index] || '';
        my $param = $params->[$index] || '';
    
        #Do not push duplicate entries
        unless(exists($detector_list->{$detector_name}{$version}{$detector_type}) and
               grep($_ eq $param, @{ $detector_list->{$detector_name}{$version}{$detector_type} })) {
        
            push @{ $detector_list->{$detector_name}{$version}{$detector_type} }, $param;
        }
        
        return $detector_list;
    };
    
    return $self->walk_tree($detector_tree, $branch_case, $leaf_case, $detector_list, $detector_type, $versions, $params);
}

sub combine_results {
    my $self = shift;
    my ($detector_tree, $detector_type, $versions, $params) = @_;
    
    #Need to combine by intersection or union
    my $branch_case = sub {
        my $self = shift;
        my ($combination, $subtrees, $branch_case, $leaf_case, $detector_type, $versions, $params) = @_;
        
        $self->error_message('Support for unions and intersections is not yet available.');
        die($self->error_message);
    };
    
    #A single detector--just find the relevant file and pass it back
    my $leaf_case = sub {
        my $self = shift;
        my ($detector_name, $index, $branch_case, $leaf_case, $detector_type, $versions, $params) = @_;
        
        #TODO When we enable the boolean expressions, the directory name will need to take into account parameters as well
        my $command_output_directory = $self->calculate_detector_output_directory($detector_name, $versions->[$index], '');
        
        #Somewhere down the line filtering should perhaps be separated from the actual detection?
        my $output_file_property = $detector_type . '_output';
        my $filtered_output_file_property = 'filtered_' . $output_file_property;
        
        my $_temp_staging_directory = $self->_temp_staging_directory;
        my $output_file = $self->$output_file_property;
        my $filtered_output_file = $self->$filtered_output_file_property;
        
        $output_file =~ s/$_temp_staging_directory/$command_output_directory/;
        $filtered_output_file =~ s/$_temp_staging_directory/$command_output_directory/;
        
        return [$output_file, $filtered_output_file];
    };
    
    return $self->walk_tree($detector_tree, $branch_case, $leaf_case, $detector_type, $versions, $params);
}
#A generic walk of the tree structure produced by the parser--takes two subs $branch_case and $leaf_case to handle the specific logic of the step
sub walk_tree {
    my $self = shift;
    my ($detector_tree, $branch_case, $leaf_case, @params) = @_;
    
    unless($detector_tree) {
        $self->error_message('No parsed detector tree provided.');
        die($self->error_message);
    }
    
    my @keys = keys %$detector_tree;
    
    #There should always be exactly one outer rule (or detector)
    unless(scalar @keys eq 1) {
        $self->error_message('Unexpected data structure encountered!  There were ' . scalar(@keys) . ' keys');
        die($self->error_message);
    }
    
    my $key = $keys[0];
    
    #Case One:  We're somewhere in the middle of the data-structure--we need to combine some results
    if($key eq 'intersect' or $key eq 'union') {
        my $value = $detector_tree->{$key};
        
        unless(ref $value eq 'ARRAY') {
            $self->error_message('Unexpected data structure encountered! I really wanted an ARRAY, not ' . (ref($value)||$value) );
            die($self->error_message);
        }
        
        return $branch_case->($self, $key, $value, $branch_case, $leaf_case, @params);
    }
    
    #Case Two: Otherwise the key should be the name of a detector, and the value is a zero-based index into the version and param arrays
    my $index = $detector_tree->{$key};
    return $leaf_case->($self, $key, $index, $branch_case, $leaf_case, @params);
}

1;
