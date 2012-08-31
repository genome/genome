package Genome::Model::Tools::SeeFive::Trial;

use strict;
use warnings;
use Genome;            
use Genome::Model::Tools::SeeFive::Rule;

class Genome::Model::Tools::SeeFive::Trial {
    is => 'Command',
    has => [ 
        rules           => { is => 'Genome::Model::Tools::SeeFive::Rule', reverse_id_by => 'trial', is_many => 1 },
        n               => { },
    ],
};

sub generate_callback_for_headers {
    my $self = shift;
    my @headers = @_;
    my $src = $self->perl_source_for_headers(@headers);
    print $src;
    my $callback = eval $src;
    unless ($callback) {
        die $@;
    }
    return $callback;
}

sub perl_source_for_headers {
    my $self = shift;
    my @headers = @_;
    my $src;
    $src .= "sub {\n";
    $src .= "    my \$data = \$_[0];\n";
    $src .= "    my (" 
            . join(",", map { '$' . $_ } @headers ) 
            . ') = @$data{@headers};'
            . "\n";

    my @rules = $self->rules;
    unless (@rules) {
        $DB::single = $DB::stopper;
        die "No rules for trial?";
    }
    my $default_value = 'G';
    $src .= '    my %answers;' . "\n";
    $src .= '    my ($answer,$prob);' . "\n";
    for my $rule (@rules) {
        my $rsrc= $rule->perl_src;
        $rsrc =~ s/\n/  /g;
        $src .= '    ($answer,$prob) = ' . $rsrc . ";\n";
        $src .= '    if(defined $answer) { $answers{$answer}+=$prob; }' . "\n\n";        
    }
    $src .= '    my $best_prob = 0; my $best_answer = "";' . "\n";
    $src .= '    for my $answer (keys %answers) { if ($answers{$answer} > $best_prob) { $best_answer = $answer; $best_prob = $answers{$answer} } }' . "\n";
    $src .= '    if (not values %answers) { $best_answer = "' . $default_value . '"' . " }\n"; 
    $src .= '    my @debug_notes = %answers;' . "\n";
    $src .= '    return ($best_answer,$best_prob,"@debug_notes");' . "\n";
    $src .= "};\n";
    
    return $src;
}

sub create_list_from_c5src {
    my $self  = shift;
    my $c5src = shift;
    my @lines = split(/\n/,$c5src);
    chomp @lines;

    my @trials;
    my $current_trial;
    my $current_rule;
    my $current_rule_lines;

    for (@lines) { 
        if (/Trial\s+(\d+)/) {
            $current_trial = Genome::Model::Tools::SeeFive::Trial->create(
                n => $1,
            );
            die unless $current_trial;
            push @trials, $current_trial;
        }
        next unless defined $current_trial;
        next if /^Rules:\s*$/;
        
        if (
            my ($trial,$n,$possible_items,$wrong_items,$lift) 
            = ($_ =~ qr{^Rule ([\d\.]+)/([\d\.]+): \(([\d\.]+)\/([\d\.]+), lift ([\d\.]+)\)\s*$}) 
        ) {
            $current_rule = Genome::Model::Tools::SeeFive::Rule->create(
                trial_id => $current_trial->id,
                n => $n,
                possible_items => $possible_items,
                wrong_items => $wrong_items,
                lift => $lift,
            );
            $current_rule_lines = [];
            $current_rule->lines($current_rule_lines);
            next;
        }
        next unless defined $current_rule;

        if (/^\s*$/) {
            $current_rule = undef;
        }
        else {
            push @$current_rule_lines, $_;
        }
    }
    
    return @trials;
}

1;

