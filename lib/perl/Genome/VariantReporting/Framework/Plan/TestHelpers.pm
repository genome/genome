# This package just defines some classes so you can create plan yaml files that are valid.
package Genome::VariantReporting::Framework::Plan::TestHelpers;

use Sub::Install qw(reinstall_sub);
use Exporter 'import';

our @EXPORT_OK = qw(
    set_what_interpreter_x_requires
);

sub set_what_interpreter_x_requires {
    my @what = @_;
    reinstall_sub( {
        into => 'Genome::VariantReporting::TestInterpreter',
        as => 'requires_annotations',
        code => sub {return @what;},
    });
}

{
    package Genome::VariantReporting::TestInterpreter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::TestInterpreter {
        is => 'Genome::VariantReporting::Framework::Component::Interpreter',
        has => [
            ix_p1 => {},
            ix_p2 => {},
        ],
    };

    sub name {
        "interpreter_x";
    }

    sub requires_annotations {
        return qw(expert_one);
    }

    sub interpret_entry {
        my $self = shift;
        my $entry = shift;
        my $passed_alleles = shift;
        my %dict;
        for my $allele (@$passed_alleles) {
            my $value = $entry->info("EXP1");
                $dict{$allele} = {
                    exp1 => $value,
                };
        }
        return %dict;
    }

    sub field_descriptions {
        return (
            exp1 => 'Expert 1 info field',
        );
    }

    1;
}

{
    package Genome::VariantReporting::AnotherTestInterpreter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::AnotherTestInterpreter {
        is => 'Genome::VariantReporting::Framework::Component::Interpreter',
        has => [
            iy_p1 => {},
            iy_p2 => {},
        ],
    };

    sub requires_annotations {
        return qw(expert_one expert_two);
    }

    sub name {
        "interpreter_y";
    }

    sub interpret_entry {
        my $self = shift;
        my $entry = shift;
        my $passed_alleles = shift;
        my %dict;
        for my $allele (@$passed_alleles) {
            $dict{$allele} = {
                chrom => $entry->{chrom},
                pos => $entry->{position},
            };
        }
        return %dict;
    }

    sub field_descriptions {
        return (
            chrom => 'Chromosome',
            pos => 'Position',
        );
    }

    1;
}

{
    package Genome::VariantReporting::WithTranslationTestInterpreter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::WithTranslationTestInterpreter {
        is => 'Genome::VariantReporting::Framework::Component::Interpreter',
        has => [
            iz_p1 => { is_translated => 1 },
            iz_p2 => {},
        ],
    };

    sub name {
        "interpreter_z";
    }

    sub interpret_entry {
        my $self = shift;
        my $entry = shift;
        my $passed_alleles = shift;
        my %dict;
        for my $allele (@$passed_alleles) {
            $dict{$allele} = {
                iz_p1 => $self->iz_p1,
            };
        }
        return %dict;
    }

    sub field_descriptions {
        return (
            iz_p1 => 'The value iz_p1',
        );
    }

    1;
}

{
    package Test::TranslatedComponent;
    use Genome;
    class Test::TranslatedComponent {
        is => 'Genome::VariantReporting::Framework::Component::Interpreter',
        has => [to_translate => {is_translated => 1, is_many => 1, is => 'Text'}],
    };
}


{
    package Genome::VariantReporting::TestReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::TestReporter {
        is => 'Genome::VariantReporting::Reporter::WithHeader',
        has => [
            ra_p1 => {},
            ra_p2 => {},
        ],
    };

    sub name {
        "reporter_alpha";
    }

    sub requires_interpreters {
        return qw(interpreter_x);
    }

    sub headers {
        return qw(exp1);
    }

    sub report {
        my $self = shift;
        my $interpretations = shift;
        while (my ($allele, $interpretations_for_allele) = each %{$interpretations->{interpreter_x}}) {
            $self->_output_fh->print(_format($interpretations_for_allele->{exp1})."\n");
        }
    }

    sub _format {
        my $string = shift;
        if (defined $string ) {
            return $string;
        }
        else {
            return "-";
        }
    }

    1;
}

{
    package Genome::VariantReporting::AnotherTestReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::AnotherTestReporter {
        is => 'Genome::VariantReporting::Framework::Component::Reporter::SingleFile',
        has => [
            rb_p1 => {},
            rb_p2 => {},
        ],
    };

    sub name {
        "reporter_beta";
    }

    sub requires_interpreters {
        return qw(interpreter_x interpreter_y);
    }

    1;
}

{
    package Genome::VariantReporting::YetAnotherTestReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::YetAnotherTestReporter {
        is => 'Genome::VariantReporting::Reporter::WithHeader',
        has => [
            rc_p1 => {},
            rc_p2 => {},
        ],

    };

    sub name {
        "reporter_gamma";
    }

    sub requires_interpreters {
        return qw(interpreter_x interpreter_y);
    }

    sub headers {
        return qw(chrom pos exp1);
    }

    sub report {
        my $self = shift;
        my $interpreters = shift;
        for my $allele (keys %{$interpreters->{interpreter_y}}) {
            my $chrom = _format($interpreters->{interpreter_y}->{$allele}->{chrom});
            my $position = _format($interpreters->{interpreter_y}->{$allele}->{pos});
            my $exp1 = _format($interpreters->{interpreter_x}->{$allele}->{exp1});
            $self->_output_fh->print(join(" ", $chrom, $position, $exp1)."\n");
        }
    }

    sub _format {
        my $string = shift;
        if (defined $string ) {
            return $string;
        }
        else {
            return "*";
        }
    }

    1;
}

{
    package Genome::VariantReporting::TestDeltaReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::TestDeltaReporter {
        is => 'Genome::VariantReporting::Framework::Component::Reporter',
        has => [
            p1 => {},
        ],
    };

    sub name {
        "reporter_delta";
    }

    sub requires_interpreters {
        return qw(interpreter_x);
    }

    sub allows_hard_filters {
        return 0;
    }

    1;
}

{
    package Genome::VariantReporting::TestEpsilonReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::TestEpsilonReporter {
        is => 'Genome::VariantReporting::Framework::Component::Reporter',
        has => [
        ],
    };

    sub name {
        "reporter_epsilon";
    }

    sub requires_interpreters {
        return qw(interpreter_z);
    }

    sub allows_hard_filters {
        return 0;
    }

    1;
}


{
    package Genome::VariantReporting::TestExpert;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::TestExpert {
        is => 'Genome::VariantReporting::Framework::Component::Expert',
        has => [
        ],
    };

    sub name {
        "expert_one";
    }

    1;
}

{
    package Genome::VariantReporting::ExpertOneAdaptor;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::ExpertOneAdaptor {
        is => 'Genome::VariantReporting::Framework::Component::Adaptor',
        has_planned_output => [
            e1_p1 => {},
            e1_p2 => {},
            e1_translated1 => {is_translated => 1},
            e1_translated2 => {is_translated => 1},
        ],
    };

    sub name {
        "expert_one";
    }

    1;
}

{
    package Genome::VariantReporting::ExpertOneRun;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::ExpertOneRun {
        is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    };

    sub name {
        "expert_one";
    }

    1;
}

{
    package Genome::VariantReporting::AnotherTestExpert;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::AnotherTestExpert {
        is => 'Genome::VariantReporting::Framework::Component::Expert',
        has => [
        ],
    };

    sub name {
        "expert_two";
    }

    1;
}

{
    package Genome::VariantReporting::ExpertTwoAdaptor;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::ExpertTwoAdaptor {
        is => 'Genome::VariantReporting::Framework::Component::Adaptor',
        has_planned_output => [
            e2_p1 => {},
            e2_p2 => {},
            e2_translated1 => {is_translated => 1},
            e2_translated2 => {is_translated => 1},
        ],
    };

    sub name {
        "expert_two";
    }

    1;
}

{
    package Genome::VariantReporting::ExpertTwoRun;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::ExpertTwoRun {
        is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    };

    sub name {
        "expert_two";
    }

    1;
}


{
    package Genome::VariantReporting::TestFilter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::TestFilter {
        is => 'Genome::VariantReporting::Framework::Component::Filter',
        has => [
            f1_p1 => {},
            f1_p2 => {},
        ],
    };

    sub name {
        'filter_one';
    }

    sub filter_entry {
        my $self = shift;
        my $entry = shift;
        my %returns;
        for my $allele (@{$entry->{alternate_alleles}}) {
            if (length $allele >= $self->f1_p1) {
                $returns{$allele} = 0;
            }
            else {
                $returns{$allele} = 1;
            }
        }
        return %returns;
    }

    sub vcf_id {
        return "FIRST_FILTER";
    }

    sub vcf_description  {
        return "The first filter";
    }

    1;
}

{
    package Genome::VariantReporting::AnotherTestFilter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::AnotherTestFilter {
        is => 'Genome::VariantReporting::Framework::Component::Filter',
        has => [
            f2_p1 => {},
            f2_p2 => {},
        ],
    };

    sub requires_annotations {
        return qw(expert_two);
    }

    sub name {
        'filter_two';
    }

    sub filter_entry {
        my $self = shift;
        my $entry = shift;
        return map{$_ => 1} @{$entry->{alternate_alleles}};
    }

    sub vcf_id {
        return "OTHER_FILTER";
    }

    sub vcf_description  {
        return "The other filter";
    }

    1;
}

# These allow the above classes to be used to create DAGs
# that are used to generate xml for testing DAG generation.
$INC{'Genome/VariantReporting/ExpertOneAdaptor.pm'} = '1',
$INC{'Genome/VariantReporting/ExpertTwoAdaptor.pm'} = '1',
$INC{'Genome/VariantReporting/ExpertOneRun.pm'} = '1',
$INC{'Genome/VariantReporting/ExpertTwoRun.pm'} = '1',

1;
