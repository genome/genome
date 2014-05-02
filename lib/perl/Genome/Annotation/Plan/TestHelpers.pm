# This package just defines some classes so you can create plan yaml files that are valid.
package Genome::Annotation::Plan::TestHelpers;

use Sub::Install qw(reinstall_sub);
use Exporter 'import';

our @EXPORT_OK = qw(
    set_what_interpreter_x_requires
);

sub set_what_interpreter_x_requires {
    my @what = @_;
    reinstall_sub( {
        into => 'Genome::Annotation::TestInterpreter',
        as => 'requires_experts',
        code => sub {return @what;},
    });
}

{
    package Genome::Annotation::TestInterpreter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestInterpreter {
        is => 'Genome::Annotation::Interpreter::Base',
        has => [
            ix_p1 => {},
            ix_p2 => {},
        ],
    };

    sub name {
        "interpreter_x";
    }

    sub requires_experts {
        return qw(expert_one);
    }

    sub available_fields {
        return qw(exp1);
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

    1;
}

{
    package Genome::Annotation::AnotherTestInterpreter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::AnotherTestInterpreter {
        is => 'Genome::Annotation::Interpreter::Base',
        has => [
            iy_p1 => {},
            iy_p2 => {},
        ],
    };

    sub requires_experts {
        return qw(expert_one expert_two);
    }

    sub name {
        "interpreter_y";
    }

    sub available_fields {
        return qw(chrom pos);
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

    1;
}

{
    package Genome::Annotation::TestReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestReporter {
        is => 'Genome::Annotation::Reporter::WithHeader',
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
        for my $allele (keys %{$interpretations->{interpreter_x}}) {
            $self->_output_fh->print(_format($interpretations->{interpreter_x}->{$allele}->{exp1})."\n");
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
    package Genome::Annotation::AnotherTestReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::AnotherTestReporter {
        is => 'Genome::Annotation::Reporter::Base',
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
    package Genome::Annotation::YetAnotherTestReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::YetAnotherTestReporter {
        is => 'Genome::Annotation::Reporter::WithHeader',
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
    package Genome::Annotation::TestExpert;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestExpert {
        is => 'Genome::Annotation::Expert::Base',
        has => [
            e1_p1 => {},
            e1_p2 => {},
        ],
    };

    sub name {
        "expert_one";
    }

    1;
}

{
    package Genome::Annotation::ExpertOneAdaptor;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::ExpertOneAdaptor {
        is => 'Genome::Annotation::AdaptorBase',
        has_planned_output => [
            e1_p1 => {},
            e1_p2 => {},
        ],
    };

    sub name {
        "expert_one";
    }

    1;
}

{
    package Genome::Annotation::ExpertOneRun;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::ExpertOneRun {
        is => 'Genome::Annotation::Expert::CommandBase',
    };

    sub name {
        "expert_one";
    }

    1;
}

{
    package Genome::Annotation::AnotherTestExpert;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::AnotherTestExpert {
        is => 'Genome::Annotation::Expert::Base',
        has => [
            e2_p1 => {},
            e2_p2 => {},
        ],
    };

    sub name {
        "expert_two";
    }

    1;
}

{
    package Genome::Annotation::ExpertTwoAdaptor;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::ExpertTwoAdaptor {
        is => 'Genome::Annotation::AdaptorBase',
        has_planned_output => [
            e2_p1 => {},
            e2_p2 => {},
        ],
    };

    sub name {
        "expert_two";
    }

    1;
}

{
    package Genome::Annotation::ExpertTwoRun;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::ExpertTwoRun {
        is => 'Genome::Annotation::Expert::CommandBase',
    };

    sub name {
        "expert_two";
    }

    1;
}


{
    package Genome::Annotation::TestFilter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestFilter {
        is => 'Genome::Annotation::Filter::Base',
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

    1;
}

{
    package Genome::Annotation::AnotherTestFilter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::AnotherTestFilter {
        is => 'Genome::Annotation::Filter::Base',
        has => [
            f2_p1 => {},
            f2_p2 => {},
        ],
    };

    sub requires_experts {
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

    1;
}

# These allow the above classes to be used to create DAGs
# that are used to generate xml for testing DAG generation.
$INC{'Genome/Annotation/ExpertOneAdaptor.pm'} = '1',
$INC{'Genome/Annotation/ExpertTwoAdaptor.pm'} = '1',
$INC{'Genome/Annotation/ExpertOneRun.pm'} = '1',
$INC{'Genome/Annotation/ExpertTwoRun.pm'} = '1',

1;
