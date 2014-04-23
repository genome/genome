# This package just defines some classes so you can create plan yaml files that are valid.
package Genome::Annotation::Plan::TestHelpers;

{
    package Genome::Annotation::TestInterpreter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestInterpreter {
        is => 'Genome::Annotation::InterpreterBase',
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

    1;
}

{
    package Genome::Annotation::AnotherTestInterpreter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::AnotherTestInterpreter {
        is => 'Genome::Annotation::InterpreterBase',
        has => [
            ix_p1 => {},
            ix_p2 => {},
        ],
    };

    sub name {
        "interpreter_y";
    }

    1;
}

{
    package Genome::Annotation::TestReporter;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestReporter {
        is => 'Genome::Annotation::ReporterBase',
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

    1;
}

{
    package Genome::Annotation::TestExpert;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestExpert {
        is => 'Genome::Annotation::ExpertBase',
        has => [
            e1_p1 => {},
            e1_p2 => {},
        ],
    };

    sub name {
        "expert_one";
    }
}

{
    package Genome::Annotation::AnotherTestExpert;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::AnotherTestExpert {
        is => 'Genome::Annotation::ExpertBase',
        has => [
            e1_p1 => {},
            e1_p2 => {},
        ],
    };

    sub name {
        "expert_two";
    }

    1;
}

1;
