package Genome::Sample::Command::Attribute::List;

class Genome::Sample::Command::Attribute::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::SubjectAttribute' 
        },
        show => {
            default_value => 'sample.id,sample.name,attribute_label,attribute_value',
        },
    ],
};

1;

