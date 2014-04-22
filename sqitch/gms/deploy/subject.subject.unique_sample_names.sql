-- Deploy subject.subject.unique_sample_names
-- requires: subject_subject

BEGIN;

CREATE UNIQUE INDEX unique_sample_name_index ON subject.subject (name) WHERE subclass_name = 'Genome::Sample';

COMMIT;
