-- Revert subject.subject.unique_sample_names

BEGIN;

DROP INDEX subject.unique_sample_name_index;

COMMIT;
