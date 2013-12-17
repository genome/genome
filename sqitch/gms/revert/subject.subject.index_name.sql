-- Revert subject.subject.index_name

BEGIN;

DROP INDEX subject.subject_name_index;

COMMIT;
