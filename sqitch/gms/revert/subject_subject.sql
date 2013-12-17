-- Revert subject_subject

BEGIN;

DROP TABLE IF EXISTS subject.subject;

COMMIT;
