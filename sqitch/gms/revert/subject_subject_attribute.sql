-- Revert subject_subject_attribute

BEGIN;

DROP TABLE IF EXISTS subject.subject_attribute;

COMMIT;
