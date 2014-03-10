-- Revert subject_misc_attribute

BEGIN;

DROP TABLE IF EXISTS subject.misc_attribute;

COMMIT;
