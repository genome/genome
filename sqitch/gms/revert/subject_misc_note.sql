-- Revert subject_misc_note

BEGIN;

DROP TABLE IF EXISTS subject.misc_note;

COMMIT;
