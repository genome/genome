-- Revert subject_role

BEGIN;

DROP TABLE IF EXISTS subject.role;

COMMIT;
