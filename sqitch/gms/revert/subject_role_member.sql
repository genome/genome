-- Revert subject_role_member

BEGIN;

DROP TABLE IF EXISTS subject.role_member;

COMMIT;
