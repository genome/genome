-- Revert subject_user

BEGIN;

DROP TABLE IF EXISTS subject."user";

COMMIT;
