-- Revert result_user

BEGIN;

DROP TABLE IF EXISTS result."user";

COMMIT;
