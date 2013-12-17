-- Revert result_input

BEGIN;

DROP TABLE IF EXISTS result.input;

COMMIT;
