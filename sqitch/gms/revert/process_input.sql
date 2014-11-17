-- Revert process_input

BEGIN;

DROP TABLE process.input;

COMMIT;
