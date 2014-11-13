-- Verify process_input

BEGIN;

SELECT * FROM process.input WHERE FALSE;

ROLLBACK;
