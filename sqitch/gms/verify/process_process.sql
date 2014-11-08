-- Verify process_process

BEGIN;

SELECT * FROM process.process WHERE FALSE;

ROLLBACK;
