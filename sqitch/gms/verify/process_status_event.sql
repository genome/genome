-- Verify process_status_event

BEGIN;

SELECT * FROM process.status_event WHERE FALSE;

ROLLBACK;
