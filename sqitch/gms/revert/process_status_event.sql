-- Revert process_status_event

BEGIN;

DROP TABLE process.status_event;

COMMIT;
