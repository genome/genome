-- Revert process.process.status_index

BEGIN;

DROP INDEX process.process_status_idx;

COMMIT;
