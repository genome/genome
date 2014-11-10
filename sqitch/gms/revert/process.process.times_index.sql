-- Revert process.process.times_index

BEGIN;

DROP INDEX process.process_created_at_idx;
DROP INDEX process.process_started_at_idx;
DROP INDEX process.process_ended_at_idx;

COMMIT;
