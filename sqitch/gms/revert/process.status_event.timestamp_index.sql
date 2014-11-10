-- Revert process.status_event.timestamp_index

BEGIN;

DROP INDEX process.process_status_event_timestamp_idx;

COMMIT;
