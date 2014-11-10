-- Verify process.status_event.timestamp_index

BEGIN;

SELECT 1/count(*) FROM pg_class
    WHERE relkind = 'i' and relname = 'process_status_event_timestamp_idx';

ROLLBACK;
