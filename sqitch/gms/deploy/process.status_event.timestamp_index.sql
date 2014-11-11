-- Deploy process.status_event.timestamp_index
-- requires: process_status_event

BEGIN;

CREATE INDEX process_status_event_timestamp_idx ON process.status_event
    USING btree (timestamp);

COMMIT;
