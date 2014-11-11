-- Deploy process.process.status_index
-- requires: process_process

BEGIN;

CREATE INDEX process_status_idx ON process.process USING btree (status);

COMMIT;
