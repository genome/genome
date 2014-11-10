-- Deploy process.process.created_by_index
-- requires: process_process

BEGIN;

CREATE INDEX process_created_by_idx ON process.process USING btree (created_by);

COMMIT;
