-- Deploy process.process.times_index
-- requires: process_process

BEGIN;

CREATE INDEX process_created_at_idx ON process.process
    USING btree (created_at);
CREATE INDEX process_started_at_idx ON process.process
    USING btree (started_at);
CREATE INDEX process_ended_at_idx ON process.process
    USING btree (ended_at);

COMMIT;
