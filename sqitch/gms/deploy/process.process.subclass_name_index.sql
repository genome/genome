-- Deploy process.process.subclass_name_index
-- requires: process_process

BEGIN;

CREATE INDEX process_subclass_name_idx ON process.process
    USING btree (subclass_name);

COMMIT;
