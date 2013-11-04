-- Deploy disk.allocation.index_status
-- requires: disk_allocation

BEGIN;

CREATE INDEX allocation_status_idx ON disk.allocation USING btree (status);

COMMIT;
