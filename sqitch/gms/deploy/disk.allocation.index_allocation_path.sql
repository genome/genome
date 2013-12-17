-- Deploy disk.allocation.index_allocation_path
-- requires: disk_allocation

BEGIN;

CREATE INDEX allocation_allocation_path_index ON disk.allocation USING btree (allocation_path);

COMMIT;
