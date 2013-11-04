-- Deploy disk.allocation.index_creation_time_reallocation_time
-- requires: disk_allocation

BEGIN;

CREATE INDEX allocation_creation_reallocation_time_index
ON disk.allocation
USING btree (creation_time, reallocation_time);

COMMIT;
