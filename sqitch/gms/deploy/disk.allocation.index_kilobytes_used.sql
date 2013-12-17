-- Deploy disk.allocation.kilobytes_used
-- requires: disk_allocation

BEGIN;

CREATE INDEX kilobytes_used_index on disk.allocation using btree (kilobytes_used);

COMMIT;
