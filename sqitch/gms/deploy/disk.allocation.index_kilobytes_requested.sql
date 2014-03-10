-- Deploy disk.allocation.kilobytes_requested
-- requires: disk_allocation

BEGIN;

CREATE INDEX kilobytes_requested_index on disk.allocation using btree (kilobytes_requested);

COMMIT;
