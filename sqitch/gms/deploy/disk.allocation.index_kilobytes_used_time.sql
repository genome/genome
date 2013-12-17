-- Deploy disk.allocation.index_kilobytes_used_time
-- requires: disk_allocation

BEGIN;

CREATE INDEX d_a_kilobytes_used_time_index ON disk.allocation USING btree (kilobytes_used_time);

COMMIT;
