-- Deploy disk.allocation.index_mount_path_kilobytes_requested
-- requires: disk_allocation

BEGIN;

CREATE INDEX d_a_multi_mount_path_kilobytes_requested_index
ON disk.allocation
USING btree (mount_path, kilobytes_requested);

COMMIT;
