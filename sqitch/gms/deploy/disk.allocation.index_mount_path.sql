-- Deploy disk.allocation.index_mount_path
-- requires: disk_allocation

BEGIN;

CREATE INDEX d_a_mount_path_index ON disk.allocation USING btree (mount_path);

COMMIT;
