-- Revert disk.allocation.index_mount_path

BEGIN;

DROP INDEX disk.d_a_mount_path_index;

COMMIT;
