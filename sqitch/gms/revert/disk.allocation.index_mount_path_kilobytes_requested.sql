-- Revert disk.allocation.index_mount_path_kilobytes_requested

BEGIN;

DROP INDEX disk.d_a_multi_mount_path_kilobytes_requested_index;

COMMIT;
