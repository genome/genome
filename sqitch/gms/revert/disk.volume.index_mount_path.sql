-- Revert disk.volume.index_mount_path

BEGIN;

DROP INDEX disk.volume_mount_path_index;

COMMIT;
