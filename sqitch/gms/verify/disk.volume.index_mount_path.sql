-- Verify disk.volume.index_mount_path

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'volume_mount_path_index';

ROLLBACK;
