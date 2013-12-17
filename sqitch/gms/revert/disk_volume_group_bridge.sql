-- Revert disk_volume_group_bridge

BEGIN;

DROP TABLE IF EXISTS disk.volume_group_bridge;

COMMIT;
