-- Verify disk_volume_group_bridge

BEGIN;

SELECT volume_id, group_id
FROM disk.volume_group_bridge
WHERE FALSE;

ROLLBACK;
