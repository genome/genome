-- Revert disk_volume_group_bridge_permission

BEGIN;

REVOKE ALL ON TABLE disk.volume_group_bridge FROM "gms-user";

COMMIT;
