-- Revert disk_volume_permission

BEGIN;

REVOKE ALL ON TABLE disk.volume FROM "gms-user";

COMMIT;
