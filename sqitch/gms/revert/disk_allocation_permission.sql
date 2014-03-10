-- Revert disk_allocation_permission

BEGIN;

REVOKE ALL ON TABLE disk.allocation FROM "gms-user";

COMMIT;
