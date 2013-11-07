-- Revert disk_group_permission

BEGIN;

REVOKE ALL ON TABLE disk."group" FROM "gms-user";

COMMIT;
