-- Revert timeline_allocation_permission

BEGIN;

REVOKE ALL ON TABLE timeline.allocation FROM "gms-user";

COMMIT;
