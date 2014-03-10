-- Revert timeline_base_permission

BEGIN;

REVOKE ALL ON TABLE timeline.base FROM "gms-user";

COMMIT;
