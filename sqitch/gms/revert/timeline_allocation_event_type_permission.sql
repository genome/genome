-- Revert timeline_allocation_event_type_permission

BEGIN;

REVOKE ALL ON TABLE timeline.allocation_event_type FROM "gms-user";

COMMIT;
