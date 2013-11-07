-- Deploy timeline_allocation_event_type_permission
-- requires: timeline_allocation_event_type

BEGIN;

REVOKE ALL ON TABLE timeline.allocation_event_type FROM PUBLIC;
REVOKE ALL ON TABLE timeline.allocation_event_type FROM genome;
GRANT ALL ON TABLE timeline.allocation_event_type TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE timeline.allocation_event_type TO "gms-user";

COMMIT;
