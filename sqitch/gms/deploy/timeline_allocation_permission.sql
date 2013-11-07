-- Deploy timeline_allocation_permission
-- requires: timeline_allocation

BEGIN;

REVOKE ALL ON TABLE timeline.allocation FROM PUBLIC;
REVOKE ALL ON TABLE timeline.allocation FROM genome;
GRANT ALL ON TABLE timeline.allocation TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE timeline.allocation TO "gms-user";

COMMIT;
